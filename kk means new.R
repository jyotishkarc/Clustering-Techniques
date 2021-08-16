kernel.km <- function(X, C, sigma, initial, ground = NULL, tolerance=1e-07)
{
  X <- as.matrix(X)
  N <- nrow(X)
  d <- ncol(X)
  centroid <- as.matrix(initial)
  Z <- matrix(0, N, C)
  kernel <- matrix(0, N, N)
  for(i in 1:N)
  {
    for(j in 1:N){
    kernel[i,j] <- exp(-sum((X[i,]-X[j,])^2)/(2*sigma^2))}
   }
  
  print("HI")
  
  #Initial Cluster Assignments
  for(i in 1:N)
  {
    dist.mat <- rep(sigma*sqrt(2*pi), C)
    for(j in 1:C)
    {
      dist.mat[j] <- dist.mat[j]*(2 - 2*dmvnorm(X[i,], centroid[j,], sigma^2 * diag(d)))
    }
    Z[i, which.min(dist.mat)] <- 1
  }
  
  #Initial Value of Objective Function
  obj.old <- 0
  for(i in 1:N)
  {
    for(j in 1:C)
    {
      obj.old <- obj.old + Z[i,j]*(2-2*dmvnorm(X[i,], centroid[j,], sigma^2 * diag(d))) 
    }
  }
  t <- 1
  while(t<=500)
  {
    #Re-Evaluation of Cluster Assignments
    new.Z <- matrix(0, N, C)
    obj.new <- 0
    for(i in 1:N)
    {
      dist.mat <- rep (sigma*sqrt(2*pi), C)
      v <- c()
      for(j in 1:C)
      {
        v <- which(Z[,j]==1)
        dist.mat[j] <- dist.mat[j]*(1-2*sum(kernel[i,v])/sum(Z[,j])+sum(kernel[v,v])/(sum(Z[,j]))^2)
      }
      new.Z[i, which.min(dist.mat)] <- 1
      obj.new <- obj.new+which.min(dist.mat)
    }
    
    #obj.new <- sum(diag(kernel)) - sum(diag(kernel%*%new.Z%*%solve(t(new.Z)%*%new.Z)%*%t(new.Z)))
    if(abs(obj.new-obj.old)<tolerance)
      break
    Z <- new.Z
    obj.old <- obj.new
    t <- t+1
  }
  
  asg.vector <- rep(0, N)
  for(i in 1:C)
    asg.vector <- asg.vector + i*Z[,i]
  
  if(FALSE)
  {
  X=as.data.frame(X)
  original.plot <- ggplot(X , aes(V1,V2)) +  
    geom_point(size = 2) + theme(legend.position = "none") 
  
  print(original.plot) 
  
  fg <- c() 
  for (k in 1:nrow(X)) { 
    fg[k] <- grDevices::rainbow(length(unique(asg.vector)))[asg.vector[k]] 
  } 
  
  original.plot <- ggplot(X , aes(V1,V2)) +  
    geom_point(size = 2) + theme(legend.position = "none") 
  
  dpm.plot <- ggplot(X , aes(V1,V2 , color = fg)) + 
    geom_point(size = 2) + theme(legend.position = "none") 
  
  print(original.plot) 
  print(dpm.plot)
  }
  
  ground <- as.numeric(ground)
  
  ARI.clus <- aricode::ARI(asg.vector, ground)
  NMI.clus <- aricode::NMI(asg.vector, ground)
  
  print(ARI.clus)
  print(NMI.clus)
  
  return(list("Z"=Z, "vec"=asg.vector, "count"=t, "ARI"=ARI.clus, "NMI"=NMI.clus))
  
  # return(list("Z"=Z, "vec"=asg.vector, "count"=t))
}

# bisection <- function(fun, left, right)
# {
#   t <- 0
#   while(t <= 500)
#   {
#     med <- (left + right) / 2
#     if(fun(med) == 0)
#       return(med)
#     else if (fun(med) > 0)
#       right <- med
#     else 
#       left <- med
#     t <- t+1
#   }
#   return(med)
# }


# centroid.init <- function(data, gt, C, sigma)
# {
#   N <- nrow(data)
#   d <- ncol(data)
#   centroid <- matrix(0, C, d)
#   for(i in 1:C)
#   {
#     vec <- which(gt == i)
#     mat <- matrix(0, length(vec), d)
#     for(j in 1:length(vec))
#     {
#       mat[j,] <- data[vec[j],] / (sigma * exp(sum((data[vec[j],])^2)/(2 * sigma^2)))
#     }
#     temp <- colMeans(mat)
#     print(temp)
#     phi.norm <- uniroot(function(y) y / exp(y^2 / 2) - sqrt(sum(temp^2)), c(1,10000))
#     phi.norm <- phi.norm$root
#     centroid[i,] <- temp * sigma * exp(phi.norm ^2 /2)
#   }
#   return(centroid)
# }
