kernel.km <- function(X, C, sigma, initial, tolerance=1e-07)
{
  N <- nrow(X)
  d <- ncol(X)
  centroid <- initial
  Z <- matrix(0, N, C)
  kernel <- matrix(0, N, N)
  for(i in 1:N)
  {
    for(j in 1:N){
    kernel[i,j] <- sigma*sqrt(2*pi)* dmvnorm(X[i,],X[j,],sigma^2 * diag(ncol(X)))}
  
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
    print(t)
  }
  asg.vector <- rep(0, N)
  for(i in 1:C)
    asg.vector <- asg.vector + i*Z[,i]
  if(d==2)
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
  return(list("Z"=Z, "vec"=asg.vector, "count"=t))
}