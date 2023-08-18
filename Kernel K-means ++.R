
kernel.km.plus.plus <- function(X, C, sigma, ground = NULL, tolerance = 1e-07)
{
  X <- as.matrix(X)
  #### Declarations
  N <- nrow(X)
  d <- ncol(X)
  Z <- matrix(0, N, C)
  centroid <- matrix(0, C, d)
  
  #### Centroid Initialization
  
  # print(X[pps1(rep(1/N, N)),])
  centroid[1,] <- X[pps1(rep(1/N, N)),]
  for(i in 1:C-1)
  {
    dist.mat <- rep(0, N)
    for(j in 1:N)
    {
      #temp <- apply(centroid[1:i,], 1, function(vec) K(X[j,], vec, sigma))
      temp <- c()
      for(l in 1:i)
      {
        temp[l] <- K(X[j,], centroid[l,], sigma)
      }
      dist.mat[j] <- 2 - min(temp)
    }
    centroid[i+1,] <- X[pps1(dist.mat),]
  }
  
  kernel <- matrix(0, N, N)
  for(i in 1:N)
  {
    for(j in 1:N){
      kernel[i,j] <- sigma*sqrt(2*pi)* dmvnorm(X[i,],X[j,],sigma^2 * diag(ncol(X)))}
    
  }
  
  print("HI")
  
  #### Initial Cluster Assignments
  for(i in 1:N)
  {
    dist.mat <- rep(sigma*sqrt(2*pi), C)
    for(j in 1:C)
    {
      dist.mat[j] <- dist.mat[j]*(2 - 2*dmvnorm(X[i,], centroid[j,], sigma^2 * diag(d)))
    }
    Z[i, which.min(dist.mat)] <- 1
  }
  
  #### Initial Value of Objective Function
  obj.old <- 0
  for(i in 1:N)
  {
    for(j in 1:C)
    {
      obj.old <- obj.old + Z[i,j]*(2-2*K(X[i,], centroid[j,], sigma)) 
    }
  }
  t <- 1
  while(t<=500)
  {
    #### Re-Evaluation of Cluster Assignments
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
  print(t)
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
  if(is.null(ground) == F){
  ground <- as.numeric(ground)
  
  ARI.clus <- aricode::ARI(asg.vector, ground)
  NMI.clus <- aricode::NMI(asg.vector, ground)
  
  print(ARI.clus)
  print(NMI.clus)
  
  return(list("Z"=Z, "vec"=asg.vector, "count"=t, "ARI"=ARI.clus, "NMI"=NMI.clus))}
  return(list("Z"=Z, "vec"=asg.vector, "count"=t))
}

K <- function(u, v, sigma)
{
  return(exp(-((u-v)^2)/(2*sigma^2)))
}
