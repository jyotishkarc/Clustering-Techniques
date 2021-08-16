kernel.DP.means <- function(X, lambda, sigma, ground = NULL, tolerance=1e-03)
{
  N <- nrow(X)
  d <- ncol(X)
  Z <- rep(1,N)
  C <- 1
  
  kernel <- matrix(0, N, N)
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      kernel[i,j] <- K(X[i,], X[j,], sigma)
    }
  }
  
  obj.old <- lambda
  for(i in 1:N)
  {
    for(j in 1:C)
    {
      obj.old <- obj.old + 1 - 2 * sum(kernel[i,])/N + sum(kernel)/N^2
    }
  }
  
  t <- 0
  while(t<=500)
  {
    for(i in 1:N)
    {
      dist.mat <- rep(0, C)
      for( j in 1:C)
      {
        v <- which(Z==j)
        dist.mat[j] <- 1 - 2 * sum(kernel[i,v])/(length(v)) + 
          sum(kernel[v,v])/((length(v))^2)
      }
      
      if(min(dist.mat) >= lambda)
      {
        C <- C+1
        Z[i] <- C
      }
      else
      {
        Z[i] <- which.min(dist.mat)
      }
    }
    
    obj.new <- lambda*C
    
    for(i in 1:N)
    {
      v <- which(Z==Z[i])
      obj.new <- obj.new + 1 - 2 * sum(kernel[i,v])/length(v) + 
        sum(kernel[v,v])/(length(v))^2
    }
    
    if(abs(obj.old-obj.new)<tolerance) {break}
    obj.old <- obj.new
    t <- t+1
    print(t)
  }
  if(d==2)
  {
    X=as.data.frame(X)
    original.plot <- ggplot(X , aes(V1,V2)) +  
      geom_point(size = 2) + theme(legend.position = "none") 
    
    print(original.plot) 
    
    fg <- c() 
    for (k in 1:nrow(X)) { 
      fg[k] <- grDevices::rainbow(length(unique(Z)))[Z[k]] 
    } 
    
    original.plot <- ggplot(X , aes(V1,V2)) +  
      geom_point(size = 2) + theme(legend.position = "none") 
    
    dpm.plot <- ggplot(X , aes(V1,V2 , color = fg)) + 
      geom_point(size = 2) + theme(legend.position = "none") 
    
    print(original.plot) 
    print(dpm.plot)
  }
  
  ari.clus <- aricode::ARI(Z, ground)
  nmi.clus <- aricode::NMI(Z, ground)
  
  print(ari.clus)
  print(nmi.clus)
  
  return(list("Z" = Z, "C" = C))
}

K <- function(u,v,sigma)
{
  return(exp(-sum((u-v)^2)/(2*sigma^2)))
}
