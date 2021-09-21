## Author : SUPRATIK BASU

sparse.km.fr.2 <- function(X, C, s, gt=NULL, tolerance = 1e-04)
{
  N <- nrow(X)
  d <- ncol(X)
  centroid <- matrix(0, C, d)
  p <- pps1(rep(1/N, N))
  centroid[1,] <- as.numeric(X[p,])
  for(i in 1:C-1)
  {
    dist.mat <- rep(0, N)
    for(j in 1:N)
    {
      #temp <- apply(centroid[1:i,], 1, function(vec) K(X[j,], vec, sigma))
      dist.mat[j] <- sum((X[j, ] - centroid[i, ])^2)
      
    }
    centroid[i+1,] <- X[pps1(dist.mat),]
  }
  Z <- matrix(0, N, C)
  obj.old <- 0
  for(i in 1:N)
  {
    dist.mat <- rep(0, C)
    for(j in 1:C)
    {
      dist.mat[j] <- sum((X[i, ] - centroid[j, ]) ^ 2)
    }
    Z[i, which.min(dist.mat)] <- 1
    obj.old <- obj.old + min(dist.mat)
  }
  
  t <- 0
  while(TRUE)
  {
    D <- matrix(0, C, d)
    new.Z <- matrix(0, N, C)
    obj.new <- 0
    r <- matrix(0, C, d)
    for(j in 1:C)
    {
      centroid[j, ] <- colMeans(X[which(Z[ ,j]==1), ])
      if(sum(Z[,j])==0)
      {
        centroid[j,] <- rep(0, d)
      }
      for(l in 1:d)
      {
        D[j,l] <- sum(Z[ ,j]) * centroid[j,l] ^ 2
      }
      r[j, ] <- d + 1 - rank(D[j, ])
    }
    for(i in 1:N)
    {
      dist.mat <- rep(0, C)
      for(j in 1:C)
      {
        u <- which(r[j, ] <= s)
        v <- which(r[j, ] > s)
        dist.mat[j] <- sum((X[i, r[j,u]] - centroid[j, r[j,u]]) ^ 2)
        if(s<d)
        {
          dist.mat[j] <- dist.mat[j] + sum((X[i, r[j,v]]) ^ 2)
        }
      }
      new.Z[i, which.min(dist.mat)] <- 1
      obj.new <- obj.new + min(dist.mat)
    }
    Z <- new.Z
    if(abs(obj.new - obj.old) < tolerance)
    {
      break
    }
    obj.old <- obj.new
    t <- t+1
  }
  asg.vector <- rep(0, N)
  for(i in 1:C)
    asg.vector <- asg.vector + i*Z[ ,i]
  
  print(NMI(asg.vector, gt))
  return(list("Z" = Z, "vec" = asg.vector, "centroids" = centroid, "NMI" = NMI(asg.vector, gt)))
}
