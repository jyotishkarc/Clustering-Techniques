## Author : SUPRATIK BASU

sparse.km.fr.2 <- function(X, C, s, initial, tolerance = 1e-07)
{
  N <- nrow(X)
  d <- ncol(X)
  centroid <- initial
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
    rank <- matrix(0, C, s)
    for(j in 1:C)
    {
      centroid[j, ] <- colMeans(X[which(Z[ ,j]==1), ])
      for(l in 1:d)
      {
        D[j,l] <- sum(Z[ ,j]) * centroid[j,l] ^ 2
      }
      rank[j, ] <- order(D[j, ])
    }
    for(i in 1:N)
    {
      dist.mat <- rep(0, C)
      for(j in 1:C)
      {
        dist.mat[j] <- sum((X[i, rank[j,1:s]] - centroid[j, rank[j,l]]) ^ 2) + sum((X[i, rank[j,s+1:N]]) ^ 2)
      }
      new.Z[i, which.min(dist.mat)] <- 1
      obj.new <- obj.new + min(dist.mat)
    }
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
  return(list("Z" = Z, "vec" = asg.vector, "centroids" = centroid))
}
