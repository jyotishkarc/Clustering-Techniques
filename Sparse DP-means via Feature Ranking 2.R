sparse.dp.fr <- function(X, s, lambda, tolerance = 1e-07)
{
  N <- nrow(X)
  d <- ncol(X)
  centroid <- colMeans(X)
  Z <- rep(1, N)
  C <- 1
  t <- 0
  obj.old <- 0
  for(i in 1:N)
  {
    obj.old <- obj.old + sum((X[i,] - centroid)^2)
  }
  centroid <- matrix(centroid, 1, d)
  while(TRUE)
  {
    D <- matrix(0, C, d)
    rank <- matrix(0, C, s)
    obj.new <- lambda * C
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
      if(min(dist.mat) > lambda)
      {
        C <- C+1
        Z[i] <- C
        centroid <- matrix(c(as.numeric(t(centroid)), X[i,]), nrow = C, byrow = TRUE)
        obj.new <- obj.new + lambda
      }
      else
      {
        Z[i] <- which.min(dist.mat)
        obj.new <- obj.new + min(dist.mat)
      }
    }
    
    if(abs(obj.new - obj.old) < tolerance)
    {
      break
    }
    obj.old <- obj.new
    t <- t+1
  }
}