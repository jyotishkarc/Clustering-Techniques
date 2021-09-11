## Author : SUPRATIK BASU

sparse.dpm.fr.2 <- function(X, s, lambda, gt = NULL, tolerance = 1e-03)
{
  N <- nrow(X)
  d <- ncol(X)
  centroid <- colMeans(X)
  Z <- rep(1, N)
  C <- 1
  t <- 0
  obj.old <- lambda
  for(i in 1:N)
  {
    obj.old <- obj.old + sum((X[i,] - centroid)^2)
  }
  centroid <- matrix(centroid, 1, d)
  while(t<=10)
  {
    D <- matrix(0, C, d)
    r <- matrix(0, C, d)
    centroid <- matrix(0,C,d)
    obj.new <- lambda * C
    for(j in 1:C)
    {
      if(length(which(Z==j))==0)
      {
        centroid[j,] <- rep(0, d)
      }
      else if(length(which(Z==j))==1)
      {
        centroid[j,] <- X[which(Z==j), ]
      }
      else
      {
        # print(length(colMeans(X[which(Z==j),])))
        centroid[j, ] <- colMeans(X[which(Z==j), ])
      }
      for(l in 1:d)
      {
        D[j,l] <- centroid[j,l] ^ 2
      }
      
      r[j, ] <- rank(D[j, ])
    }
    for(i in 1:N)
    {
      dist.mat <- rep(0, C)
      for(j in 1:C)
      {
        u <- which(r[j, ] <= s)
        if(s<d){
          v <- which(r[j, ] > s)}
        dist.mat[j] <- sum((X[i, r[j,u]] - centroid[j, r[j,u]]) ^ 2) 
        if(s < d)
        {
          dist.mat[j] <- dist.mat[j] + sum((X[i, r[j,v]]) ^ 2)
        }
      }
      if(min(dist.mat) > lambda)
      {
        C <- C+1
        Z[i] <- C
        centroid <- matrix(c(as.vector(t(centroid)), X[i,]), nrow = C, byrow = TRUE)
        D <- matrix(c(as.vector(D),centroid[C,]^2),C,d,byrow=T)
        r <- matrix(c(as.vector(r),rank(X[i,]^2)),C,d,byrow=T)
        obj.new <- obj.new + lambda
      }
      else
      {
        Z[i] <- which.min(dist.mat)
        obj.new <- obj.new + min(dist.mat)
      }
    }
    v <- sort(unique(Z))
    
    u <- Z
    for(i in 1:(length(v)))
    {
      u[which(Z==v[i])]=i
    }
    Z <- u
    C <- max(Z)
    if(abs(obj.new / obj.old - 1) < tolerance)
    {
      break
    }
    #print(obj.new)
    obj.old <- obj.new
    print(C)
    t <- t+1
    if(t>500)
    {
      break
    }
  }
  print(C)
  if(is.null(gt)==F)
  {
    print(aricode::NMI(Z, gt))
  }
  return(list("C"=C,"Z"=Z,"NMI"=aricode::NMI(Z,gt)))
}
