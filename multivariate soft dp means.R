dp.soft.multi <- function(X, lambda, tolerance = 1e-03)
{
  N <- nrow(X)
  d <- ncol(X)
  mu <- t(as.matrix(colMeans(X)))
  C <- 1
  prob <- 1
  sigma <- (1-1/N) * cov(X)
  sigma <- list(sigma)
  obj.old <- 100
  
  while(TRUE)
  {
    for(i in 1:N)
    {
      dist.mat <- rep(0, N)
      for(j in 1:C)
      {
        print(sigma)
        dist.mat[i] <- dist.mat[i] + prob[j] * exp(-(X[i,] - mu[j,])%*%solve(sigma[[j]])%*%(X[i,] - mu[j,]))
      }
      if(dist.mat[i]<exp(-lambda))
      {
        C <- C+1
        temp2 <- km.pp(X, C)
        print(temp2[[2]])
        temp3 <- GMM(X, temp2[[1]], temp2[[2]])
        Z <- temp3[[1]]
        prob <- temp3[[4]]
        mu <- temp3[[2]]
        sigma <- temp3[[3]]
        break
      }
    }
    obj.new <- -lambda * C
    for(i in 1:N)
    {
      for(j in 1:C)
      {
        obj.new <- obj.new + Z[i,j]*(log(prob[j])-(X[i,] - mu[j,])%*%solve(sigma[[j]])%*%(X[i,] - mu[j,]))
      }
    }
    if(abs(obj.old-obj.new)<tolerance){
      break}
    
    obj.old <- obj.new
    
  }
  return(list("Z"=Z,"mu"=mu,"prob"=prob))
}