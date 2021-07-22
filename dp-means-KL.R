dphi=function(u,v)
{
  return(sum(u*log(u/v,2)-log(exp(1),2)*(u-v)))
}

DP.means.KL <- function(X, lambda, ground = NULL, epsilon = 1e-6){
  
  n <- nrow(X)
  C <- 1
  mu <- matrix(0, 1, ncol(X))
  mu[1,] <- t(colMeans(X))
  Z <- rep(1, n)
  obj.old <- lambda
  
  t <- 1
  while(t<=n)
  {
    obj.old <- obj.old + dphi(X[t,],mu[1,])
    t=t+1
  }
  
  count <- 0
  
  while(TRUE)
  {
    for (i in 1:n)
    {
      dist.vec <- c()
      for (h in 1:C) {dist.vec[h] <- dphi(X[i,],mu[h,])}
      
      if(min(dist.vec) > lambda)
      {
        Z[i] <- C <- C+1
        mu <- matrix(c(as.numeric(t(mu)), X[i,]), nrow = C, byrow = TRUE)
      }
      
      else {Z[i] <- which.min(dist.vec)}
#      print(X[Z==2,])
      for(j in 1:C) {
        if(length(which(Z==j))>1)
        mu[j,] <- colMeans(X[Z==j,])
        else
          mu[j,]<- X[Z==j,]}
    }
    
    print(C)
    count <- count + 1
    
    R <- matrix(0, n, ncol(X))
    for (j in 1:n) {R[j,] <- unlist(mu[Z[j],])}
    
    obj.new <-  lambda * C
    t=1
    while(t<=n)
    {
      obj.new=obj.new+dphi(X[t,],R[t,])
      t=t+1
    }
    if((obj.new - obj.old)^2 < epsilon) {break}
    obj.old <- obj.new
  }
  
  if(is.null(ground) == FALSE){
    
    return(list("Z" = Z , "Number-of-Clusters" = C,
                "No. of Iterations" = count,
                "ARI" = aricode::ARI(Z, ground), 
                "NMI" = aricode::NMI(Z, ground)))
  }
  
  return(list("Z" = Z , "mu" = mu, "Number-of-Clusters" = C, 
              "No. of Iterations" = count))
}

