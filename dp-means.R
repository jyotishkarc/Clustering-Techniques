
DP.means <- function(X, lambda, epsilon){
   
   n <- nrow(X)
   C <- 1
   mu[1,] <- t(colMeans(X))
   Z <- rep(1, n)
   
   obj.old <- sum((X - matrix(rep(mu[1,], n), n, ncol(X), byrow = TRUE))^2) + lambda
   
   while(TRUE)
   {
      for (i in 1:n)
      {
         for (h in 1:C) {dist.vec[h] <- sum((X[i,] - mu[h])^2)}
         
         if(min(dist.vec) > lambda)
         {
            Z[i] <- C <- C+1
            mu[C,] <- t(X[i,])
         }
         
         else {Z[i] <- which.min(dist.vec)}
         
         for(j in 1:C) {mu[j,] <- colMeans(X[Z==j,])}
      }
      
      
      R <- matrix(0,n,ncol(X))
      
      for (j in 1:C) {R[Z==j,] <- mu[j,]}
      
      obj.new <- sum((X-R)^2) + lambda * C
      if((obj.new - obj.old)^2 < epsilon) {break}
      obj.old <- obj.new
   }
   
   return(list("Z" = Z , "mu" = mu))
}
