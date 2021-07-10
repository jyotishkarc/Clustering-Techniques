
DP.means <- function(X, lambda, epsilon){
   
   n <- nrow(X)
   C <- 1
   mu <- matrix(0, 1, ncol(X))
   mu[1,] <- t(colMeans(X))
   Z <- rep(1, n)
   
   obj.old <- sum((X - matrix(rep(mu[1,], n), n, ncol(X), byrow = TRUE))^2) + lambda
   
   count <- 0
   
   while(TRUE)
   {
      for (i in 1:n)
      {
         dist.vec <- c()
         for (h in 1:C) {dist.vec[h] <- sum((X[i,] - mu[h,])^2)}
         
         if(min(dist.vec) > lambda)
         {
            Z[i] <- C <- C+1
            mu <- matrix(c(as.numeric(t(mu)), X[i,]), nrow = C, byrow = TRUE)
         }
         
         else {Z[i] <- which.min(dist.vec)}
         
         for(j in 1:C) {mu[j,] <- colMeans(X[Z==j,])}
      }
      
      print(C)
      count <- count + 1
      
      R <- matrix(0, n, ncol(X))
      for (j in 1:n) {R[j,] <- unlist(mu[Z[j],])}
      
      obj.new <- sum((X-R)^2) + lambda * C
      if((obj.new - obj.old)^2 < epsilon) {break}
      obj.old <- obj.new
   }
   
   return(list("Z" = Z , "mu" = mu, 
               "Number-of-Clusters" = C, "Iterations" = count))
}

