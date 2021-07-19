
DP.means <- function(X, lambda, ground = NULL, epsilon = 1e-6){
   
   n <- nrow(X)
   d <- ncol(X)
   C <- 1
   mu <- matrix(0, 1, d)
   mu[1,] <- t(colMeans(X))
   Z <- rep(1, n)
   
   obj.old <- sum((X - matrix(rep(mu[1,], n), n, d, byrow = TRUE))^2) + lambda
   
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
      
      count <- count + 1
      
      R <- matrix(0, n, d)
      for (j in 1:n) {R[j,] <- unlist(mu[Z[j],])}
      
      obj.new <- sum((X-R)^2) + lambda * C
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

