
bregman.soft <- function(X, dphi, pdensity, k){
   
   d <- ncol(X)
   n <- nrow(X)
   
   Z <- matrix(0, n, k)
   mu <- matrix(0, k, d)
   
   prob <- rep(1/k, k)
   obj.old <- count <- 0
   
   while (count < 1000) {
      #### E-Step
      for (i in 1:n) {
         c <- prob * apply(mu, 1, function(vec) exp(-dphi(X[i,], vec)))
         Z[i,] <- c / sum(c)
      }
      
      #### M-Step
      prob <- colMeans(Z)
      
      for (i in 1:k) {
         mu[h,] <- t(Z[,h]) %*% X / sum(Z[,h])
      }
      
      obj.new <- 1
      
      for (i in 1:n) {
         S <- 0
         for (j in 1:k) {
            S <- S + prob[h] * pdensity(X[i,], mu[h,])
         }
         obj.new <- obj.new * S
      }
      
      if (abs(obj.new - obj.old) < epsilon) {break}
      
      count <- count + 1
      
      return(list("Z" = Z, "mu" = mu))
   }
}

