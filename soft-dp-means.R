
dp.soft <- function(X, f, lambda, dphi, epsilon = 1e-3){
   n <- nrow(X)
   d <- ncol(X)
   mu <- t(colMeans(X))
   soft.Z <- as.matrix(rep(1,n))
   hard.Z <- rep(1,n)
   prob <- 1
   count <- 0
   C <- 1
   obj.old <- (-1)
   
   while (count <= 1000) {
      
      indicator <- 0
      
      for (i in 1:n) {
         # durotto <- apply(mu, 1, function(vec) dist.F(X[i,], vec))
         durotto <- apply(mu, 1, function(vec) dphi(X[i,], vec))
         
         if (min(durotto) < lambda) {
            hard.Z[i] <- which.min(durotto)
         }
         else {
            C <- C + 1
            indicator <- indicator + 1
            mu <- matrix(c(as.numeric(t(mu)), X[i,]), C, d, byrow = TRUE)
            hard.Z[i] <- C
         }
      }
      
      
      if (indicator > 0) {
         prob <- as.vector(table(hard.Z))/n
      }
      
      res <- bregman.soft.clus(X, dphi, C, mu, prob, f)
      
      prob <- res$Pi
      soft.Z <- res$Z
      mu <- res$mu
      
      # obj.new <- C * log(sum(prob * apply(mu, 1, 
      #                                     function(vec) f(F.inv(lambda), vec))))
      
      obj.new <- 0
      
      for (i in 1:n) {
         S <- 0
         for (j in 1:C){
            S <- S + prob[j] * f(X[i,], mu[j,])
         }
         obj.new <- obj.new + S
      }
      
      if (abs(obj.new - obj.old) < epsilon) {break}
      
      obj.old <- obj.new
      count <- count + 1
      
   }
   
   return(list("mu" = mu, "Pi" = prob, "Z" = soft.Z, 
               "C" = C))
}


# ############### Binomial
# dphi <- function(u,v)
# {
#    if(u==0)
#       return(N*log(N/(N-v)))
#    else if(u==N)
#       return(N*log(N/v))
#    else
#       return(u*log(u/v)+(N-u)*log((N-u)/(N-v)))
# }
#
# f <- function(u,v)
# {
#    return(choose(N,u)*((v/N)^u)*((N-v)/N)^(N-u))
# }

############### Gaussian
dphi <- function(u,v)
{
   return(((u-v)^2)/2)
}

f <- function(u,v)
{
   return(exp(-((u-v)^2)/2)/sqrt(2*pi))
}


###############
bregman.soft.clus <- function(X, dphi, k, mu, prob, f, epsilon = 10^(-2)){
   
   n <- nrow(X)
   Z <- matrix(rep(0,n*k), n, k)
   
      #### E-Step
      for (i in 1:n) {
         c <- prob * apply(mu, 1, function(vec) exp(-dphi(X[i,], vec)))
         Z[i,] <- c / sum(c)
      }
   
      #### M-Step
      prob <- colMeans(Z)
      for (h in 1:k) {
         mu[h,] <- t(Z[,h]) %*% X / sum(Z[,h])
      }
      
      
   return(list("Z" = Z, "mu" = mu,"Pi"=prob))
}


