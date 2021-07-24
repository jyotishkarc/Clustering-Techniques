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
# pdensity <- function(u,v)
# {
#    return(choose(N,u)*((v/N)^u)*((N-v)/N)^(N-u))
# }


dphi <- function(u,v)
{
   return(((u-v)^2)/2)
}


pdensity <- function(u,v)
{
   return(exp(-((u-v)^2)/2)/sqrt(2*pi))
}

bregman.soft <- function(X, dphi, pdensity, k, epsilon = 10^(-2)){
   
   d <- ncol(X)
   n <- nrow(X)
   
   Z <- matrix(rep(0,n*k), n, k)
   mu <- matrix(0, k, d)
   g <- sample(1:n,k,replace=FALSE)
   for(i in 1:k){
      mu[i,]=X[g[i],]
   }
   prob <- rep(1/k, k)
   obj.old <- count <- 0
   obj.old <- -1
   while (count < 1000) {
      #### E-Step
      for (i in 1:n) {
         c <- prob * apply(mu, 1, function(vec) exp(-dphi(X[i,], vec)))
         Z[i,] <- c / sum(c)
      }
      #### M-Step
      prob <- colMeans(Z)
      print(prob)
      for (h in 1:k) {
         mu[h,] <- t(Z[,h]) %*% X / sum(Z[,h])
      }
      
      obj.new <- 0
      
      for (i in 1:n) {
         S <- 0
         for (h in 1:k) {
            S <- S + prob[h] * pdensity(X[i,], mu[h,])
         }
         obj.new <- obj.new + log(S)
      }
      if (abs(obj.new - obj.old) < epsilon) {break}
      obj.old <- obj.new
      count <- count + 1
   }
   print(count)
   return(list("Z" = Z, "mu" = mu,"T"=count,"Pi"=prob))
}

