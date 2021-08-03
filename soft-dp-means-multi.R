
dp.soft.multi <- function(X, lambda, epsilon = 1e-3){
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   mu <- t(colMeans(X))
   sig <- 1
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
         durotto <- apply(mu, 1, function(vec) dphi.multi(X[i,], vec))
         
         if (min(durotto) < lambda) {
            hard.Z[i] <- which.min(durotto)
         }
         else {
            C <- C + 1
            indicator <- indicator + 1
            mu <- matrix(c(as.numeric(t(mu)), X[i,]), C, d, byrow = TRUE)
            hard.Z[i] <- C
         }
         
          if(i %% 500 == 0) { print(i) }
      }
      
      
      if (indicator > 0) { prob <- as.vector(table(hard.Z))/n }
      
      res <- bregman.soft.clus(X, C, mu, prob)
      
      prob <- res$Pi
      soft.Z <- res$Z
      mu <- res$mu
      
      # obj.new <- C * log(sum(prob * apply(mu, 1, 
      #                             function(vec) f.multi(F.inv(lambda), vec))))
      
      obj.new <- 0
      
      for (i in 1:n) {
         S <- 0
         for (j in 1:C){
            S <- S + prob[j] * f.multi(X[i,], mu[j,])
         }
         obj.new <- obj.new + S
      }
      
      if (abs(obj.new - obj.old) < epsilon) {break}
      
      obj.old <- obj.new
      count <- count + 1
      
      print("Iteration")
   }
   
   return(list("mu" = mu, "Pi" = prob, "Z" = soft.Z, 
               "C" = C))
}


# ############### Binomial
# dphi.multi <- function(u,v)
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
dphi.multi <- function(u, v, sig = 1){
   return(distance.sq(u,v)/(2 * sig^2))
}

f.multi <- function(u, v, sig = 1){
   d <- length(u)
   return(exp(-distance.sq(u,v)/(2 * sig^2))/(sqrt(2*pi))^d)
}


###############
bregman.soft.clus <- function(X, k, mu, prob, epsilon = 1e-2){
   
   n <- nrow(X)
   Z <- matrix(rep(0,n*k), n, k)
   
   #### E-Step
   for (i in 1:n) {
      g <- prob * apply(mu, 1, function(vec) exp(-dphi.multi(X[i,], vec)))
      Z[i,] <- g / sum(g)
   }
   
   #### M-Step
   prob <- colMeans(Z)
   
   for (h in 1:k) {
      mu[h,] <- t(Z[,h]) %*% X / sum(Z[,h])
      
      for (i in 1:N) {
         # new.mu.matrix[,h] <- new.mu.matrix[,h] + gamma.prob[i,h] * dm[i,]
         
         new.Sigma.list[[h]] <- new.Sigma.list[[h]] + Z[i,h] * 
            (as.matrix(X[i,]) - as.matrix(mu[,h])) %*% 
            t(X[i,] - mu[,h])
      }
      
      new.Sigma.list[[h]] <- new.Sigma.list[[h]] / sum(Z[,h])
   }
   
   
   return(list("Z" = Z, "mu" = mu, "Sigma" = new.Sigma.list, "Pi"=prob))
}


