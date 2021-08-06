
dp.soft.multi <- function(X, lambda, epsilon = 1e-3){
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   mu <- t(colMeans(X))
   Sigma <- (t(X) - colMeans(X)) %*% t(t(X) - colMeans(X))/n
   Sigma <- list(Sigma)
   
   print(Sigma)
   
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
         # durotto <- lapply(list(mu,Sigma), function(vec, mat) 
         #                                     dphi.multi(X[i,], vec, ))
         
         durotto <- c()
         for (j in 1:C) {
            durotto[j] <- dphi.multi(X[i,], mu[j,], Sigma[[j]])
         }
         
         if (min(durotto) < lambda) {
            hard.Z[i] <- which.min(durotto)
         }
         else {
            C <- C + 1
            indicator <- indicator + 1
            mu <- matrix(c(as.numeric(t(mu)), X[i,]), C, d, byrow = TRUE)
            hard.Z[i] <- C
            Sigma[[]]
         }
         
          if(i %% 500 == 0) { print(i) }
      }
      
      
      if (indicator > 0) { prob <- as.vector(table(hard.Z))/n }
      
      res <- bregman.soft.clus(X, C, mu, Sigma, prob)
      
      prob <- res$Pi
      soft.Z <- res$Z
      mu <- res$mu
      Sigma <- res$Sigma
      
      # obj.new <- C * log(sum(prob * apply(mu, 1, 
      #                             function(vec) f.multi(F.inv(lambda), vec))))
      
      obj.new <- 0
      
      for (i in 1:n) {
         S <- 0
         for (j in 1:C){
            S <- S + prob[j] * f.multi(X[i,], mu[j,], Sigma[[j]])
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
dphi.multi <- function(vec, mu, sigma.mat){
   return((vec-mu) * solve(sigma.mat) * (vec-mu)/2)
}

f.multi <- function(vec, mu, sigma.mat){
   dmvnorm(vec, mu, sigma.mat)
}

# f.multi <- function(u, v, mat){
#    d <- length(u)
#    return(exp(-(u-v) * solve(mat) * (u-v)/2)/(sqrt(2*pi))^d)
# }


###############
bregman.soft.clus <- function(X, k, mu, Sigma.list, prob, epsilon = 1e-2){
   
   n <- nrow(X)
   Z <- matrix(rep(0,n*k), n, k)
   
   #### E-Step
   for (i in 1:n) {
      g <- c()
      for (j in 1:k) {
         g[j] <- prob[j] * exp(-dphi.multi(X[i,], mu[j,], Sigma.list[[j]]))
      }
      
      Z[i,] <- g / sum(g)
   }
   
   # for (i in 1:n) {
   #    g <- prob * apply(mu, 1, function(vec) exp(-dphi.multi(X[i,], vec)))
   #    Z[i,] <- g / sum(g)
   # }
   
   
   #### M-Step
   prob <- colMeans(Z)
   Sigma.list <- list()
   for (h in 1:k) {
      mu[h,] <- t(Z[,h]) %*% X / sum(Z[,h])
      
      for (i in 1:N) {
         print(dim(as.matrix(X[i,])))
         print(dim(as.matrix(mu[,h])))
         Sigma.list[[h]] <- Sigma.list[[h]] + Z[i,h] * 
            (as.matrix(X[i,]) - as.matrix(mu[,h])) %*% 
            t(X[i,] - mu[,h])
      }
      
      Sigma.list[[h]] <- Sigma.list[[h]] / sum(Z[,h])
   }
   
   
   return(list("Z" = Z, "mu" = mu, "Sigma" = Sigma.list,  "Pi"=prob))
}

