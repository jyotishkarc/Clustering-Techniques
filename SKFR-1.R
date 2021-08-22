## Author : JYOTISHKA RAY CHOUDHURY

sparse.km.fr.1 <- function(X, k, s, initial.mu, ground = NULL, tolerance = 1e-3){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   Z <- matrix(0, nrow = n, ncol = k)
   Z <- c()
   mu <- initial.mu
   dist.mat <- matrix(0, nrow = n, ncol = k)
   
   for (j in 1:k) {
      dist.mat[,j] <- apply(X, 1, function(val) distance.sq(val, initial.mu[j,]))
   }
   
   Z <- apply(dist.mat, 1, which.min)
   # print(Z)
   obj.old <- 0
   for(i in 1:n)
   {
      dist.mat <- rep(0, k)
      for(j in 1:k)
      {
         dist.mat[j] <- sum((X[i, ] - mu[j, ])^2)
      }
      obj.old <- obj.old + min(dist.mat)
   }
   
   count <- 0
   
   while (TRUE) {
      
      mu <- matrix(0, k, d)
      obj.new <- 0
      
      for (j in 1:k) {
         mu[j,] <- colMeans(X[which(Z==j),])
      }
      
      ranks <- matrix(sapply(1:k, function(val) length(which(Z==val))), 1, k) %*% mu^2
      ranking <- d + 1 - rank(ranks)
      
      L <- which(ranking <= s)
      notL <- which(ranking > s)
      
      for (i in 1:n) {
         dist.vec <- sapply(1:k, function(val) asg.func(X,mu,i,val,L,notL))
         Z[i] <- which.min(dist.vec)
         obj.new <- obj.new + min(dist.vec)
      }
      
      if ((obj.new - obj.old)^2 < tolerance) {break}
      
      obj.old <- obj.new
      count <- count + 1
   }
   
   if (is.null(ground) == FALSE) {
      # ari.clus <- aricode::ARI(ground, Z)
      nmi.clus <- aricode::NMI(ground, Z)
      
      # print(ari.clus)
      print(nmi.clus)
      
      return(list(Z, "ARI" = ari.clus, "NMI" = nmi.clus))
   }
   
   return(list(Z))
}


################
asg.func <- function(X, mu,i,val,L,notL){
   return(sum((X[i,L] - mu[val,L])^2) + sum(X[i,notL]^2))
}

################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}


