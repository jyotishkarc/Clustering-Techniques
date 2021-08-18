
sparse.km <- function(X, k, s, initial.mu, ground = NULL, tolerance = 1e-3){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   Z <- matrix(0, nrow = n, ncol = k)
   Z <- c()
   
   dist.mat <- matrix(0, nrow = n, ncol = k)
   
   for (j in 1:k) {
      dist.mat[,j] <- apply(X, 1, function(val) distance.sq(val, initial.mu[j,]))
   }
   
   Z <- apply(dist.mat, 1, which.min)
   # print(Z)
   count <- 0
   
   while (TRUE) {
      
      mu <- new.mu <- matrix(0, k, d)
      
      for (j in 1:k) {
         mu[j,] <- colMeans(X[which(Z==j),])
      }
      
      ranks <- matrix(sapply(1:k, function(val) length(which(Z==val))), 1, k) %*% mu^2
      ranking <- d - rank(ranks)
      
      L <- which(ranking <= s)
      notL <- which(ranking > s)
      
      for (i in 1:n) {
         Z[i] <- which.min(sapply(1:k, function(val) asg.func(X,mu,i,val,L,notL)))
      }
      
      #print(Z)
      
      for (j in 1:k) {
         new.mu[j,] <- colMeans(X[which(Z==j),])
      }
      
      if (sum((mu - new.mu)^2) < tolerance) {break}
      
      mu <- new.mu
      count <- count + 1
      print(count)
   }
   
   if (is.null(ground) == FALSE) {
      ari.clus <- aricode::ARI(ground, Z)
      nmi.clus <- aricode::NMI(ground, Z)
      
      return(list(Z, new.mu, "ARI" = ari.clus, "NMI" = nmi.clus))
   }
   
   return(list(Z, new.mu))
}


################
asg.func <- function(X, mu,i,val,L,notL){
   return(sum((X[i,L] - mu[val,L])^2) + sum(X[i,notL]^2))
}

################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}


