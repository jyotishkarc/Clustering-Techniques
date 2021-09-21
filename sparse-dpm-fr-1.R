## Author : JYOTISHKA RAY CHOUDHURY

sparse.dpm.fr.1 <- function(X, s, lambda, ground = NULL, tolerance = 1e-3){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   mu <- colMeans(X)
   Z <- rep(1, n)
   C <- 1
   t <- 0
   
   obj.old <- sum((X - matrix(rep(mu, n), n, d, byrow = TRUE))^2) + lambda
   
   while(t<=50){
      
      mu <- matrix(0, C, d)
      obj.new <- lambda * C
      
      for (j in 1:C) {
         if(length(which(Z==j))==0)
         {
            mu[j,] <- rep(0, d)
         }
         else if(length(which(Z==j))==1)
         {
            mu[j,] <- X[which(Z==j), ]
         }
         else
         {
            mu[j, ] <- colMeans(X[which(Z==j), ])
         }
      }
      
      
      
      ranks <- matrix(sapply(1:C, function(val) length(which(Z==val))), 1, C) %*% mu^2
      ranking <- d + 1 - rank(ranks)
      
      L <- which(ranking <= s)
      notL <- which(ranking > s)
      
      for (i in 1:n)
      {
         # dist.vec <- c()
         # for (h in 1:C) {dist.vec[h] <- sum((X[i,] - mu[h,])^2)}
         
         dist.vec <- sapply(1:C, function(val) asg.func(X,mu,i,val,L,notL))
         
         if(min(dist.vec) > lambda){
            Z[i] <- C <- C+1
            mu <- matrix(c(as.numeric(t(mu)), X[i,]), nrow = C, byrow = TRUE)
            obj.new <- obj.new + lambda
         }
         else {
            Z[i] <- which.min(dist.vec)
            obj.new <- obj.new + min(dist.vec)
         }
         
         # for(j in 1:C) {
         #    if (length(which(Z == j)) > 1) {
         #       mu[j,] <- colMeans(X[which(Z==j), ])
         #    }
         #    else {mu[j,] <- X[which(Z==j), ]}
         # }
      }
      
      
      v <- sort(unique(Z))
      
      u <- Z
      for(i in 1:(length(v)))
      {
         u[which(Z==v[i])]=i
      }
      Z <- u
      C <- max(Z)
      if(abs(obj.new/obj.old - 1) < tolerance) {break}
      
      obj.old <- obj.new
      t <- t+1
      print(C)
   }
   
   if (is.null(ground) == FALSE) {
      # ari.clus <- aricode::ARI(ground, Z)
      nmi.clus <- aricode::NMI(ground, Z)
      
      # print(ari.clus)
      print(C)
      print(nmi.clus)
      
      return(list("Z" = Z, 
                  "C" = C, 
                  # "ARI" = ari.clus, 
                  "NMI" = nmi.clus))
   }
   
   return(list("Z" = Z, "No of Clusters" = C))
}



################
asg.func <- function(X, mu,i,val,L,notL){
   return(sum((X[i,L] - mu[val,L])^2) + sum(X[i,notL]^2))
}


################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}




