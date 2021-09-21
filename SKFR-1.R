## Author : JYOTISHKA RAY CHOUDHURY

sparse.km.fr.1 <- function(X, k, s, initial.mu, ground = NULL, tolerance = 1e-3){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   asg.vec <- matrix(0, nrow = n, ncol = k)
   asg.vec <- c()
   Z <- matrix(0, n, k)
   mu <- initial.mu
   dist.mat <- matrix(0, nrow = n, ncol = k)
   
   for (j in 1:k) {
      dist.mat[,j] <- apply(X, 1, function(val) distance.sq(val, initial.mu[j,]))
   }
   
   asg.vec <- apply(dist.mat, 1, which.min)
   # print(asg.vec)
   obj.old <- 0
   for(i in 1:n)
   {
      # dist.mat <- rep(0, k)
      # for(j in 1:k)
      # {
      #    dist.mat[j] <- sum((X[i, ] - mu[j, ])^2)
      # }
      obj.old <- obj.old + min(dist.mat)
   }
   
   count <- 0
   
   while (count<=500) {
      
      mu <- matrix(0, k, d)
      obj.new <- 0
      
      for (j in 1:k) {
         if(length(which(asg.vec==j))==1)
            mu[j,] <- X[which(asg.vec==j),]
         else if(length(which(asg.vec==j))==0)
            mu[j,] <- rep(0, d)
         else
            mu[j,] <- colMeans(X[which(asg.vec==j),])
      }
      
      ranks <- matrix(sapply(1:k, function(val) length(which(asg.vec==val))), 1, k) %*% mu^2
      ranking <- d+1-rank(as.numeric(ranks))
      
      L <- which(ranking <= s)
      notL <- which(ranking > s)
      
      for (i in 1:n) {
         dist.vec <- sapply(1:k, function(val) asg.func(X,mu,i,val,L,notL))
         asg.vec[i] <- which.min(dist.vec)
         obj.new <- obj.new + min(dist.vec)
      }
      
      if (abs(obj.new /obj.old- 1) < tolerance) {break}
      
      obj.old <- obj.new
      count <- count + 1
      #print(obj.new)
   }
   
   for(i in 1:n){
      Z[i, asg.vec[i]] <- 1
   }
   # print(asg.vec)
   if (is.null(ground) == FALSE) {
      # ari.clus <- aricode::ARI(ground, asg.vec)
      nmi.clus <- aricode::NMI(ground, asg.vec)
      
      # print(ari.clus)
      print(nmi.clus)
      
      return(list("Z" = Z, "vec" = asg.vec, "NMI" = nmi.clus, "obj"=obj.new))
   }
   
   return(list("Z" = Z, "vec" = asg.vec, "obj"=obj.new,"features"=L))
}


################
asg.func <- function(X, mu,i,val,L,notL){
   return(sum((X[i,L] - mu[val,L])^2) + sum(X[i,notL]^2))
}

################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}


