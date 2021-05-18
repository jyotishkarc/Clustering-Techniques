# Author : Jyotishka Ray Choudhury

minmax <- function(data , k , initial , p = 0.5){
   dm <- as.matrix(data)
   m <- nrow(dm)
   n <- ncol(dm)
   centroid <- initial
   Z <- matrix(0, nrow = k , ncol = n)
   
   dist.mat <- matrix(0 , nrow = k , ncol = n)
   for (j in 1:n) {
      for (i in 1:k) {
         dist.mat[i,j] <- distance.sq(dm[,j] , centroid[,i])
      }
      Z[which.min(dist.mat[,j]) , j] <- 1
   }
   
   V <- rep(0,k)
   for(i in 1:k){
      V[i] <- sum(Z[i,] * dist.mat[i,])
   }
   
#   W <- V^(1/(1-p))/sum(V^(1/(1-p)))
   W <- rep(1/k , k)
   print(W)
   
   Z <- matrix(0, nrow = k , ncol = n)
   for (j in 1:n) {
      bb <- c()
      for (i in 1:k) {
         bb[i] <- W[i]^p * dist.mat[i,j]
      }
      Z[which.min(bb) , j] <- 1
   }
   
   new.centroid <- matrix(0, nrow = m , ncol = k)
   for (i in 1:k) {
      d <- which(Z[i,]==1)
      qq <- as.matrix(dm[,d])
      new.centroid[,i] <- rowMeans(qq)
   }
   
   count <- 0
   tolerance <- 1E-4
   
   while (distance.sq(centroid , new.centroid) > tolerance) {
      centroid <- new.centroid
      
      for (j in 1:n) {
         for (i in 1:k) {
            dist.mat[i,j] <- distance.sq(dm[,j] , centroid[,i])
         }
         Z[which.min(dist.mat[,j]) , j] <- 1
      }
      
      for(i in 1:k){
         V[i] <- sum(Z[i,] * dist.mat[i,])
      }
      
      W <- V^(1/(1-p))/sum(V^(1/(1-p)))
      print(W)
      
      Z <- matrix(0, nrow = k , ncol = n)
      for (j in 1:n) {
         bb <- c()
         for (i in 1:k) {
            bb[i] <- W[i]^p * dist.mat[i,j]
         }
         Z[which.min(bb) , j] <- 1
      }
      
      for (i in 1:k) {
         d <- which(Z[i,]==1)
         qq <- as.matrix(dm[,d])
         new.centroid[,i] <- rowMeans(qq)
      }
      
      count <- count + 1
   }
   
   final.cluster <- which(Z == 1) - k * 0:(ncol(data)-1)
   
   result <- list(rowSums(Z),Z,final.cluster, count)
   return(result)
}




################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}