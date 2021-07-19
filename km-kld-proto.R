f=function(x,y)
{
  s=sum(x*log(x/y,2)-log(exp(1),2)*(x-y))
  return(s)
}

km.kl <- function(data , k , initial){
   dm <- as.matrix(data)
   m <- nrow(dm)
   n <- ncol(dm)
   centroid <- initial
   Z <- matrix(0, nrow = k , ncol = n)
   
   for (i in 1:n) {
      dist.mat <- rep(0,k)
      for (j in 1:k) {
         dist.mat[j] <-  distance.sq(dm[,i] , centroid[,j])
      }
      Z[which.min(dist.mat) , i] <- 1
   }
   
   new.centroid <- matrix(0, nrow = m , ncol = k)
   for (i in 1:k) {
      d <- which(Z[i,]==1)
      qq <- as.matrix(dm[,d])
      new.centroid[,i] <- rowMeans(qq)
   }
   
   count <- 0
   tolerance <- 1E-8
   while (distance.sq(centroid , new.centroid) > tolerance) {
      centroid <- new.centroid
      
      for (i in 1:n) {
         dist.mat <- rep(0,k)
         for (j in 1:k) {
            dist.mat[j] <- distance.sq(dm[,i] , centroid[,j])
         }
         
         Z[,i] <- rep(0,k)
         Z[which.min(dist.mat) , i] <- 1
      }
      
      for (i in 1:k) {
         d <- which(Z[i,]==1)
         qq <- as.matrix(dm[,d])
         new.centroid[,i] <- rowMeans(qq)
      }
      
      count <- count + 1
   }
   
   result <- list(rowSums(Z), Z, count, new.centroid)
   return(result)
}


################
distance.sq <- function(x,y){
   return(sum(x*log(x/y,2)-log(exp(1),2)*(x-y)))
}

