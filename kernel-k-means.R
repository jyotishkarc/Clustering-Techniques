kernel.km <- function(data , k , initial, tolerance = 1e-6){
   X <- as.matrix(data)
   n <- nrow(X)
   d <- ncol(X)
   centroid <- initial
   Z <- matrix(0, nrow = n , ncol = k)
   
   for (i in 1:n) {
      dist.mat <- rep(0,k)
      for (j in 1:k) {
         dist.mat[j] <-  distance.sq(X[i,] , centroid[j,])
      }
      
      Z[i, which.min(dist.mat)] <- 1
   }
   
   new.centroid <- matrix(0, nrow = k , ncol = d)
   
   for (i in 1:k) {
      h <- which(Z[,i] == 1)
      qq <- as.matrix(X[h,])
      new.centroid[i,] <- colMeans(qq)
   }
   
   count <- 0
   obj.old <- (norm(X - Z %*% centroid, type = "F"))^2
   
   while (TRUE) {
      
      for (i in 1:n) {
         dist.mat <- rep(0,k)
         for (j in 1:k) {
            dist.mat[j] <- distance.sq(X[i,] , centroid[j,])
         }
         
         Z[i,] <- rep(0,k)
         Z[i, which.min(dist.mat)] <- 1
      }
      
      for (i in 1:k) {
         h <- which(Z[,i]==1)
         qq <- as.matrix(X[h,])
         new.centroid[i,] <- colMeans(qq)
      }
      
      obj.new <- (norm(X - Z %*% new.centroid, type = "F"))^2
      
      if (abs(obj.new - obj.old) < tolerance) {break}
      
      obj.old <- obj.new
      centroid <- new.centroid
      count <- count + 1
   }
   
   cluster.asg.vector <- apply(Z, 1, function(val) which(val == 1))
   
   result <- list(cluster.asg.vector, colSums(Z), count, new.centroid)
   return(result)
}


################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}


################
