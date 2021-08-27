SKFR.perm <- function(X, k, B){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   X.perm <- list()
   
   for(b in 1:B){
      X.perm[[b]] <- apply(X, 2, 
            function(column) sample(column, length(column), replace = FALSE))
   }
   
   gap.stat <- c()
   # initial.mu <- t(initial(t(X), k))
   
   for(s in 1:d){
      initial.mu <- t(initial(t(X), k))
      
      O.stat.X <- O.stat(X, k, s, initial.mu)
      print(O.stat.X)
      O.stat.X.perm <- sapply(1:B, function(b) O.stat(X.perm[[b]], k, s, initial.mu))
      if(s == 1) print(O.stat.X.perm)
      # print(O.stat.X.perm)
      gap.stat[s] <- log(O.stat.X) - mean(log(O.stat.X.perm))
   }
   
   print(gap.stat)
   return(which.max(gap.stat))
}



O.stat <- function(X, k, s, initial.mu){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   X.bar <- matrix(rep(colMeans(X), n), n, d, byrow = TRUE)
   
   result <- sparse.km.fr.1(X, k, s, initial.mu, ground = NULL, tolerance = 1e-3)
   Z <- result[[1]]
   asg.vec <- result[[2]]
   centroid <- matrix(0, k, d)
   
   for(i in 1:k){
      if(length(which(asg.vec == i)) == 1){
         centroid[i,] <- X[i,]
      }
      else if(length(which(asg.vec == i)) == 1){
         centroid[i,] <- rep(0, d)
      }
      else centroid[i,] <- colMeans(X[which(asg.vec == i),])
   }
   
   return(sum((X - X.bar)^2) - sum((X - Z %*% centroid)^2))
}
