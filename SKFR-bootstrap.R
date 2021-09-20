SKFR.boot <- function(X, k, B){
   
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   X.boot <- list()
   
   for(b in 1:B){
      X.boot[[b]] <- apply(X, 2, 
                           function(column){sample(column, 
                                                   length(column), 
                                                   replace = F)})
   }
   
   gap.stat <- c()
   # initial.mu <- t(initial(t(X), k))
   
   for(s in 1:d){
      print(s)
      
      initial.mu <- t(initial(t(X), k))
      O.stat.X <- O.stat(X, k, s, initial.mu)
      print(O.stat.X)
      O.stat.X.boot <- sapply(1:B, function(b) O.stat(X.boot[[b]], k, s, initial.mu))
      gap.stat[s] <- log(O.stat.X) - mean(log(O.stat.X.boot))
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
   L=result$features
   mu=matrix(0,n,d)
   for(i in 1:n)
      mu[i,L]=colMeans(X)[L]
   Z <- result[[1]]
   asg.vec <- result[[2]]
   centroid <- matrix(0, k, d)
   
   for(i in 1:k){
      if(length(which(asg.vec == i)) == 1){
         # centroid[i,] <- X[i,]
         centroid[i,] <- X[which(asg.vec == i),]
      }
      else if(length(which(asg.vec == i)) == 0){
         centroid[i,] <- rep(0, d)
      }
      else centroid[i,] <- colMeans(X[which(asg.vec == i),])
   }
   
   return(sum((X-mu)^2)-result$obj)
}
