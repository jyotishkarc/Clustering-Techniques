## Author : JYOTISHKA RAY CHOUDHURY

GMM <- function(data, K){
   dm <- as.matrix(data)
   initial <- initial.clusters(dm, K)
   dm <- t(dm)
   N <- nrow(dm)
   
   E <- E.step(dm, K, initial[[1]], initial[[2]], initial[[3]])
   Q <- E[[2]]
   difference <- 1
   log.likelihood <- c(Q)
   while (difference > 0.2) {
      M <- M.step(dm, K, E)
      Q <- M[[4]]
      
      new.E <- E.step(dm, K, M[[1]], M[[2]], M[[3]])
      
      new.Q <- new.E[[2]]
      log.likelihood <- c(log.likelihood, new.Q)
      difference <- abs(Q - new.Q)
      print(new.Q)
      print(difference)
      E <- new.E
   }
   
   gamma.prob <- matrix(0, nrow = N, ncol = K)
   pi.prob <- E[[5]]
   mu.matrix <- E[[3]]
   Sigma.list <- E[[4]]
   
   for (i in 1:N) {
      S <- 0
      for (j in 1:K) {
         S <- S + pi.prob[j] * normal.density(dm[i,], mu.matrix[,j], 
                                              Sigma.list[[j]])
      }
      
      for (j in 1:K) {
         gamma.prob[i,j] <- pi.prob[j] * normal.density(dm[i,], mu.matrix[,j],
                                                        Sigma.list[[j]]) / S
      }
   }
   
   final <- list(gamma.prob, pi.prob, mu.matrix, Sigma.list, log.likelihood)
   
   return(final)
}



########## Initialization

initial.clusters <- function(data,K){
   dm <- as.matrix(data)
   d <- nrow(dm)
   n <- ncol(dm)
   centroid.mat <- matrix(0,nrow = d,ncol = K)
   centroid.mat[,1] <- dm[,sample(1:n , 1)]
   
   for (h in 2:K) {
      dist.mat <- matrix(0 , nrow = n , ncol = h-1)
      
      for (i in 1:n) {
         for (j in 1:h-1) {
            dist.mat[i,j] <- distance.sq(dm[,i] , centroid.mat[,j])
         }
      }
      
      min.mat <- c()
      for(q in 1:n){
         min.mat[q] <- min(dist.mat[q,])
      }
      
      selection.prob <- min.mat / sum(min.mat)
      ss <- sample(1:n , 1 , prob = selection.prob)
      centroid.mat[,h] <- dm[,ss]
   }
   
   Z <- matrix(0, nrow = K , ncol = n)
   
   for (i in 1:n) {
      dist.mat <- rep(0,K)
      for (j in 1:K) {
         dist.mat[j] <-  distance.sq(dm[,i] , centroid.mat[,j])
      }
      Z[which.min(dist.mat) , i] <- 1
   }
   
   initial.pi.prob <- rowSums(Z)/sum(Z)
   
   Sigma.list <- list()
   
   for (i in 1:K) {
      mark <- which(Z[i,] == 1)
      avg <- as.matrix(rowMeans(dm[,mark]))
      mat <- as.matrix(dm[,mark]) - matrix(rep(avg,length(mark)), 
                                           ncol = length(mark))
      J <- 0
      
      for (h in 1:length(mark)) {
         J <- J + as.matrix(mat[,h]) %*% t(mat[,h])
      }
      J <- J / length(mark)
      
      Sigma.list[[i]] <- J
   }
   
   initial.parameters <- list(initial.pi.prob, centroid.mat, Sigma.list)
   return(initial.parameters)
}



########## Expectation Step

E.step <- function(data, K, pi.prob, mu.matrix, Sigma.list){
   dm <- as.matrix(data)
   N <- nrow(dm)
   D <- ncol(dm)
   gamma.prob <- matrix(0, nrow = N, ncol = K)
   
   for (i in 1:N) {
      S <- 0
      for (j in 1:K) {
         S <- S + pi.prob[j] * normal.density(dm[i,], mu.matrix[,j], 
                                           Sigma.list[[j]])
      }
      
      for (j in 1:K) {
         gamma.prob[i,j] <- pi.prob[j] * normal.density(dm[i,], mu.matrix[,j],
                                                       Sigma.list[[j]]) / S
      }
      
   }
   
   Q <- C <- 0
   for (i in 1:N) {
      Q <- 0
      for (j in 1:K) {
#         Q <- Q + log((pi.prob[j] * 
#      normal.density(dm[i,], mu.matrix[,j], Sigma.list[[j]]))^gamma.prob[i,j])
      
         Q <- Q + pi.prob[j] * 
            normal.density(dm[i,], mu.matrix[,j], Sigma.list[[j]])
      }
      
      C <- C + log(Q)
   }
   
   result <- list(gamma.prob , C , mu.matrix , Sigma.list, pi.prob)
   return(result)
}



########## Maximization Step

M.step <- function(data, K, res){
   dm <- as.matrix(data)
   N <- nrow(dm)
   D <- ncol(dm)
   gamma.prob <- res[[1]]
   mu.matrix <- res[[3]]
   S <- colSums(gamma.prob)
   
   new.pi.prob <- numeric(K)
   new.mu.matrix <- matrix(0, nrow = D, ncol = K)
   new.Sigma.list <- list()
   
   for (j in 1:K) {
      new.Sigma.list[[j]] <- 0
   }
   
   for (j in 1:K) {
      new.pi.prob[j] <- S[j] / N
      
      for (i in 1:N) {
         new.mu.matrix[,j] <- new.mu.matrix[,j] + gamma.prob[i,j] * dm[i,]
         
         new.Sigma.list[[j]] <- new.Sigma.list[[j]] + gamma.prob[i,j] * 
            (as.matrix(dm[i,]) - as.matrix(mu.matrix[,j])) %*% 
            t(dm[i,] - mu.matrix[,j])
      }
      new.mu.matrix[,j] <- new.mu.matrix[,j] / S[j]
      new.Sigma.list[[j]] <- new.Sigma.list[[j]] / S[j]
   }
   
   result <- list(new.pi.prob, new.mu.matrix, new.Sigma.list, res[[2]])
   return(result)
}



################ Euclidean Norm

distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}



########## Multivariate Normal Density Function

normal.density <- function(vec, mu, Sigma){
   D <- length(vec)
   vec <- as.matrix(vec)
   mu <- as.matrix(mu)
   Sigma <- as.matrix(Sigma)
   
   f <- exp(-1/2 * t(vec - mu) %*% solve(Sigma) %*% (vec - mu)) /
      sqrt((2 * pi)^D * det(Sigma))
   return(as.numeric(f))
}



########### Plotting the Ellipses

if(FALSE){
   data <-                                            # DATASET NAME
   data_ground <- as.numeric(as.matrix())             # GROUND TRUTH
   
   RES <- GMM(t(data), K)                             # GMM Implementation
   
   plot(data[,1], data[,2], col="black", pch=20)      # Plot original dataset
   
   i=1
   while(i <= K)
   {
      mixtools::ellipse(RES[[3]][,i], RES[[4]][[i]], 0.01, col="blue", lwd = 2)
      i=i+1
   }
   
   ARI_cluster <- fclust::ARI.F(data_ground, RES[[1]])
   print(ARI_cluster)
}



