
DP.means <- function(X, lambda, ground = NULL, epsilon = 1e-3){
   
   n <- nrow(X)
   d <- ncol(X)
   C <- 1
   mu <- matrix(0, 1, d)
   mu[1,] <- t(colMeans(X))
   Z <- rep(1, n)
   
   obj.old <- sum((X - matrix(rep(mu[1,], n), n, d, byrow = TRUE))^2) + lambda
   
   count <- 0
   
   while(TRUE)
   {
      for (i in 1:n)
      {
         dist.vec <- c()
         for (h in 1:C) {dist.vec[h] <- sum((X[i,] - mu[h,])^2)}
         
         if(min(dist.vec) > lambda)
         {
            Z[i] <- C <- C+1
            mu <- matrix(c(as.numeric(t(mu)), X[i,]), nrow = C, byrow = TRUE)
         }
         
         else {Z[i] <- which.min(dist.vec)}
         
         for(j in 1:C) {
            if (length(which(Z == j)) > 1) {
               mu[j,] <- colMeans(X[which(Z==j), ])
            }
            else {mu[j,] <- X[which(Z==j), ]}
         }
      }
      
      count <- count + 1
      print(count)
      
      R <- matrix(0, n, d)
      for (j in 1:n) {R[j,] <- unlist(mu[Z[j],])}
      
      obj.new <- sum((X-R)^2) + lambda * C
      if((obj.new - obj.old)^2 < epsilon) {break}
      obj.old <- obj.new
   }
   
   
   ############ 2D plot
   if(ncol(X) == 2){
      print("HELLO")
      
      # original.plot <- ggplot(X , aes(V1,V2)) + 
      #    geom_point(size = 2) + theme(legend.position = "none")
      # 
      # print(original.plot)
      
      fg <- c()
      for (k in 1:nrow(X)) {
         fg[k] <- grDevices::rainbow(length(unique(Z)))[Z[k]]
      }
      
      dpm.plot <- ggplot(X , aes(V1,V2 , color = fg)) +
                  geom_point(size = 2) + theme(legend.position = "none")
      
      print(dpm.plot)
   }
   
   
   if(is.null(ground) == FALSE){
      print("Done")
      
      return(list("Z" = Z , "Number-of-Clusters" = C,
                  "No. of Iterations" = count,
                  "ARI" = aricode::ARI(Z, ground), 
                  "NMI" = aricode::NMI(Z, ground)))
   }
   
   print("Done")

   return(list("Z" = Z , "mu" = mu, "Number-of-Clusters" = C, 
               "No. of Iterations" = count))
}

