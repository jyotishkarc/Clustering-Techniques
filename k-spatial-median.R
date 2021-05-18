### Author : Jyotishka Ray Choudhury
### Designation : Student, B.Stat 2nd year


require("ICSNP")
max.iter <- 200

bis.k.spatial <- function(k, data){
   dm <- as.matrix(data)
   m <- nrow(dm)
   # n <- ncol(dm)
   # C <- spatial.median(dm)
   # C.L <- dm[,runif(1, min = 1, max = ncol(dm))]
   # C.R <- 2 * C - C.L
   
   # Z <- matrix(0, nrow = k , ncol = n)
   
   S <- list(dm)
   
   for(i in 1:(k-1)){
      
      RAD.mat <- c()
      for (f in 1:i) {
            RAD.mat[f] <- RAD(S[[f]])
      }
      
#      print(RAD.mat)
      
      U <- which.max(RAD.mat)
      J <- subc.create(S[[U]])
      S[[U]] <- NULL
      S <- c(S , list(J[[1]]), list(J[[2]]))
      
      print(S)
   }
   
   return(S)
}

################
distance.sq <- function(x,y){
   return(sum((x-y)*(x-y)))
}

distance.L1 <- function(x,y){
   return(sqrt(distance.sq(x,y)))
}


################
SPD <- function(v, X){
   
   v <- as.matrix(v)
   X <- as.matrix(X)
   sum.SSF <- c()
   
   for (i in 1:ncol(X)) {
      
      X[,i] <- as.matrix(X[,i])
      if (v == X[,i]) {next}
      
      sum.SSF <- sum.SSF + (v - X[,i])/distance.L1(0, v - X[,i])
   }
   
   return(1 - distance.L1(0 , sum.SSF/ncol(X)))
}


#################
subc.create <- function(C.mat){
   
   C.mat <- as.matrix(C.mat)
   n <- ncol(C.mat)
   C <- spatial.median(C.mat)
   C.L <- C.mat[,sample(1:ncol(C.mat) , size = 1)]
   C.R <- 2 * C - C.L
   
   for (j in 1 : max.iter) {
      
      dist.L <- dist.R <- 0
      
      for(h in 1:n){
         dist.L[h] <- distance.sq(C.mat[,h] , C.L)
         dist.R[h] <- distance.sq(C.mat[,h] , C.R)
      }
      
      alpha <- which(dist.L < dist.R)
      beta <- which(dist.L >= dist.R)
      
      M.L <- spatial.median(C.mat[,alpha])
      M.R <- spatial.median(C.mat[,beta])
      
      
      if (distance.sq(M.L , C.L) < tol && distance.sq(M.R , C.R) < tol) {
         J.left <- C.mat[,alpha]
         J.right <- C.mat[,beta]
         J <- list(left = J.left , right = J.right)
         return(J)
         break
      }
      
      C.L <- M.L
      C.R <- M.R
   }
   
   J.left <- C.mat[,alpha]
   J.right <- C.mat[,beta]
   J <- list(J.left , J.right)
   return(J)
}


#################
RAD <- function(C.mat){
   
   C.mat <- as.matrix(C.mat)
   
   J <- subc.create(C.mat)
   C1 <- J[[1]]
   C2 <- J[[2]]
   
   sum.W.1 <- sum.B.1 <- sum.W.2 <- sum.B.2 <- 0
   
   for (i in 1:ncol(C1)) {
      sum.W.1 <- sum.W.1 + SPD(C1[,i] , C1)
      sum.B.1 <- sum.B.1 + SPD(C1[,i] , C2)
   }
   
   for (i in 1:ncol(C2)) {
      sum.W.2 <- sum.W.2 + SPD(C2[,i] , C2)
      sum.B.2 <- sum.B.2 + SPD(C2[,i] , C1)
   }
   
   return((sum.W.1 - sum.B.1)/ncol(C1) + (sum.W.2 - sum.B.2)/ncol(C2))
}

