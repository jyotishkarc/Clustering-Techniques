
DP.means <- function(X, lambda,epsilon){
   
   n <- nrow(X)
   K <- 1
   mu <- t(colMeans(X))
   z <- rep(1, n)
   fold<-0
   
   while(TRUE)
   {
      fnew<-0
      
      for (i in 1:n)
      {
         m <- matrix(rep(, ), 3, 4, byrow = TRUE)
         m<-m%*%t(m)
         if(min(diag(m))>lambda)
         {
            K<-K+1
            Z[i]<-K
            mu<-matrix(c(as.vector(t(m)),X[i,]),K,byrow<-T)
         }
         else
            Z[i]<-which.min(diag(m))
         
         j<-1
         while(j<=K)
         {
            mu[j,]<-colMeans(X[Z==j,])
            j<-j+1
         }
      }
      
      R<-matrix(0,n,ncol(X))
      j<-1
      while(j<=K)
      {
         R[Z==j,]<-mu[j,]
         j<-j+1
      }
      fnew<-sum(diag((X-R)%*%t(X-R)))
      if(abs(fnew-fold)<epsilon)
         break
      fold<-fnew
   }
   return(list("z"<-Z,"m"<-mu))
   
   
}

