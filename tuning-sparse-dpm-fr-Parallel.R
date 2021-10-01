
library(magrittr)
library(doParallel)
library(parallelDist)

no.cores <- round(detectCores() * 0.75)
cl <- makeCluster(spec = no.cores, type = 'PSOCK')
registerDoParallel(cl)

clusterExport(cl, ls(), envir = environment())

tuning.sparse.dpm.fr <- function(X, B){
   
   X <- X %>% as.matrix()
   N <- X %>% nrow()
   d <- X %>% ncol()
   
   U <- V <- matrix(0, nrow = B, ncol = N)
   
   for(i in 1:B){
      U[i,] <- sample(1:N, N, replace = TRUE)
      V[i,] <- sample(1:N, N, replace = TRUE)
   }
   
   p <- as.numeric(parDist(X, upper = TRUE, 
                           diag = TRUE, 
                           threads = no.cores))
   
   # i <- low
   low <- min(p[which(p!=0)])
   high <- max(diag(X %*% t(X)))
   incr <- (high - low)/999
   
   lambda.seq <- seq(low, high, by = incr)
   M <- matrix(0, nrow = 1000, ncol = d)
   
   clusterExport(cl, c('X','U','V','s','B','N','d','M',
                       'sparse.dpm.fr.1'),
                 envir = environment())
   
   for (s in 1:d) {
      print(s)
      
      parSapply(cl, lambda.seq, function(i){
         v <- rep(0,B)
         
         for (j in 1:B){
            res1 <- sparse.dpm.fr.1(X[U[j,],], s, i)
            res2 <- sparse.dpm.fr.1(X[V[j,],], s, i)
            
            for (h in 1:N) {
               for (g in 1:N) {
                  fU <- fV <- 0
                  
                  if (length(which(U[j,] == h)) > 0 &&
                      length(which(U[j,] == g)) > 0) {
                     if (res1$Z[min(which(U[j,] == h))] == 
                         res1$Z[min(which(U[j,] == g))] ) {
                        
                        fU <- 1
                     }
                  }
                  
                  if (length(which(V[j,] == h)) > 0 &&
                      length(which(V[j,] == g)) > 0) {
                     if (res2$Z[min(which(V[j,] == h))] == 
                         res2$Z[min(which(V[j,] == g))] ) {
                        
                        fV <- 1
                     }
                  }
                  
                  v[j] <- v[j] + abs(fU - fV)
               }
            }
            v[j] <- v[j] / N^2
         }
         
         M[which(lambda.seq == i),s] <- mean(v)
      })
   }
   
   T.ind <- which.min(M)
   
   lambda <- lambda.seq[T.ind %% 1000]
   S <- ceiling(T.ind / 1000)
   
   stopCluster(cl)
   gc()
   
   return(list(lambda, S))
}
