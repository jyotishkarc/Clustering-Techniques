library(doParallel)

no.cores = round(detectCores() * 0.75)
cl = makeCluster(spec = no.cores, type = 'PSOCK')
registerDoParallel(cl)

clusterExport(cl, ls(), envir = environment())

SKFR.perm.upd <- function(X, k, B){
   
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   X.boot <- list()
   
   for(b in 1:B){
      X.boot[[b]] <- apply(X, 2, 
                           function(column){sample(column, 
                                                   length(column), 
                                                   replace = FALSE)})
   }
   
   gap.stat <- c()
   # initial.mu <- t(initial(t(X), k))
   
   for(s in 1:d){
      print(s)
      
      initial.mu <- t(initial(t(X), k))
      O.stat.X <- O.stat.upd(X, k, s, initial.mu)
      print(O.stat.X)
      
      # O.stat.X.boot <- sapply(1:B, function(b) O.stat.upd(X.boot[[b]],k,s,initial.mu))
      
      clusterExport(cl, c('X','k','s','initial.mu'), envir = environment())
      O.stat.X.boot <- parSapply(cl, 1:B, function(b){
                                          O.stat.upd(X.boot[[b]],k,s,initial.mu)})
      
      gap.stat[s] <- log(O.stat.X) - mean(log(O.stat.X.boot))
   }
   
   # stopCluster(cl)
   # gc()
   
   print(gap.stat)
   return(which.max(gap.stat))
}



O.stat.upd <- function(X, k, s, initial.mu){
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   # X.bar <- matrix(rep(colMeans(X), n), n, d, byrow = TRUE)
   
   result <- sparse.km.fr.1(X, k, s, initial.mu, ground = NULL, tolerance = 1e-3)
   L=result$features
   mu=matrix(0,n,d)
   for(i in 1:n)
      mu[i,L]=colMeans(X)[L]
   
   return(sum((X-mu)^2)-result$obj)
}

