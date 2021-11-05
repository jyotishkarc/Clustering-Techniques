
library(doParallel)

no.cores <- round(detectCores() * 0.75)
cl <- makeCluster(spec = no.cores, type = 'PSOCK')
registerDoParallel(cl)

# clusterExport(cl, ls(), envir = environment())


SDPFR.boot <- function(X, gt.X, B){
   
   start.time <- proc.time()[3]
   
   R <- grid.search(X, 1, gt.X)
   low <- R$low
   high <- R$high
   inc <- R$inc
   
   "Grid 1 done"
   
   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)
   
   X.boot <- list()
   
   for(b in 1:B){
      X.boot[[b]] <- apply(X, 2, 
                           function(column){sample(column, 
                                                   length(column), 
                                                   replace = TRUE)})
   }
   
   mat <- matrix(0, floor((high-low)/inc + 1), d)
   
   lambda.seq <- seq(low, high, by = inc)
   
   print("STARTING...")
   
   mat <- foreach(u = 1:length(lambda.seq), .combine = rbind,
                  .export = c('O.stat','sparse.dpm.fr.1',
                              'asg.func','distance.sq')) %dopar%
   {
      gap.stat <- c()
      
      for(s in 1:d){
         #tic()
         # start.time <- proc.time()
         print(s)
         
         O.stat.X <- O.stat(X, lambda.seq[u], s)
         print(O.stat.X)
         O.stat.X.boot <- sapply(1:B, function(b) O.stat(X.boot[[b]], 
                                                         lambda.seq[u], s))
         gap.stat[s] <- log(O.stat.X) - mean(log(O.stat.X.boot))
         # total.time
         #toc()
      }
      
      # t <- t + 1
      print(gap.stat)
      return(gap.stat)
      # mat[t, ] <- gap.stat
      # i <- i + inc
   }
   
   T.ind <- which.max(mat)
   # l <- seq(low, high, by = inc)
   lambda <- lambda.seq[T.ind %% nrow(mat)]
   s <- ceiling(T.ind / nrow(mat))
   
   print(proc.time()[3] - start.time)
   return(list(lambda,s))
}


O.stat <- function(X, lambda, s){

   X <- as.matrix(X)
   n <- nrow(X)
   d <- ncol(X)

   X.bar <- matrix(rep(colMeans(X), n), n, d, byrow = TRUE)

   result <- sparse.dpm.fr.1(X, s, lambda, ground = NULL, tolerance = 1e-3)
   L=result$features
   mu=matrix(0,n,d)
   for(i in 1:n)
      mu[i,L]=colMeans(X)[L]
   Z <- result[[1]]
   asg.vec <- result[[2]]

   return(sum((X-mu)^2)+lambda*result[[2]]-result$obj)
}

