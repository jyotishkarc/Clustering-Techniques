
library(doParallel)

grid.search <- function(X, s = NULL, gt){
  
  start.time <- proc.time()
  
  no.cores = round(detectCores()*0.75)
  clus = makeCluster(spec = no.cores, type = 'PSOCK')
  registerDoParallel(clus)
  
  N <- nrow(X)
  # d <- ncol(X)
  mat <- matrix(0, N, N)
  
  for(i in 1:N){
    for(j in 1:N){
      mat[i, j] <- sum((X[i, ] - X[j, ]) ^ 2)
    }
  }
  
  mat2 <- c()
  for(l in 1:N){
    mat2[l] <- sum(X[l,]^2)
  }
  
  mat <- as.numeric(mat)
  dist.mat <- mat[mat != 0]
  
  low <- min(dist.mat)
  high <- max(mat2)
  print(c(low, high))
  c.int <- (high - low) / 99
  
  seq.lambda <- seq(low, high, by = c.int)
  
  print("HI")
  
  # source('~/R/R Codes/Clustering Techniques/dp-means-jrc.R')
  source('~/R/R Codes/Clustering Techniques/sparse-dpm-fr-1.R')
  # source('~/R/R Codes/Clustering Techniques/sparse-dpm-fr-2.R')
  
  # clusterExport(clus, c('X','s','gt','DP.means'), 
  #               envir = environment())
  
  clusterExport(clus, c('X','s','gt','sparse.dpm.fr.1','asg.func'),
                envir = environment())
  
  # clusterExport(clus, c('X','s','gt','sparse.dpm.fr.2'), 
  #               envir = environment())
  
  result <- parSapply(clus, seq.lambda, function(lambda.val){
    res <- sparse.dpm.fr.1(X, s, lambda.val, gt)
    # res <- sparse.dpm.fr.2(X, s, lambda.val, gt)
    # res <- DP.means(X, lambda.val, gt)
    
    return(res$NMI)
  })
  
  exec.time <- proc.time() - start.time
  print(exec.time)
  
  
  stopCluster(clus)
  gc()
  
  final.result <- list("Maximum NMI" = max(result), 
                       "Corresponding lambda" = seq.lambda[which.max(result)])
  
  print(final.result)
  return(final.result)
}
