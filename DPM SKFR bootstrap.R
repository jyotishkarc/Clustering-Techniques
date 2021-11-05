
library(tictoc)

SDPFR.boot <- function(X, low, high, inc, B){
  
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
  
  mat <- matrix(0, floor((high-low)/inc +1), d)
  i <- low
  t <- 1
  while(i <= high)
  {
  gap.stat <- c()
  # initial.mu <- t(initial(t(X), k))

  for(s in 1:d){
    tic()
    print(s)
    
    O.stat.X <- O.stat(X, i, s)
    print(O.stat.X)
    O.stat.X.boot <- sapply(1:B, function(b) O.stat(X.boot[[b]], i, s))
    gap.stat[s] <- log(O.stat.X) - mean(log(O.stat.X.boot))
    toc()
  }
  
  print(gap.stat)
  mat[t, ] <- gap.stat
  i <- i + inc
  t <- t + 1
  }
  T.ind <- which.max(mat)
  l <- seq(low, high, by = inc)
  lambda <- l[T.ind %% nrow(mat)]
  s <- ceiling(T.ind / nrow(mat))
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

