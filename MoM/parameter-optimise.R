parameter.optimise <- function(X, ground)
{
  distances <- dist(X)
  d <- c()
  for(i in 1 : nrow(X))
    d[i] <- sum(X[i, ] ^ 2)
  low <- 0
  # for(i in 1 : nrow(X))
  #   low <- low + sum((X[i, ] - colMeans(X)) ^ 2)
  low <- sum(distances) / (n * (n - 1))
  low <- min(distances)
  high <- max(d)
  temp <- max(low,high)
  temp2 <- min(low,high)
  low <- temp2
  high <- temp
  c <- (high - low) / 100
  lambda <- seq(low,high,c)
  nmi.clus <- c()
  ari.clus <- c()
  val[[i]] <- list()
  # print(lambda)
  # return()
  for(i in 1 : length(lambda))
  {
    val[[i]] <- optimise.partitions(X, lambda[i], ground)
    nmi.clus[i] <- val[[i]][[4]]
    ari.clus[i] <- val[[i]][[5]]
    print("YES")
    print(i)
  }
  j <- which.max(nmi.clus)
  r <- which.max(ari.clus)
  return(list(lambda[j], val[[j]][[3]], val[[j]][[2]], nmi.clus[j], ari.clus[r]))
}

optimise.partitions <- function(X, lambda, ground)
{
  n <- nrow(X)
  qt <- floor(n / 2)
  #lambda <- (ncol(X) + 1) * sum(dist(X)) / (n * (n - 1))
  #print(lambda)
  # B <- seq(3, (qt - ((qt + 1) %% 2)), 2)
  b <- floor(2 * nrow(X)/5) - floor(2 * nrow(X)/5) %% 2 - 1
  B <- seq(b,(qt - ((qt + 1) %% 2)),2)
  val <- list()
  nmi.clus <- c()
  ari.clus <- c()
  for(i in 1 : length(B))
  {
    val[[i]] <- dpm.mom(X, lambda, 0.05, B[i], 1, 10, ground)
    nmi.clus[i] <- val[[i]]$NMI
    ari.clus[i] <- val[[i]]$ARI
    if(val[[i]]$C == 1 | val[[i]]$C >= (n / 3))
      break
  }
  j <- which.max(nmi.clus)
  r <- which.max(ari.clus)
  return(list(lambda, B[j], val[[j]]$C, nmi.clus[j], ari.clus[r]))
}

detect.outlier <- function(X)
{
  #med <- spatial.median(X, maxiter = 20000)
  #med <- TukeyMedian(X)
  # dist <- c()
  # for(i in 1 : nrow(X))
  # {
  #   dist[i] <- sum(X[i, ] - med) ^ 2
  # }
  outlier <- 0
  mu <- colMeans(X)
  sigma <- sqrt((1 / nrow(X)) * sum(X ^ 2) - sum(mu ^ 2))
  Z <- c()
  for(i in 1 : nrow(X))
  {
    Z[i] <- sqrt(sum((X[i, ] - mu) ^ 2)) / sigma
    if(abs(Z[i]) > 2)
      outlier <- outlier + 1
  }
  #print(Z)
  print(outlier)
  return((outlier / nrow(X)))
}
