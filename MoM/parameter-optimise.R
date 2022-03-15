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
  c <- (high - low) / 10
  lambda <- seq(low,high,c)
  print(lambda)
  nmi.clus <- c()
  ari.clus <- c()
  val <- list()
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
  c <- c / 10
  if(j == 11)
    lambda.second <- seq(lambda[j - 1], lambda[j], c)
  else if(j == 1)
    lambda.second <- seq(lambda[j], lambda[j + 1], c)
  else
    lambda.second <- seq(lambda[j - 1], lambda[j + 1], c)
  val.second <- list()
  nmi.clus.second <- c()
  ari.clus.second <- c()
  for(i in 1 : length(lambda.second))
  {
    val.second[[i]] <- optimise.partitions(X, lambda.second[i], ground)
    nmi.clus.second[i] <- val.second[[i]][[4]]
    ari.clus.second[i] <- val.second[[i]][[5]]
    print("YES")
    print(i)
  }
  j <- which.max(nmi.clus.second)
  r <- which.max(ari.clus.second)
  c <- c / 10
  print(j)
  lambda.third <- seq(lambda.second[j - 1], lambda.second[j + 1], c)
  val.third <- list()
  nmi.clus.third <- c()
  ari.clus.third <- c()
  for(i in 1 : length(lambda.third))
  {
    val.third[[i]] <- optimise.partitions(X, lambda.third[i], ground)
    nmi.clus.third[i] <- val.third[[i]][[4]]
    ari.clus.third[i] <- val.third[[i]][[5]]
    print("YES")
    print(i)
  }
  j <- which.max(nmi.clus.third)
  r <- which.max(ari.clus.third)
  return(list(lambda.third[j], val.third[[j]][[3]], val.third[[j]][[2]], nmi.clus.third[j], ari.clus.third[r]))
}

optimise.partitions <- function(X, lambda, ground)
{
  n <- nrow(X)
  qt <- floor(n / 2)
  #lambda <- (ncol(X) + 1) * sum(dist(X)) / (n * (n - 1))
  #print(lambda)
  #B <- seq(1501, (qt - ((qt + 1) %% 2)), 2)
  B <- 1500
  val <- list()
  nmi.clus <- c()
  ari.clus <- c()
  for(i in 1 : length(B))
  {
    val[[i]] <- dpm.mom(X, lambda, 1, B[i], 1.02, 10, ground)
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

spatial <- function(X)
{
  mu <- colMeans(X)
  t <- 0
  while(t <= 100){
    s <- 0
    for(i in 1 : nrow(X))
    {
      s <- s + (1 / sqrt(sum(X[i, ] - mu)))
    }
    m <- rep(0, ncol(X))
    for(i in 1 : nrow(X))
    {
      m <- m + X[i, ] / sqrt(sum(X[i, ] - mu))
    }
    mu <- m / s
    t <- t + 1}
  return(mu)
}
