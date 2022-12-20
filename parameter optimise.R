parameter.optimise <- function(X, n.0, ground)
{
  distances <- matrix(0, nrow(X), nrow(X))
  for(i in 1 : nrow(X))
  {
    for(j in i : nrow(X))
    {
      distances[i, j] <- sum((X[i, ] - X[j, ]) ^ 2)
      distances[j, i] <- distances[i, j]
    }
  }
  distances <- as.vector(distances)
  distances <- distances[which(distances != 0)]
  d <- c()
  d = (as.matrix(dist(X))) ^ 2
  low <- min(d[d > 0])
  high = max(d)
  # for(i in 1 : nrow(X))
  #   low <- low + sum((X[i, ] - colMeans(X)) ^ 2)
  #low <- sum(distances) / (n * (n - 1))
  #low <- min(distances)
  #high <- max(d)
  #temp <- max(low,high)
  #temp2 <- min(low,high)
  #low <- temp2
  #high <- temp
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
    val[[i]] <- optimise.partitions(X, lambda[i], n.0, ground)
    nmi.clus[i] <- val[[i]][[4]]
    ari.clus[i] <- val[[i]][[5]]
    print("YES")
    print(i)
  }
  print(ari.clus)
  j <- which.max(nmi.clus)
  r <- which.max(ari.clus)
  c <- c / 10
  if(r == 11)
    lambda.second <- seq(lambda[r - 1], lambda[r], c)
  else if(r == 1)
    lambda.second <- seq(lambda[r], lambda[r + 1], c)
  else
    lambda.second <- seq(lambda[r - 1], lambda[r + 1], c)
  val.second <- list()
  nmi.clus.second <- c()
  ari.clus.second <- c()
  for(i in 1 : length(lambda.second))
  {
    val.second[[i]] <- optimise.partitions(X, lambda.second[i], n.0, ground)
    nmi.clus.second[i] <- val.second[[i]][[4]]
    ari.clus.second[i] <- val.second[[i]][[5]]
    print("YES")
    print(i)
  }
  print(ari.clus.second)
  j <- which.max(nmi.clus.second)
  r <- which.max(ari.clus.second)
  c <- c / 10
  print(j)
  if(r == 21)
    lambda.third <- seq(lambda.second[r - 1], lambda.second[r], c)
  else if(r == 1)
    lambda.third <- seq(lambda.second[r], lambda.second[r + 1], c)
  else
    lambda.third <- seq(lambda.second[r - 1], lambda.second[r + 1], c)
  val.third <- list()
  nmi.clus.third <- c()
  ari.clus.third <- c()
  for(i in 1 : length(lambda.third))
  {
    val.third[[i]] <- optimise.partitions(X, lambda.third[i], n.0, ground)
    nmi.clus.third[i] <- val.third[[i]][[4]]
    ari.clus.third[i] <- val.third[[i]][[5]]
    print("YES")
    print(i)
  }
  print(lambda.third)
  print(ari.clus.third)
  j <- which.max(nmi.clus.third)
  r <- which.max(ari.clus.third)
  return(list(lambda.third[r], val.third[[r]][[3]], val.third[[r]][[2]], nmi.clus.third[r], ari.clus.third[r]))
}

optimise.partitions <- function(X, lambda, n.0, ground)
{
  n <- nrow(X)
  qt <- floor(n / 2)
  #lambda <- (ncol(X) + 1) * sum(dist(X)) / (n * (n - 1))
  #print(lambda)
  B <- seq(3, (qt - ((qt + 1) %% 2)), 2)
  #B <- 1500
  val <- list()
  nmi.clus <- c()
  ari.clus <- c()
  for(i in 1 : length(B))
  {
    val[[i]] <- dp.mom(X, lambda, 1, B[i], 1000, 10, n.0, ground)
    nmi.clus[i] <- val[[i]]$NMI
    ari.clus[i] <- val[[i]]$ARI
    if(val[[i]]$C == 1 | val[[i]]$C >= (n / 3))
      break
  }
  #print(ari.clus)
  j <- which.max(nmi.clus)
  r <- which.max(ari.clus)
  plot(X[, 1], X[, 2], col = val[[r]]$Z)
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