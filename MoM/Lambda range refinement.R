search.space <- function(X)
{
  mu <- spatial.median(X)
  D <- c()
  for(i in 1 : nrow(X))
  {
    D[i] <- sqrt(sum((X[i, ] - mu) ^ 2))
  }
  S <- D / mean(D)
  high.2 <- max((D[S < 4]) ^ 2)
  distances <- matrix(0, nrow(X), nrow(X))
  for(i in 1 : nrow(X))
  {
    for(j in i : nrow(X))
    {
      distances[i, j] <- sum((X[i, ] - X[j, ]) ^ 2)
      distances[j, i] <- distances[i, j]
    }
  }
  low.2 <- c()
  for(i in 1 : nrow(X))
    low.2[i] <- min(X[i, X[i, ] != 0])
  low.2 <- median(low.2)
  low.1 <- min(distances[distances != 0])
  high.1 <- max(distances)
  return(list(low.1, high.1, low.2, high.2))
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