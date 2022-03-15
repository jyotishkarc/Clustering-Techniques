binomial.mom <- function(N, n, p)
{
   data <- rbinom(N, n, p)
   p.hat <- 1 - (1 - 1 / N) * var(data) / mean(data)
   n.hat <- mean(data) / p.hat
   return(list(n.hat, p.hat))
}

multinormal.mle <- function(N, mu,sigma)
{
   data <- matrix(0, N, length(mu))
   for(i in 1 : N)
      data[i, ] <- rmvnorm(1, mu, sigma)
   mu.hat <- colMeans(data)
   sigma.hat <- matrix(0, length(mu), length(mu))
   for(i in 1 : N)
   {
      sigma.hat <- sigma.hat + (data[i, ] - mu.hat) %*% t(data[i, ] - mu.hat)
   }
   sigma.hat <- sigma.hat / N
   return(list(mu.hat, sigma.hat))
}
