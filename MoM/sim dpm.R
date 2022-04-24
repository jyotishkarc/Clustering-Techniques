library(SimDesign)

simulate = function(k, out.lim, sigma)
{
  cd = seq(-1, 1, 2 / (k - 1))
  X = matrix(0, 650, 2)
  centroid = matrix(0, k^2, 2)
  l = 1
  Z = c()
  for(i in 1:k)
  {
    for(j in 1:k){
      centroid[l, ] = c(cd[i],cd[j])
      l = l+1}
  }
  for(i in 1:500)
  {
    s = sample(1:(k^2), 1)
    X[i, ] = rmvnorm(1, centroid[s, ], diag(sigma^2, 2))
    Z[i] = s
  }
  for(i in 501:650)
    X[i, ] <- replicate(1, runif(2,-out.lim, out.lim))
  return(list(X, Z))
}
