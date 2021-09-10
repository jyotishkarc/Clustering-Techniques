grid.search <- function(X, s, gt)
{
  N <- nrow(X)
  d <- ncol(X)
  mat <- matrix(0, N, N)
  t <- 1
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      mat[i, j] <- sum((X[i, ] - X[j, ]) ^ 2)
    }
  }
  mat2=c()
  for(l in 1:N)
  {
    mat2[l]=sum(X[l,]^2)
  }
  mat=as.numeric(mat)
  dist.mat <- mat[mat != 0]
  dist.mat <- sort(dist.mat)
  low <- min(dist.mat)
  high <- max(mat2)
  res <- list()
  i <- low
  c <- (high - low) / 1000
  t <- 1
  j <- c()
  C <- c()
  nmi <- c()
  print("HI")
  while(TRUE)
  {
    
    #res[[t]] <- DP.means(X, i, gt)
    res[[t]] <- sparse.dp.fr(X, s, i, gt)
    #res[[t]] <- sparse.dpm.fr.1(X,s,i,gt)
    print(i)
    C[t] <- res[[t]]$C
    nmi[t] <- res[[t]]$NMI
    j[t] <- i
    if(C[t]==1)
      break
    t <- t + 1
    i <- i + c
  }
  m <- C[-1] - C[-length(C)]
  m <- m / c
  angle <- c()
  for(i in 1:(length(m)-1))
  {
    angle[i] <- -abs((m[i+1] - m[i]) / (m[i+1] * m[i] + 1))
    #angle[i] <- atan(angle[i])
    # if(angle[i] < 0)
    #   {angle[i] <- angle[i] + pi}
  }
  print(j[which.max(nmi)])
  print(max(nmi))
  return(nmi)
}