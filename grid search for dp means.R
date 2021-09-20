
grid.search <- function(X, s, gt)
{
  N <- nrow(X)
  d <- ncol(X)
  mat <- matrix(0, N, N)
  t <- 1
  mat=dist(X)
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
  print(low)
  print(high)
  i <- high
  c <- (high - low) / 1000
  t <- 1
  j <- c()
  C <- c()
  ob=c()
  nmi <- c()
  sil=c()
  print("HI")
  while(i>=low)
  {
    print("yay")
    #res[[t]] <- DP.means(X, i, gt)
    res[[t]] <- sparse.dpm.fr.2(X, s, i, gt)
    #res[[t]] <- sparse.dpm.fr.1(X,s,i,gt)
    print(i)
    Z=res[[t]]$Z
    C[t] <- res[[t]]$C
    nmi[t] <- res[[t]]$NMI
    ob[t]=res[[t]]$obj
    j[t] <- i
    
    print(t)
    if(C[t]>sqrt(N))
      break
    t <- t + 1
    i <- i - c
  }
  print("Bye")
  print(j[which.max(nmi)])#dss
  print(max(nmi))#bgdbsb
  return(max(nmi))
  
  #return(list("nmi"=nmi,"angle"=angle,"lambda"=j1,"clus"=v))
}
