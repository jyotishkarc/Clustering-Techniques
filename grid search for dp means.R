


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
  # print(low)
  # print(high)
  i <- low
  c <- (high - low) / 100
  t <- 1
  j <- c()
  C <- c()
  ob=c()
  nmi <- c()
  print("HI")
  print(low)
  while(i<=high)
  {
    print("yay")
    #res[[t]] <- DP.means(X, i, gt)
    # res[[t]] <- sparse.dpm.fr.2(X, s, i, gt)
    res[[t]] <- sparse.dpm.fr.1(X,s,i,gt)
    print(i)
    # Z=res[[t]]$Z
    # if(max(Z)==1)
    #  sil[t]=10
    # else
    #  sil[t]=index.DB(X,res[[t]]$Z)$DB
    C[t] <- res[[t]]$C
    nmi[t] <- res[[t]]$NMI
    ob[t]=res[[t]]$obj
    j[t] <- i
    if(C[t]==1)
      break
    print(t)
    t <- t + 1
    i <- i + c
  }
  print("Bye")
  
  print(j[which.max(nmi)])#dss
  print(max(nmi))#bgdbsb
  
  # print(t)
  # print(low)
  # print(i)
  # print(c)
  # return(max(nmi))
  
  return(list("nmi" = nmi, "lambda" = i, "low" = low, "high" = high, "inc" = c,
              "vec" = j[C < 2*N/3]))
  # return(list("nmi"=nmi,"angle"=angle,"lambda"=j1,"clus"=v))
}
