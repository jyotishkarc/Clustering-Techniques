
start.time <- proc.time()

library(doParallel)

grid.search <- function(X, s, gt)
{
  no.cores = round(detectCores()*0.75)
  clus = makeCluster(spec = no.cores, type = 'PSOCK')
  registerDoParallel(clus)
  
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
  # dist.mat <- sort(dist.mat)
  low <- min(dist.mat)
  high <- max(mat2)
  # res <- list()
  # i <- low
  c.int <- (high - low) / 999
  # t <- 1
  # j <- c()
  # C <- c()
  # ob=c()
  # nmi <- c()
  print("HI")
  
  seq.lambda <- seq(low, high, by = c.int)
  
  clusterExport(clus, c('X','s','gt','sparse.dpm.fr.1','asg.func'), 
                envir = environment())
  
  result <- parSapply(clus, seq.lambda, function(lambda.val){
    res <- sparse.dpm.fr.1(X, s, lambda.val, gt)
    #res <- DP.means(X, i, gt)
    #res <- sparse.dp.fr(X, s, i, gt)
    return(res$NMI)
  })
  
  exec.time <- proc.time() - start.time
  print(exec.time)
  
  stopCluster(clus)
  
  return(list("Maximum NMI" = max(result), 
              "Corresponding lambda" = seq.lambda[which.max(result)]))
  
  
  # while(TRUE)
  # {
  #   
  #   #res[[t]] <- DP.means(X, i, gt)
  #   #res[[t]] <- sparse.dp.fr(X, s, i, gt)
  #   res[[t]] <- sparse.dpm.fr.1(X,s,i,gt)
  #   print(i)
  #   C[t] <- res[[t]]$C
  #   nmi[t] <- res[[t]]$NMI
  #   ob[t]=res[[t]]$obj
  #   j[t] <- i
  #   if(C[t]==1)
  #     break
  #   t <- t + 1
  #   i <- i + c.int
  # }
  # print(j[which.max(nmi)])
  # print(max(nmi))
  
  
  
  # t=1
  # v=c()
  # j1=c()
  # moja=0
  # fun=c()
  # for(i in 1:length(C))
  # {
  #   if(C[i] != moja)
  #   {
  #     v[t]=C[i]
  #     fun[t]=nmi[i]
  #     moja=v[t]
  #     j1[t]=j[i]
  #     t=t+1
  #   }
  # }
  # 
  # u=which(C<=ceiling(sqrt(N)))
  # C1=C[u]
  # j1=j[u]
  # ob1=ob[u]
  # nmi1=nmi[u]
  # # m <- (v[-1] - v[-length(v)])/(j1[-1]-j1[-length(j1)])
  # # angle <- c()
  # # for(i in 1:(length(m)-1))
  # # {
  # #   angle[i] <- atan(-abs((m[i+1] - m[i]) / (m[i+1] * m[i] + 1)))
  # #   #angle[i] <- atan(angle[i])
  # #   # if(angle[i] < 0)
  # #   #   {angle[i] <- angle[i] + pi}
  # # }
  # u=diff(ob1)
  # print(C1[which.max(u)])
  # print(nmi1[which.max(u)])
  #return(list("nmi"=nmi,"angle"=angle,"lambda"=j1,"clus"=v))
}
