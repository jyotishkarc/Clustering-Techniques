


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
  sil=c()
  nmi <- c()
  sil=c()
  print("HI")
  print(low)
  while(i<=high)
  {
    print("yay")
    #res[[t]] <- DP.means(X, i, gt)
    #res[[t]] <- sparse.dp.fr(X, s, i, gt)
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
  #plot(j,sil,type="o")
  # print(j[which.min(sil)])
  # print(nmi[which.min(sil)])
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
  # print(length(ob))
  # C1=unique(C)
  # j1=rep(0,length(C1))
  # ob1=rep(0,length(C1))
  # moja=rep(0,length(C1))
  # for(i in 1:length(C1))
  # {
  #   j1[i]=mean(j[which(C==C1[i])])
  #   ob1[i]=ob[median(which(C==C1[i]))]
  # }
  # ob2=ob1-j1*C1
  # a=rep(0,length(C1)-3)
  # for(i in 2:length(C1)-2)
  # {
  #   model1=lm(C1[1:i]~log(j1[1:i]))
  #   model2=lm(C1[-(1:i)]~log(j1[-(1:i)]))
  #   a[i-1]=deviance(model1)+deviance(model2)
  # }
  # print(j1[which.min(a)+1])#fdsgdsd
  # a=rep(0,length(C1)-3)
  # for(i in 2:length(C1)-2)
  # {
  #   model1=lm(log(ob2[1:i])~C1[1:i])
  #   model2=lm(log(ob2[-(1:i)])~C1[-(1:i)])
  #   a[i-1]=deviance(model1)+deviance(model2)
  # }
  # print(j1[which.min(a)+1])#dsbsdbsd
  # # j1=rev(j1)
  # # C1=rev(C1)
  # # ob1=rev(ob1)
  # # ob2=rev(ob2)
  # slope=(log(ob2[-1])-log(ob2[-length(C1)]))/(C1[-1]-C1[-length(j1)])
  # angle=c()
  # print(slope)
  # #plot(C1,log(ob2),type="o")
  # for(i in 1:(length(slope)-1))
  # {
  #   if(slope[i]<=slope[i+1])
  #     angle[i]=-abs((slope[i]-slope[i+1])/(1+slope[i]*slope[i+1]))
  #   else
  #     angle[i]=0
  # }
  # # m=1
  # # for(i in 1:length(angle))
  # # {
  # #   if(angle[i]>=(pi/2) & angle[i]<angle[m])
  # #   {
  # #     m=i
  # #   }
  # # }
  # print(angle)
  # print(j1)
  # print(C1)
  # print(j1[which.min(angle)+1])
  print(j[which.max(nmi)])#dss
  print(max(nmi))#bgdbsb
  # print(t)
  # print(low)
  # print(i)
  # print(c)
  # return(max(nmi))
  return(list("nmi" = nmi, "lambda" = i, "low" = low, "high" = high, "inc" = c))
  # return(list("nmi"=nmi,"angle"=angle,"lambda"=j1,"clus"=v))
}
