library(foreach)
library(ICSNP)

no.cores <- round(detectCores() * 0.83)
cl <- makeCluster(spec = no.cores, type = 'PSOCK')
registerDoParallel(cl)
parameter.optimise.par <- function(X, ground){
  
  start.time <- proc.time()
  n <- nrow(X)
  distances <- matrix(0, nrow(X), nrow(X))
  for(i in 1 : nrow(X)){
    for(j in i : nrow(X)){
      distances[i, j] <- sum((X[i, ] - X[j, ]) ^ 2)
      distances[j, i] <- distances[i, j]
    }
  }
  med <- spatial.median(X)
  d <- c()
  for(i in 1 : nrow(X))
    d[i] <- sum((X[i, ] - med) ^ 2)
  
  # H <- hist(d)
  #high <- H$breaks[which.max(H$counts) + 1]
  low <- min(as.numeric(distances[1:(n*(6/8)), 1:(n*(6/8))])[as.numeric(distances[1:(n*(6/8)), 1:(n*(6/8))]) != 0])
  high <- max(as.numeric(distances[1:(n*(6/8)), 1:(n*(6/8))])[as.numeric(distances[1:(n*(6/8)), 1:(n*(6/8))]) != 0])
  #low <- min(min(distances[distances!=0]), which.max(H$counts))
  #high <- mean(d)
  #low <- min(distances[distances!=0])
  #high <- max(distances[distances!=0])
  c <- (high - low) / 10
  lambda <- seq(low,high,c)
  print(lambda)
  
  print("Grid 1 starting...")
  
  #################### 1st Grid
  
  res.1 <- foreach(i = 1:length(lambda), .combine = rbind,
                   .export = c('optimise.partitions','dpm.mom')) %dopar%
    {
      val <- optimise.partitions(X, lambda[i], ground)
      nmi.clus <- val[[4]]
      ari.clus <- val[[5]]
      
      return(c(nmi.clus,ari.clus))
    }
  
  print(res.1)
  
  j <- which.max(res.1[,1])
  r <- which.max(res.1[,2])
  c <- c / 10
  
  if(r == 11){
    lambda.second <- seq(lambda[r - 1], lambda[r], c)
  }else if(r == 1){
    lambda.second <- seq(lambda[r], lambda[r + 1], c)
  }else{
    lambda.second <- seq(lambda[r - 1], lambda[r + 1], c)
  }
  
  print("Grid 2 starting...")
  
  #################### 2nd Grid
  
  res.2 <- foreach(i = 1:length(lambda.second), .combine = rbind,
                   .export = c('optimise.partitions','dpm.mom'),
                   .packages = c('aricode')) %dopar%
    {
      val.second <- optimise.partitions(X, lambda.second[i], ground)
      nmi.clus.second <- val.second[[4]]
      ari.clus.second <- val.second[[5]]
      
      return(c(nmi.clus.second,ari.clus.second))
    }
  
  print(res.2)
  
  j <- which.max(res.2[,1])
  r <- which.max(res.2[,2])
  c <- c / 10
  
  if(r == 21 | (length(lambda.second == 11 & r == 11))){
     lambda.third <- seq(lambda.second[r - 1], lambda.second[r], c)
  }else if(r == 1){
     lambda.third <- seq(lambda.second[r], lambda.second[r + 1], c)
  }else{
     lambda.third <- seq(lambda.second[r - 1], lambda.second[r + 1], c)
  }
  
  print("Grid 3 starting...")
  
  #################### 3rd Grid
  
  res.3 <- foreach(i = 1:length(lambda.third), .combine = rbind,
                   .export = c('optimise.partitions','dpm.mom')) %dopar%
    {
      val.third <- optimise.partitions(X, lambda.third[i], ground)
      nmi.clus.third <- val.third[[4]]
      ari.clus.third <- val.third[[5]]
      
      return(c(nmi.clus.third,ari.clus.third))
    }
  
  print(res.3)
  
  j <- which.max(res.3[,1])
  r <- which.max(res.3[,2])
  
  val.final <- optimise.partitions(X, lambda.third[r], ground) %>% invisible()
  
  print(proc.time() - start.time)
  return(list(lambda.third[r], val.final[[3]], 
              val.final[[2]], val.final[[4]], val.final[[5]], val.final[[6]]))
}


optimise.partitions <- function(X, lambda, ground){
  n <- nrow(X)
  qt <- floor(n / 2)
  B <- seq(3, (qt - ((qt + 1) %% 2)), 2)
  B <- 9
  val <- list()
  nmi.clus <- c()
  ari.clus <- c()
  
  for(i in 1 : length(B)){
    val[[i]] <- dpm.mom(X, lambda, 1.05, B[i], 2, 10, ground)
    nmi.clus[i] <- val[[i]]$NMI
    ari.clus[i] <- aricode::ARI(val[[i]]$Z[1:600], ground[1:600])
    if(val[[i]]$C == 1 | val[[i]]$C >= (n / 3))
      break
  }
  j <- which.max(nmi.clus)
  r <- which.max(ari.clus)
  return(list(lambda, B[j], val[[j]]$C, nmi.clus[j], max(ari.clus), val[[r]]$Z))
}
