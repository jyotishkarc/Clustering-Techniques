#### Uncomment the next 3 lines if necessary
# install.packages("aricode")
# install.packages("pbapply")
# install.packages("gtools")

library(aricode)
library(pbapply)
library(gtools)


median.0 = function(u) return((sort(u))[(length(u) + 1) / 2])


DP.MoM = function(X, lambda, eps, L, eta, n.0, ground = NULL, 
                  index = NULL, tol = 1e-03, tmax = 51)
{
  n = nrow(X)
  p = ncol(X)
  K = floor(n / L)
  rem = n %% L
  B = list()
  mu = colMeans(X)
  
  if(is.null(index) == FALSE)
    indices = index
  else
    indices = gtools::permute(c(1 : n))
  for(i in 1 : L) ## Partitioning of the data points
  {
    if(i <= rem)
      B[[i]] = X[indices[((i - 1) * (K + 1) + 1) : (i * (K + 1))], ]
    else
      B[[i]] = X[indices[((i - 1) * K + 1 + rem) : (i * K + rem)], ]
  }
  C = 1
  Z = rep(1, n)
  centroids = t(colMeans(X))
  t = 1
  f.old = lambda
  for(i in 1 : L)
  {
    dist = rep(0, L)
    for(j in 1 : nrow(B[[i]]))
    {
      dist[i] = dist[i] + sum((B[[i]][j, ] - centroids[1, ]) ^ 2)
    }
    dist[i] = dist[i] / nrow(B[[i]])
  }
  f.old = f.old + median.0(dist) ## Computation of the initial value of the objective function
  G = list()
  while(t <= tmax)
  {
    G[[t]] = matrix(0, n, p)
    temp = MoM(X, B, centroids, L)
    f = temp[[1]]
    mom = temp[[2]] ## Determining the bucket corresponding to the Median-of-Means estimator
    
    temp = AdaGrad(X, B, G, centroids, L, C, eps, eta, mom, t)
    G[[t]][1 : nrow(centroids), ] = temp[[2]]
    centroids = temp[[1]] ## Recomputation of the centroids using the AdaGrad algorithm
    
    for(i in 1 : n)
    {
      dist=c()
      for(j in 1 : C)
      {
        dist[j] = sum((X[i, ] - centroids[j, ]) ^ 2)
      }
      y = min(dist)
      if(y <= lambda) ## Checking for the need to create a new cluster
      {
        Z[i] = which.min(dist)
      }
      else
      {
        C = C + 1
        Z[i] = C
        centroids = rbind(centroids, X[i, ])
      }
    }
    
    C = length(unique(Z))
    centroids = centroids[sort(unique(Z)), ]
    v = sort(unique(Z))
    u = Z
    for(i in 1 : length(v))
    {
      u[which(Z == v[i])] = i
    }
    Z = u
    if(length(unique(Z)) == 1)
      centroids = t(centroids)
    
    for(j in 1 : n) ## Assigning data points to nearest cluster centroid
    {
      dist = c()
      for(j in 1 : C)
      {
        dist[j] = sum((X[i, ] - centroids[j, ]) ^ 2)
      }
      Z[i] = which.min(dist)
    }
    
    C = length(unique(Z))
    centroids = centroids[sort(unique(Z)), ]
    v = sort(unique(Z))
    u = Z
    for(i in 1 : length(v))
    {
      u[which(Z == v[i])] = i
    }
    Z = u
    if(length(unique(Z)) == 1)
      centroids = t(centroids)
    
    temp = MoM(X, B, centroids, L)
    f.theta = temp[[1]]
    
    f.new = median.0(f.theta) + lambda * C ## Computation of the new value of the objective function
    
    if(f.new == 0)
    {
      t = t + 1
      next
    }
    if(abs(f.new / f.old - 1) < tol) ## Checking the stopping condition
    {
      break
    }
    f.old = f.new
    t = t+1
  }
  v = sort(unique(Z))
  centroids = centroids[v,]
  u = Z
  for(i in 1 : length(v))
  {
    u[which(Z == v[i])] = i
  }
  Z = u
  C = max(Z)
  
  counts = as.numeric(table(Z))
  if(min(counts)<3)
  {
    ind.new = which(counts < 3)
    for(i in 1 : length(ind.new))
    {
      for(j in 1 : n)
      {
        if(Z[j] == ind.new[i])
        {
          dist = c()
          for(k in 1 : C)
          {
            dist[k] = sum((centroids[k, ] - X[j, ])^2)
          }
          Z[j] = which.min(dist[-ind.new])
        }
      }
    }
  }
  C = length(unique(Z))
  centroids = t(t(centroids))
  f.old = C * lambda
  
  for(i in 1 : L)
  {
    dist = rep(0, L)
    for(j in 1 : nrow(B[[i]]))
    {
      d = c()
      for(l in 1:C)
        d[l] = sum((B[[i]][j, ] - centroids[l, ]) ^ 2)
      dist[i] = dist[i] + min(d)
    }
    dist[i] = dist[i] / nrow(B[[i]])
  }
  #print(dist)
  f.old = f.old + median(dist)
  centroids = centroids[unique(Z), ]
  
  # lam.est = 0
  # distances = matrix(0, n, C)
  # for(i in 1 : n)
  # {
  #   for(j in 1 : C)
  #   {
  #     distances[i, j] = sum((X[i, ] - centroids[j, ]) ^ 2)
  #   }
  # }
  # temp = rep(0, n)
  # for(i in 1 : n)
  # {
  #   temp[i] = min(distances[i, -which.min(distances[i, ])])
  # }
  # lam.est = (lambda - (max(rowMin(distances)) + min(temp))/2)^2
  # if(C==1)
  #   lam.est = 0
  
  if(is.null(ground) == FALSE)
  {
    nmi.clus = aricode::NMI(ground, Z[1 : (n - n.0)])
    ari.clus = aricode::ARI(ground, Z[1 : (n - n.0)])
    return(list("Z" = Z, "C" = C, "NMI" = nmi.clus, "ARI" = ari.clus,
                "objective" = f.old, "Time" = t, "Indices" = indices, 
                "Centroids" = centroids))
  }
  return(list("Z" = Z, "C" = C, "objective" = f.new, "Centroids" = centroids))
}


AdaGrad = function(X, B, G, centroids, L, C, eps, eta, mom, t)
{
  n = nrow(X)
  p = ncol(X)
  K = floor(n / L)
  rem = n %% L
  #C = nrow(centroids)
  
  grad = matrix(0, C, p)
  dist = matrix(0, nrow(B[[mom]]), C)
  
  for(i in 1 : nrow(B[[mom]]))
  {
    for(j in 1 : C)
    {
      dist[i, j] = sum((B[[mom]][i, ] - centroids[j, ]) ^ 2)
    }
  }
  for(i in 1 : C)
  {
    for(j in 1 : nrow(B[[mom]]))
    {
      if((which.min(dist[j, ]))[1] == i)
        grad[i, ] = grad[i, ] + 2 * (centroids[i, ] - B[[mom]][j, ])
    }
    grad[i, ] = grad[i, ] / nrow(B[[mom]])
  }
  G[[t]][c(1 : C), ] = grad
  gsum = rep(eps, C)
  for(i in 1 : C)
  {
    for(j in 1 : t)
    {
      gsum[i] = gsum[i] + sum((G[[j]][i, ])^2)
    }
    centroids[i, ] = centroids[i, ] - eta * grad[i, ] / sqrt(gsum[i])
  }
  return(list(centroids, grad))
}


MoM = function(X, B, centroids, L)
{
  n = nrow(X)
  p = ncol(X)
  K = floor(n / L)
  rem = n %% L
  C = nrow(centroids)
  f = c()
  for(i in 1 : L)
  {
    f[i] = 0
    for(j in 1 : nrow(B[[i]]))
    {
      dist = c()
      for(k in 1 : C)
      {
        dist[k] = sum((B[[i]][j, ] - centroids[k, ]) ^ 2)
      }
      f[i] = f[i] +  min(dist)
    }
    f[i] = f[i] / nrow(B[[i]])
  }
  mom.est = median.0(f)
  L.t = (which(f == mom.est))[1]
  return(list(f, L.t))
}


distances = function(X) ##Computation of optiumum value of learning rate eta
{
  d <- c()
  d = (as.matrix(dist(X))) ^ 2
  low <- min(d[d > 0])
  high = max(d)
  return(10 ^ (ceiling(log10(high)/2) ))
}


parameter.optimise.par <- function(X, n.0, ground, eta, index = NULL){
  d <- c()
  d = (as.matrix(dist(X))) ^ 2
  low <- min(d[d > 0])
  high = max(d)
  
  buckets = c()
  c <- (high - low) / 10
  lambda <- seq(low,high,c) ## Vector of lambda values for first run
  lambda = lambda[-c(1, length(lambda))]
  print(lambda)
  nmi.clus <- c()
  ari.clus <- c()
  ob.clus <- c()
  dist.clus <- c()
  val <- list()
  for(i in 1 : length(lambda))
  {
    if(is.null(index) == TRUE)
      val[[i]] <- optimise.partitions(X, lambda[i], n.0, ground, eta)
    else
      val[[i]] <- optimise.partitions(X, lambda[i], n.0, ground, eta, index)
    nmi.clus[i] <- val[[i]][[4]]
    ari.clus[i] <- val[[i]][[5]]
    ob.clus[i] <- val[[i]][[7]]
    print("YES")
    print(i)
  }
  print(ari.clus)
  print(ob.clus)
  j.1 <- which.max(nmi.clus)
  r.1 <- which.max(ari.clus)
  #r.1 = find_shallowest_local_minima(ob.clus)
  #r.1 = which.min((diff(ob.clus)))
  x.1 <- which.min(ob.clus)
  clus.1 <- val[[r.1]][[6]]
  buckets = c(buckets, val[[r.1]][[2]])
  
  c <- c / 10
  if(r.1 == 9) ## Selecting vector of lambda values for second run
    lambda.second <- seq(lambda[r.1 - 1], lambda[r.1], c)
  else if(r.1 == 1)
    lambda.second <- seq(lambda[r.1], lambda[r.1 + 1], c)
  else
    lambda.second <- seq(lambda[r.1 - 1], lambda[r.1 + 1], c)
  val.second <- list()
  nmi.clus.second <- c()
  ari.clus.second <- c()
  ob.clus.second <- c()
  for(i in 1 : length(lambda.second))
  {
    if(is.null(index) == TRUE)
      val.second[[i]] <- optimise.partitions(X, lambda.second[i], n.0, ground, eta)
    else
      val.second[[i]] <- optimise.partitions(X, lambda.second[i], n.0, ground, eta, index)
    nmi.clus.second[i] <- val.second[[i]][[4]]
    ari.clus.second[i] <- val.second[[i]][[5]]
    ob.clus.second[i] <- val.second[[i]][[7]]
    buckets = c(buckets, val.second[[i]][[2]])
    print("YES")
    print(i)
  }
  print(ari.clus.second)
  print(ob.clus.second)
  j.2 <- which.max(nmi.clus.second)
  r.2 <- which.max(ari.clus.second)
  #r.2 = find_shallowest_local_minima(ob.clus.second)
  #r.2 = which.min((diff(ob.clus.second)))
  x.2 <- which.min(ari.clus.second)
  clus.2 <- val.second[[r.2]][[6]]
  
  c <- c / 10
  # print(j)
  
  if(length(lambda.second) == 11) ## Selecting vector of lambda values for second run
  {
    if(r.2 == 11)
      lambda.third <- seq(lambda.second[r.2 - 1], lambda.second[r.2], c)
    else if(r.2 == 1)
      lambda.third <- seq(lambda.second[r.2], lambda.second[r.2 + 1], c)
    else
      lambda.third <- seq(lambda.second[r.2 - 1], lambda.second[r.2 + 1], c)
  }
  else
  {
    if(r.2 == 21)
      lambda.third <- seq(lambda.second[r.2 - 1], lambda.second[r.2], c)
    else if(r.2 == 1)
      lambda.third <- seq(lambda.second[r.2], lambda.second[r.2 + 1], c)
    else
      lambda.third <- seq(lambda.second[r.2 - 1], lambda.second[r.2 + 1], c)
  }
  print(min(lambda.third))
  print(max(lambda.third))
  print(min(buckets))
  print(max(buckets))
  val.third <- list()
  nmi.clus.third <- c()
  ari.clus.third <- c()
  ob.clus.third <- c()
  for(i in 1 : length(lambda.third))
  {
    if(is.null(index) == TRUE)
      val.third[[i]] <- optimise.partitions(X, lambda.third[i], n.0, ground, eta)
    else
      val.third[[i]] <- optimise.partitions(X, lambda.third[i], n.0, ground, eta, index)
    nmi.clus.third[i] <- val.third[[i]][[4]]
    ari.clus.third[i] <- val.third[[i]][[5]]
    ob.clus.third[i] <- val.third[[i]][[7]]
    print("YES")
    print(i)
  }
  print(lambda.third)
  print(ari.clus.third)
  print(ob.clus.third)
  j.3 <- which.max(nmi.clus.third)
  r.3 <- which.max(ari.clus.third)
  #r.3 = find_shallowest_local_minima(ob.clus.third)
  #r.3 = which.min((diff(ob.clus.third)))
  x.3 <- which.min(ob.clus.third)
  clus.3 <- val.third[[r.3]][[6]]
  m = max(max(ari.clus), max(ari.clus.second), max(ari.clus.third))
  if(m == max(ari.clus))
    return(list('lambda' = lambda[r.1], 'Number of Buckets' = val[[r.1]][[2]], 
                'Estimated number of clusters' = val[[r.1]][[3]], 
                'NMI' = nmi.clus[r.1], 'ARI' = ari.clus[r.1],
                'Clster assignment vector' = val[[r.1]][[6]], 
                'Index' = val[[r.1]][[8]], 'Centroids' = val[[r.1]][[9]]))
  else if(m == max(ari.clus.second))
    return(list('lambda' = lambda.second[r.2], 'Number of Buckets' = val.second[[r.2]][[2]], 
                'Estimated number of clusters' = val.second[[r.2]][[3]], 
                'NMI' = nmi.clus.second[r.2], 'ARI' = ari.clus.second[r.2], 
                'Clster assignment vector' = val.second[[r.2]][[6]], 
                'Index' = val.second[[r.2]][[8]], 'Centroids' = val.second[[r.2]][[9]]))
  else return(list('lambda' = lambda.third[r.3], 'Number of Buckets' = val.third[[r.3]][[2]], 
                   'Estimated number of clusters' = val.third[[r.3]][[3]], 
                   'NMI' = nmi.clus.third[r.3], 'ARI' = ari.clus.third[r.3], 
                   'Clster assignment vector' = val.third[[r.3]][[6]], 
                   'Index' = val.third[[r.3]][[8]], 'Centroids' = val.third[[r.3]][[9]]))
}


optimise.partitions <- function(X, lambda, n.0, ground, eta, index = NULL){
  n <- nrow(X)
  qt <- floor((n - n.0) / 3)
  B <- c(3 : qt)
  val <- list()
  nmi.clus <- c()
  ari.clus <- c()
  ob.clus <- c()
  n_iter = length(B)
  pb = txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for(i in 1 : length(B))
  {
    ind = index.select(X, B[i])
    if(is.null(index) == TRUE)
      val[[i]] <- DP.MoM(X, lambda, 1, B[i], eta, n.0, ground, ind)
    else
      val[[i]] <- DP.MoM(X, lambda, 1, B[i], eta, n.0, ground, index)
    nmi.clus[i] <- val[[i]]$NMI
    ari.clus[i] <- val[[i]]$ARI
    ob.clus[i] <- val[[i]]$objective
    setTxtProgressBar(pb, i)
    if(val[[i]]$C == 1 | val[[i]]$C >= (n / 3))
      break
  }
  j <- which.max(nmi.clus)
  r <- which.max(ari.clus)
  #r <- which.min(ob.clus)
  #r = find_first_local_minima(ob.clus)
  if(is.null(r) == TRUE)
    r = 2 + which.min(ob.clus[-c(1,2)])
  x <- which.min(ob.clus)
  return(list(lambda, B[r], val[[r]]$C, nmi.clus[r], ari.clus[r], 
              val[[r]]$Z, min(ob.clus), val[[r]]$Indices, val[[r]]$Centroids))
}

index.select = function(X, L)
{
  n = nrow(X)
  q = floor(n / L)
  r = n %% L
  
  index = rep(0, n)
  k = 1
  for(i in 1 : L)
  {
    if(i <= r){
      if(i==1)
        index[((i - 1) * (q + 1) + 1) : (i * (q + 1))] = (1:n)[km.pp(X, q+1)]
      else
        index[((i - 1) * (q + 1) + 1) : (i * (q + 1))] = ((1:n)[-index[index>0]])[km.pp(X[-index[index>0],], q+1)]
    }
    else{
      if(i==1)
        index[((i - 1) * q + 1 + r) : (i * q + r)] = (1:n)[km.pp(X, q)]
      else
        index[((i - 1) * q + 1 + r) : (i * q + r)] = ((1:n)[-index[index>0]])[km.pp(X[-index[index>0],], q)]
    }
  }
  return(index)
}


km.pp = function(X, k)
{
  n = nrow(X)
  index = rep(0, k)
  index[1] = base::sample(c(1:n), 1)
  for(i in 1 : (k-1))
  {
    D = rep(0, n)
    for(j in 1 : n)
    {
      dist = rep(0, i)
      for(l in 1 : i)
      {
        dist[l] = sum((X[j, ] - X[index[l]])^2)
      }
      D[j] = min(dist)
    }
    pr = D / sum(D)
    index[i+1] = base::sample(c(1:n), 1, prob = pr)
  }
  return(index)
}