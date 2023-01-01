dp.mom <- function(X, lambda, eps, L, eta, tmax, n.0, ground = NULL, tol = 1e-03)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- floor(n / L)
  rem <- n %% L
  B <- list()
  mu = colMeans(X)
  #d = as.matrix(dist(X))
  # j = 1
  # vec = as.numeric(c(1 : n))
  # ind = c(nrow(X))
  # print(class(ind))
  # for(i in 2 : L)
  #   if(i <= rem)
  #   {
  #     print(vec[-ind])
  #     ind[j] = permute(vec[-ind])[1]
  #     print(ind)
  #     for(k in 1 : K)
  #       ind[j + k] = which.min(d[ind[j], -ind])
  #     j = j + K + 1
  #   }
  #   else
  #   {
  #     ind[j] = permute(vec)[1]
  #     for(k in 1 : K - 1)
  #       ind[j + k] = which.min(d[ind[j], -ind])
  #     j = j + K
  #   }
  # dis = c()
  # for(i in 1 : nrow(X))
  #   dis[i] = sum((X[i, ] - mu) ^ 2)
  # vec = as.numeric(c(1 : n))
  # ind = 1
  # if(1 <= rem)
  #   ind = c(ind, (order(d[ind[1], -ind]))[1 : K])
  # else
  #   ind = c(ind, (order(d[ind[1], -ind]))[1 : (K - 1)])
  # for(i in 2 : L)
  # {
  #   ind = c(ind, (vec[-ind])[1])
  #   if(i <= rem)
  #     ind = c(ind, (order(d[ind[length(ind)], -ind]))[1 : K])
  #   else
  #     ind = c(ind, (order(d[ind[length(ind)], -ind]))[1 : (K - 1)])
  # }
  # indices = ind
  # indices = c(35, 34, 13, 30, 40, 38, 37,  2, 23, 36,  9, 41, 17,  8, 43, 26,  3, 24, 11, 20, 45, 14, 27, 15, 31, 39,
  #  32, 29, 19, 33, 12, 28, 21,  5 ,47, 42, 46,  1,  7, 10, 25, 18, 22,  6, 16,  4, 44)
   indices = permute(c(1 : n))
   #indices = order(dis)
  # indices = c(1,22,43,2,23,44,3,24,45,4,25,46,5,26,47,6,27,7,28,8,29,9,30,10,31,11,32,12,33,13,34,14,35,15,36,16,37,17,38,18,39,19,40,20,41,21,42)
  #print(indices)
  for(i in 1 : L)
  {
    if(i <= rem)
      B[[i]] = X[indices[((i - 1) * (K + 1) + 1) : (i * (K + 1))], ]
    else
      B[[i]] = X[indices[((i - 1) * K + 1 + rem) : (i * K + rem)], ]
  }
  C <- 1
  Z <- rep(1, n)
  centroid <- t(colMeans(X))
  t <- 1
  f.old <- lambda
  for(i in 1 : n)
  {
    f.old <- f.old + sum((X[i, ] - centroid[1, ]) ^ 2)
  }
  G <- list()
  while(t <= tmax)
  {
    G[[t]] <- matrix(0, n, p)
    f.theta <- c()
    for(i in 1 : L)
    {
      f.theta[i] <- 0
      for(j in 1 : nrow(B[[i]]))
      {
        dist <- c()
        for(l in 1 : C)
        {
          dist[l] <- sum((B[[i]][j, ] - centroid[l, ]) ^ 2)
        }
        f.theta <- f.theta + min(dist)
      }
      f.theta <- f.theta / nrow(B[[i]])
    }
    mom <- which(f.theta == median(f.theta))[1]
    grad <- matrix(0, C, p)
    dist <- matrix(0, nrow(B[[mom]]), C)
    for(j in 1 : nrow(B[[mom]]))
    {
      for(l in 1 : C)
      {
        dist[j, l] <- sum((B[[mom]][j, ] - centroid[l, ]) ^ 2)
      }
    }
    for(i in 1 : C)
    {
      for(j in 1 : nrow(B[[mom]]))
      {
        if(which.min(dist[j, ])==i)
          grad[i, ] <- grad[i, ] + 2 * (centroid[i, ] - B[[mom]][j, ])
      }
      grad[i, ] <- grad[i, ] / nrow(B[[mom]])
      # centroid[i, ] <- centroid[i, ] - eta * grad[i, ] / sqrt(eps + sum(grad[i, ] ^ 2))
      #centroid[i, ] <- centroid[i, ] - eta * grad[i, ]
    }
    G[[t]][1 : C, ] <- grad
    gsum <- rep(eps, C)
    for(i in 1 : C)
    {
      for(j in 1 : t)
      {
        gsum[i] <- gsum[i] + sum((G[[j]][i, ])^2) 
      }
      centroid[i, ] <- centroid[i, ] - eta * grad[i, ] / gsum[i]
    }
    #print(grad)
    f.new <- lambda * C
    
    for(i in 1 : n)
    {
      dist=c()
      for(j in 1 : C)
      {
        dist[j] <- sum((X[i, ] - centroid[j, ]) ^ 2)
      }
      y <- min(dist)
      if(y <= lambda)
      {
        Z[i] <- which.min(dist)
        f.new <- f.new + min(dist)
      }
      else
      {
        C <- C + 1
        Z[i] <- C
        #print(centroid)
        centroid <- matrix(c(as.vector(t(centroid)), X[i, ]), ncol = ncol(X), byrow = T)
        #print(X[i,])
        #print(centroid)
        #print(ncol(centroid))
        f.new <- f.new + lambda
      }
    }
    v <- sort(unique(Z))
    u <- Z
    for(i in 1 : length(v))
    {
      u[which(Z == v[i])] <- i
    }
    Z <- u
    C <- max(Z)
    if(f.new == 0)
    {
      t <- t + 1
      # print(t)
      next
    }
    if(abs(f.new / f.old - 1) < tol)
    {
      break
    }
    f.old <- f.new
    #print(f.new)
    t <- t+1
    #print(C)
  }
  
  counts= as.numeric(table(Z))
  cts.1 = which(counts >= 3)
  cts.2 = which(counts < 3)
  # print(counts)
  # print(cts.2)
  
  if(length(cts.1) < length(counts) & length(cts.2) < length(counts))
  {
    centroid = centroid[-cts.2, ]
    if(length(cts.1) == 1)
      centroid = t(centroid)
    # print(dim(centroid))
    #centroid = rbind(centroid, m)
    #print(dim(centroid))
    C = length(cts.1)

    for(i in 1 : n)
    {
      dist = c()
      for(j in 1 : C)
      {
        dist[j] = sum((X[i, ] - centroid[j, ]) ^ 2)
      }
      Z[i] = which.min(dist)
    }}
  #print("BALLE")
  if(is.null(ground) == F)
  {
    nmi.clus <- aricode::NMI(ground, Z[1 : (n - n.0)])
    ari.clus <- aricode::ARI(ground, Z[1 : (n - n.0)])
    # print(C)
    # print(nmi.clus)
    # print(ari.clus)
    return(list("Z" = Z, "C" = C, "NMI" = nmi.clus, "ARI" = ari.clus, 
                "objective" = f.new))
  }
  return(list("Z" = Z, "C" = C, "objective" = f.new))
}
