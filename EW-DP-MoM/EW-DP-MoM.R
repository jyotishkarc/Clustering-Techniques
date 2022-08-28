ew.dpm.mom = function(X, lambda.k, lambda.w, eps, L, eta, T.max, ground = NULL, tol = 1e-03)
{
  n = nrow(X) #Initialization of the Data
  p = ncol(X)
  K = floor(n / L)
  rem = n %% L
  w = rep(1 / p, p)
  C = 1
  Z = rep(1, n)
  W = diag(w)
  centroid = t(colMeans(X))
  G = list()
  
  B <= list() #Partitioning of the data into buckets
  indices = permute(1 : n)
  for(i in 1 : L)
  {
    if(i <= rem)
      B[[i]] = X[indices[((i - 1) * (K + 1) + 1) : (i * (K + 1))], ]
    else
      B[[i]] = X[indices[((i - 1) * K + 1 + rem) : (i * K + rem)], ]
  }
  
  f.old = rep(lambda.k - lambda.w * log(p)) #Computation of initial objective function 
  for(i in 1 : L)
  {
    for(j in 1 : nrow(B[[i]]))
    {
      f.old[i] = f.old[i] + (1 / nrow(B[[i]])) * as.numeric(t(B[[i]][j, ] - centroid[1, ]) %*% W %*% (B[[i]][j, ] - centroid[1, ]))
    }
  }
  fun.old = median(f.old)
  
  t = 1
  while(t <= T.max)
  {
    G[[t]] <- matrix(0, n, p)
    for(i in 1 : n) #Updating cluster assignments
    {
      dist = c()
      for(j in 1 : C)
      {
        dist[j] = as.numeric(t(X[i, ] - centroid[j, ]) %*% W %*% (X[i, ] - centroid[j, ]))
      }
      if(dist[j] > lambda.k)
      {
        centroid = rbind(centroid, X[i, ])
        C = C + 1
        Z[i] = C
      }
      else
        Z[i] = which.min(dist)
    }
    
    f.new = rep((lambda.k * C + lambda.w * sum(w * log(w))), L)
    if(rem > 0){
    for(i in 1 : rem)
    {
      for(j in 1 : nrow(B[[i]]))
      {
        for(k in 1 : C)
        {
          if(Z[indices[(i - 1) * (L + 1) + j]] == k)
            f.new[i] = f.new[i] + (1 / nrow(B[[i]])) * as.numeric(t(B[[i]][j, ] - centroid[k, ]) %*% W %*% (B[[i]][j, ] - centroid[k, ]))
        }
      }
    }}
    for(i in (rem + 1) : L)
    {
      for(j in 1 : nrow(B[[i]]))
      {
        for(k in 1 : C)
        {
          if(Z[indices[(rem) * (L + 1) + (i - rem - 1) * L + j]] == k)
            f.new[i] = f.new[i] + (1 / nrow(B[[i]])) * as.numeric(t(B[[i]][j, ] - centroid[k, ]) %*% W %*% (B[[i]][j, ] - centroid[k, ]))
        }
      }
    }
    L.t = which.median(f.new)
    grad <- matrix(0, C, p)
    dist <- matrix(0, nrow(B[[L.t]]), C)
    for(j in 1 : nrow(B[[L.t]]))
    {
      for(l in 1 : C)
      {
        dist[j, l] <- sum(w * (B[[L.t]][j, ] - centroid[l, ]) ^ 2)
      }
    }
    D = rep(0, p)
    for(i in 1 : C)
    {
      for(j in 1 : nrow(B[[L.t]]))
      {
        if(which.min(dist[j, ])==i)
        {
          grad[i, ] <- grad[i, ] + 2 * (centroid[i, ] - B[[L.t]][j, ])
          D = D + (centroid[i, ] - B[[L.t]][j, ]) ^ 2
        }
      }
      grad[i, ] <- grad[i, ] / nrow(B[[L.t]])
      # centroid[i, ] <- centroid[i, ] - eta * grad[i, ] / sqrt(eps + sum(grad[i, ] ^ 2))
      #centroid[i, ] <- centroid[i, ] - eta * grad[i, ]
    }
    G[[t]][1 : C, ] <- grad
    gsum <- rep(eps, C)
    for(i in 1 : C)
    {
      for(j in 1 : t)
      {
        gsum[i] <- gsum[i] + sum((G[[t]][i, ])^2) 
      }
      centroid[i, ] <- centroid[i, ] - eta * grad[i, ] / gsum[i]
    }
    w = exp(-D / lambda.w)
    w = w / sum(w)
    fun.new = median(f.new)
    v <- sort(unique(Z))
    u <- Z
    for(i in 1 : length(v))
    {
      u[which(Z == v[i])] <- i
    }
    Z <- u
    C <- max(Z)
    if(fun.new == 0)
    {
      t <- t + 1
      # print(t)
      next
    }
    if(abs(fun.new / fun.old - 1) < tol)
      break
    fun.old = fun.new
    t = t + 1
  }
  if(is.null(ground) == F)
  {
    nmi.clus <- aricode::NMI(ground, Z)
    ari.clus <- aricode::ARI(ground, Z)
    # print(C)
    # print(nmi.clus)
    # print(ari.clus)
    return(list("Z" = Z, "C" = C, "NMI" = nmi.clus, "ARI" = ari.clus, 
                "objective" = fun.new))
  }
  return(list("Z" = Z, "C" = C, "objective" = fun.new))
}