ew.dpm.mom = function(X, lambda.k, lambda.w, eps, L, eta, n.0, T.max, ground = NULL, tol = 1e-03)
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
  
  B = list() #Partitioning of the data into buckets
  indices = permute(1 : n)
  for(i in 1 : L)
  {
    if(i <= rem)
      B[[i]] = X[indices[((i - 1) * (K + 1) + 1) : (i * (K + 1))], ]
    else
      B[[i]] = X[indices[((i - 1) * K + 1 + rem) : (i * K + rem)], ]
  }

  f.old = rep((lambda.k - lambda.w * log(p)) / n, L) #Computation of initial objective function 
  for(i in 1 : L)
  {
    for(j in 1 : nrow(B[[i]]))
    {
      f.old[i] = f.old[i] + (1 / nrow(B[[i]])) * as.numeric(t(B[[i]][j, ] - centroid[1, ]) %*% W %*% (B[[i]][j, ] - centroid[1, ]))
    }
  }
  fun.old = median(f.old)
  print(f.old)
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
      if(min(dist) > lambda.k)
      {
        centroid = rbind(centroid, X[i, ])
        C = C + 1
        Z[i] = C
        for(x in 1 : 3)
        {
          f.new = rep((lambda.k * C + lambda.w * sum(w * log(w))), L)
          if(rem > 0){
            for(r in 1 : rem)
            {
              for(s in 1 : nrow(B[[r]]))
              {
                for(k in 1 : C)
                {
                  if(Z[indices[(r - 1) * (K + 1) + s]] == k)
                    f.new[r] = f.new[r] + (1 / nrow(B[[r]])) * as.numeric(t(B[[r]][s, ] - centroid[k, ]) %*% W %*% (B[[r]][s, ] - centroid[k, ]))
                }
              }
            }}
          for(r in (rem + 1) : L)
          {
            for(s in 1 : nrow(B[[r]]))
            {
              for(k in 1 : C)
              {
                if(Z[indices[(rem) * (K + 1) + (r - rem - 1) * K + s]] == k)
                  f.new[r] = f.new[r] + (1 / nrow(B[[r]])) * as.numeric(t(B[[r]][s, ] - centroid[k, ]) %*% W %*% (B[[r]][s, ] - centroid[k, ]))
              }
            }
          }
          print("Iteration")
          print(i)
          #print(f.new)
          L.t = which(f.new == median(f.new))
          grad <- matrix(0, C, p)
          dist <- matrix(0, nrow(B[[L.t]]), C)
          for(v in 1 : nrow(B[[L.t]]))
          {
            for(l in 1 : C)
            {
              dist[v, l] <- sum(w * (B[[L.t]][v, ] - centroid[l, ]) ^ 2)
            }
          }
          D = rep(0, p)
          for(r in 1 : C)
          {
            for(s in 1 : nrow(B[[L.t]]))
            {
              if(which.min(dist[s, ])==r)
              {
                grad[r, ] <- grad[r, ] + 2 * (centroid[r, ] - B[[L.t]][s, ])
                D = D + (centroid[r, ] - B[[L.t]][s, ]) ^ 2
              }
            }
            grad[r, ] <- grad[r, ] / nrow(B[[L.t]])
            # centroid[i, ] <- centroid[i, ] - eta * grad[i, ] / sqrt(eps + sum(grad[i, ] ^ 2))
            #centroid[i, ] <- centroid[i, ] - eta * grad[i, ]
          }
          G[[t]][1 : C, ] <- grad
          gsum <- rep(eps, C)
          for(r in 1 : C)
          {
            for(s in 1 : t)
            {
              gsum[r] <- gsum[r] + sum((G[[s]][r, ])^2) 
            }
            centroid[r, ] <- centroid[r, ] - eta * grad[r, ] / gsum[r]
          }
          for(y in 1 : n)
          {
            dist.2 = c()
            for(a in 1 : C)
            {
              dist.2[a] = as.numeric(t(X[y, ] - centroid[a, ]) %*% W %*% (X[y, ] - centroid[a, ]))
            }
            Z[y] = which.min(dist.2)
          } 
        }
      }
      else
        Z[i] = which.min(dist)
    }
    print(Z)
    f.new = rep((lambda.k * C + lambda.w * sum(w * log(w))), L)
    if(rem > 0){
    for(i in 1 : rem)
    {
      for(j in 1 : nrow(B[[i]]))
      {
        for(k in 1 : C)
        {
          if(Z[indices[(i - 1) * (K + 1) + j]] == k)
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
          if(Z[indices[(rem) * (K + 1) + (i - rem - 1) * K + j]] == k)
            f.new[i] = f.new[i] + (1 / nrow(B[[i]])) * as.numeric(t(B[[i]][j, ] - centroid[k, ]) %*% W %*% (B[[i]][j, ] - centroid[k, ]))
        }
      }
    }
    #print(f.new)
    L.t = which(f.new == median(f.new))
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
        gsum[i] <- gsum[i] + sum((G[[j]][i, ])^2) 
      }
      centroid[i, ] <- centroid[i, ] - eta * grad[i, ] / gsum[i]
    }
    w = exp(-D / lambda.w)
    w = w / sum(w)
    W = diag(w)
    f.new = rep((lambda.k * C + lambda.w * sum(w * log(w))) / n, L)
    print(f.new)
    if(rem > 0){
      for(i in 1 : rem)
      {
        for(j in 1 : nrow(B[[i]]))
        {
          for(k in 1 : C)
          {
            if(Z[indices[(i - 1) * (K + 1) + j]] == k)
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
          if(Z[indices[(rem) * (K + 1) + (i - rem - 1) * K + j]] == k)
            f.new[i] = f.new[i] + (1 / nrow(B[[i]])) * as.numeric(t(B[[i]][j, ] - centroid[k, ]) %*% W %*% (B[[i]][j, ] - centroid[k, ]))
        }
      }
    }
    fun.new = median(f.new)
    print(f.new)
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
    print(t)
  }
  print("HI")
  counts= as.numeric(table(Z))
  cts.1 = which(counts > 2)
  cts.2 = which(counts <= 2)
  if(length(cts.1) < length(counts))
  {centroid = centroid[-cts.2, ]
  print(centroid)
  C = length(cts.1)
  for(i in 1 : n)
  {
    dist = c()
    for(j in 1 : C)
    {
      dist[j] = sum(w * (X[i, ] - centroid[j, ]) ^ 2)
    }
    Z[i] = which.min(dist)
  }}
  if(is.null(ground) == F)
  {
    nmi.clus <- aricode::NMI(ground, Z[1:(n-n.0)])
    ari.clus <- aricode::ARI(ground, Z[1:(n-n.0)])
    # print(C)
    # print(nmi.clus)
    # print(ari.clus)
    return(list("Z" = Z, "C" = C, "NMI" = nmi.clus, "ARI" = ari.clus, 
                "objective" = fun.new, "W" = w))
  }
  return(list("Z" = Z, "C" = C, "objective" = fun.new))
}