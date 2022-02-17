dpm.mom <- function(X, lambda, eps, L, eta, tmax, ground = NULL, tol = 1e-03)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- floor(n / L)
  rem <- n %% L
  B <- list()
  for(i in 1 : L)
  {
    if(i <= rem)
      B[[i]] <- X[((i - 1) * (K + 1) + 1) : (i * (K + 1)), ]
    else
      B[[i]] <- X[((i - 1) * K + 1 + rem) : (i * K + rem), ]
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
        gsum[i] <- gsum[i] + sum((G[[t]][i, ])^2) 
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
    t <- t+1
    #print(C)
  }
  if(is.null(ground) == F)
  {
    nmi.clus <- aricode::NMI(ground, Z)
    ari.clus <- aricode::ARI(ground, Z)
    # print(C)
    # print(nmi.clus)
    # print(ari.clus)
    return(list("Z" = Z, "C" = C, "NMI" = nmi.clus, "ARI" = ari.clus, 
                "objective" = f.new))
  }
  return(list("Z" = Z, "C" = C, "objective" = f.new))
}
