DP.means.new = function(X, lambda, ground = NULL, tmax,tol = 1e-3)
{
  n = nrow(X)
  p = ncol(X)
  Z = rep(1, n)
  centroid = colMeans(X)
  obj.old = lambda
  for(i in 1 : n)
    obj.old = obj.old + sum((X[i, ] - centroid) ^ 2)
  C = 1
  t = 1
  centroid = t(centroid)
  
  
  while(t <= tmax)
  {
    for(i in 1 : n)
    {
      dist = c()
      for(j in 1 : C)
        dist[j] = sum((X[i, ] - centroid[j, ]) ^ 2)
      if(min(dist) > lambda)
      {
        C = C + 1
        centroid = rbind(centroid, X[i, ])
        for(k in 1 : 3)
        {
          for(i in 1 : n){
            dist = c()
          for(j in 1 : C)
            dist[j] = sum((X[i, ] - centroid[j, ]) ^ 2)
          Z[i] = which.min(dist)}
          for(i in 1 : C)
          {
            if(length(which(Z == i)) == 1)
              centroid[i] = X[which(Z==i)]
            else if(length(which(Z == i)) == 0)
              centroid[i] = rep(0, p)
            else
              centroid[i] = colMeans(X[which(Z == i), ])
          }
        }
      }
      else
        Z[i] = which.min(dist)
    }
    for(i in 1 : C)
    {
      if(length(which(Z == i)) == 1)
        centroid[i] = X[which(Z==i)]
      else if(length(which(Z == i)) == 0)
        centroid[i] = rep(0, p)
      else
        centroid[i] = colMeans(X[which(Z == i), ])
    }
    obj.new = lambda * C
    for(i in 1 : n)
      obj.new = obj.new + sum((X[i, ] - centroid[Z[j], ]) ^ 2)
    v <- sort(unique(Z))
    
    u <- Z
    for(i in 1:(length(v)))
    {
      u[which(Z==v[i])]=i
    }
    Z <- u
    C <- max(Z)
    if(abs(obj.new / obj.old - 1) < tol)
      break
    t = t + 1
    obj.old = obj.new
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
  if(is.null(ground)==F){
    nmi=NMI(Z,ground)
    ari=ARI(Z,ground)}
  print(C)
  print(nmi)
  return(list("Z" = Z , "mu" = centroid,"NMI"=nmi, 
              "C" = C, "ARI"=ari, "objective" = obj.new))
}