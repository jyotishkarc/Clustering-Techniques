
library(aricode)
library(pbapply)
library(gtools)

## Following is a demo to run the following code:
## When the initial permutation/indexing is not mentioned:
## parameter.optimise(data, outliers, ground, distances(data)) 
## parameter.optimise(data, outliers, ground, 10 * distances(data))
##
## When the initial permutation/indexing is mentioned:
## parameter.optimise(data, outliers, ground, distances(data), index) 
## parameter.optimise(data, outliers, ground, 10 * distances(data), index)
## 
## data     =  dataset in matrix form with n (number of data points) rows and p (number of features) columns
## outliers =  number of noisy data points added to the dataset
## ground   =  the ground truth corresponding to the dataset 'data' excluding the 'outliers'
## index    =  the random permutation of the data points

median.0 = function(u) return((sort(u))[(length(u) + 1) / 2])

DP.MoM = function(X, lambda, eps, L, eta, n.0, ground = NULL, 
                  index = NULL, tol = 1e-03, tmax = 21)
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
    u = Z
    for(i in 1 : length(v))
    {
        u[which(Z == v[i])] = i
    }
    Z = u
    C = max(Z)
    
    if(is.null(ground) == F)
    {
        nmi.clus = aricode::NMI(ground, Z[1 : (n - n.0)])
        ari.clus = aricode::ARI(ground, Z[1 : (n - n.0)])
        return(list("Z" = Z, "C" = C, "NMI" = nmi.clus, "ARI" = ari.clus,
                    "objective" = f.new, "Time" = t, "Indices" = indices, 
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
    return(10 ^ (ceiling(2 * log10(high)) / 2 - 1))
}

parameter.optimise <- function(X, n.0, ground, eta, index = NULL)
{
    d <- c()
    d = (as.matrix(dist(X))) ^ 2
    low <- min(d[d > 0])
    high = max(d)
    
    buckets = c()
    c <- (high - low) / 10
    lambda <- seq(low,high,c) ## Vector of lambda values for first run
    print(lambda)
    nmi.clus <- c()
    ari.clus <- c()
    ob.clus <- c()
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
    j.1 <- which.max(nmi.clus)
    r.1 <- which.max(ari.clus)
    x.1 <- which.min(ob.clus)
    clus.1 <- val[[r.1]][[6]]
    buckets = c(buckets, val[[r.1]][[2]])
    
    c <- c / 10
    if(r.1 == 11) ## Selecting vector of lambda values for second run
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
    j.2 <- which.max(nmi.clus.second)
    r.2 <- which.max(ari.clus.second)
    x.2 <- which.min(ari.clus.second)
    clus.2 <- val.second[[r.2]][[6]]
    
    c <- c / 10
    print(j)
    
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
    j.3 <- which.max(nmi.clus.third)
    r.3 <- which.max(ari.clus.third)
    x.3 <- which.min(ob.clus.third)
    clus.3 <- val.third[[r.3]][[6]]
    m = max(max(ari.clus), max(ari.clus.second), max(ari.clus.third))
    if(m == max(ari.clus))
        return(list(lambda[r.1], val[[r.1]][[3]], val[[r.1]][[2]], nmi.clus[r.1], ari.clus[r.1],
                    val[[r.1]][[6]], val[[r.1]][[8]]))
    else if(m == max(ari.clus.second))
        return(list(lambda.second[r.2], val.second[[r.2]][[3]], val.second[[r.2]][[2]], 
                    nmi.clus.second[r.2], ari.clus.second[r.2], val.second[[r.2]][[6]], 
                    val.second[[r.2]][[8]]))
    else return(list(lambda.third[r.3], val.third[[r.3]][[3]], val.third[[r.3]][[2]], 
                    nmi.clus.third[r.3], ari.clus.third[r.3], val.third[[r.3]][[6]], 
                    val.third[[r.3]][[8]]))
}

optimise.partitions <- function(X, lambda, n.0, ground, eta, index = NULL)
{
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
        if(is.null(index) == TRUE)
            val[[i]] <- DP.MoM(X, lambda, 1, B[i], eta, n.0, ground)
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
    x <- which.min(ob.clus)
    return(list(lambda, B[r], val[[r]]$C, nmi.clus[r], ari.clus[r], 
                val[[r]]$Z, min(ob.clus), val[[r]]$Indices))
}
