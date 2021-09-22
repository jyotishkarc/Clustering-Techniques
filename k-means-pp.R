## Author : JYOTISHKA RAY CHOUDHURY

km.pp <- function(data, k, ground = NULL){ # Input transposed data matrix
  
  data <- t(as.matrix(data))
  KM.PP <- km(data, k, initial(data,k))
  final.Z <- KM.PP[[2]]
  final.centroids <- KM.PP[[4]]
  final.cluster <- which(final.Z == 1) - k * 0:(ncol(data)-1)
  
  if(is.null(ground) == FALSE){
    ground <- as.numeric(as.matrix(ground))
    
    #ARI.clus <- aricode::ARI(final.cluster, ground)
    NMI.clus <- aricode::NMI(final.cluster, ground)
    
    # print(ARI.clus)
    # print(NMI.clus)
    
    cat(#"ARI =",ARI.clus,
        "\nNMI =",NMI.clus,"\n")
    
    return(list(final.cluster, 
                t(final.centroids),
                #"ARI" = ARI.clus, 
                "NMI" = NMI.clus))
  }
  
  return(list(final.cluster, t(final.centroids)))
}



###############
initial <- function(data,k){
  dm <- as.matrix(data)
  m <- nrow(dm)
  n <- ncol(dm)
  centroid.mat <- matrix(0,nrow = m,ncol = k)
  centroid.mat[,1] <- dm[,sample(1:n , 1)]
  
  for (h in 2:k) {
    dist.mat <- matrix(0 , nrow = n , ncol = h-1)
    
    for (i in 1:n) {
      for (j in 1:h-1) {
        dist.mat[i,j] <- distance.sq(dm[,i] , centroid.mat[,j])
      }
    }
    
    min.mat <- c()
    for(q in 1:n){
      min.mat[q] <- min(dist.mat[q,])
    }
    
    selection.prob <- min.mat / sum(min.mat)
    ss <- sample(1:n , 1 , prob = selection.prob)
    centroid.mat[,h] <- dm[,ss]
  }
  
  return(centroid.mat)
}


###############
km <- function(data , k , initial){
  dm <- as.matrix(data)
  m <- nrow(dm)
  n <- ncol(dm)
  centroid <- initial
  Z <- matrix(0, nrow = k , ncol = n)
  
  for (i in 1:n) {
    dist.mat <- rep(0,k)
    for (j in 1:k) {
      dist.mat[j] <-  distance.sq(dm[,i] , centroid[,j])
    }
    Z[which.min(dist.mat) , i] <- 1
  }
  
  new.centroid <- matrix(0, nrow = m , ncol = k)
  for (i in 1:k) {
    d <- which(Z[i,]==1)
    qq <- as.matrix(dm[,d])
    new.centroid[,i] <- rowMeans(qq)
  }
  
  count <- 0
  tolerance <- 1E-8
  while (distance.sq(centroid , new.centroid) > tolerance) {
    centroid <- new.centroid
    
    for (i in 1:n) {
      dist.mat <- rep(0,k)
      for (j in 1:k) {
        dist.mat[j] <- distance.sq(dm[,i] , centroid[,j])
      }
      
      Z[,i] <- rep(0,k)
      Z[which.min(dist.mat) , i] <- 1
    }
    
    for (i in 1:k) {
      d <- which(Z[i,]==1)
      qq <- as.matrix(dm[,d])
      new.centroid[,i] <- rowMeans(qq)
    }
    
    count <- count + 1
  }
  
  result <- list(rowSums(Z), Z, count, new.centroid)
  return(result)
}


################
distance.sq <- function(x,y){
  return(sum((x-y)*(x-y)))
}




################ Unbalance Dataset
if(FALSE) {
  
  library("ggplot2")
  original.plot <- ggplot(unbalance , aes(V1,V2)) + 
    geom_point(size = 2) + theme(legend.position = "none")
  
  print(original.plot)
  
  fg <- c()
  for (k in 1:6500) {
    fg[k] <- grDevices::rainbow(8)[W[[1]][k]]
  }
  
  km.pp.plot <- ggplot(rbind(unbalance , t(W[[4]])) , 
                       aes(V1,V2 , color = c(fg , rep("#000000",8)))) +
    geom_point(size = 2) + theme(legend.position = "none")
  
  print(km.pp.plot)
  
}

