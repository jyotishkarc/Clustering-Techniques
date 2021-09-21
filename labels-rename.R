labels.rename <- function(X){
   
   X <- as.matrix(X)
   
   if (length(setdiff(unique(X[,1]), 1:length(unique(X[,1])))) == 0) {
      return(list("TRAIN" = X.train, "TEST" = X.test))
   }
   
   originial.labels <- X[,ncol(X)] %>% as.character()
   new.label.names <- 1 : length(unique(originial.labels))
   
   X[,ncol(X)] <- new.label.names[as.factor(originial.labels)]
   
   return(X)
}

