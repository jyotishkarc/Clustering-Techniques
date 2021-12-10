
library(magrittr)

labels.rename <- function(X){
   
   X <- as.matrix(X)
   
   if (length(setdiff(unique(X[,ncol(X)]), 1:length(unique(X[,ncol(X)])))) == 0) {
      return(X)
   }
   
   original.labels <- X[,ncol(X)] %>% as.character()
   new.label.names <- 1 : length(unique(original.labels))
   
   X[,ncol(X)] <- new.label.names[as.factor(original.labels)]
   
   return(X)
}

