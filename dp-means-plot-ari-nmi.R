dpm.any.data <- function(X){
   
   var.X <- sum(apply(X[,-3], 2, var))
   lambda.X <- seq(var.X - 9*0.3, var.X + 9*0.3, by = 0.3)
   
   dpm.X <- lapply(lambda.X, function(x) {
      DP.means(X[,-3], x, X[,3])})
   
   cluster.X <- ari.X <- nmi.X <- c()
   
   for(i in 1:length(lambda.X)){
      print(i)
      cluster.X[i] <- dpm.X[[i]]$`Number-of-Clusters`
      ari.X[i] <- dpm.X[[i]]$ARI
      nmi.X[i] <- dpm.X[[i]]$NMI
   }
   
   plot(lambda.X, nmi.X, type = "o", ylim = c(min(nmi.X)-0.1,max(nmi.X)+1),
        lwd=2, main="Some Dataset : Plot for lambda vs ARI and NMI value", 
        cex.main = 1.5, xaxt="none", yaxt="none", 
        xlab = " ", ylab = " ", col = "blue")
   lines(lambda.X, ari.X, type = "o", lwd = 2, col = "red")
   
   text(lambda.X, nmi.X, labels = cluster.X, cex = 1, pos = 3, font = 2)
   
   axis(1, at = seq(var.X-2.7,var.X+2.7, by = 0.9), font = 2)
   axis(2, at=seq(min(nmi.X)-0.1,max(nmi.X)+0.1,0.05), las = 2, font = 2)
   
   abline(h=seq(min(nmi.X)-0.1,max(nmi.X)+0.1,0.05),v=lambda.X, 
          lty=3, col="gray")
   abline(v=var.X, lty=2, lwd = 3, col="gray")
   
   mtext(side=2, line=2.2, "ARI / NMI value\n", col="black", font=2,cex=1)
   mtext(side=1, line=3, "\nlambda (Penalty parameter)", 
         col="black", font=2,cex=1)
   
   legend(1.8, 0.8, legend=c("ARI","NMI"),
          col=c("red", "blue"), lty=1, cex=1.2, lwd = 3, text.font = 17)
   
   return(list(cluster.X, ari.X, nmi.X))
}
