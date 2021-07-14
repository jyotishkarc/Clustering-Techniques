v.iris <- sum(apply(iris[,-5], 2, var))
lambda.iris <- seq(v.iris - 9*0.3,v.iris + 9*0.3, by = 0.3)

dpm.iris <- lapply(lambda.iris, function(x) {
   DP.means.proto(iris.1, x, rep(1:3, each = 50))})

C.iris <- ari.iris <- nmi.iris <- c()

for(i in 1:19){
   C.iris[i] <- dpm.iris[[i]]$`Number-of-Clusters`
   ari.iris[i] <- dpm.iris[[i]]$ARI
   nmi.iris[i] <- dpm.iris[[i]]$NMI
}

plot(lambda.iris, nmi.clus, type = "o", ylim = c(0.45,0.8), lwd=2, 
     main="Iris Dataset : Plot for lambda vs ARI and NMI value", 
     cex.main = 1.5, xaxt="none", yaxt="none", 
     xlab = " ", ylab = " ", col = "blue")
lines(lambda.iris, ari.clus, type = "o", lwd = 2, col = "red")

text(lambda.iris, nmi.clus, labels = C.iris, cex = 1, pos = 3, font = 2)
axis(1, at = seq(v.iris-2.7,v.iris+2.7, by = 0.9), font = 2)
axis(2, at=seq(0.45,0.8,0.05), las = 2, font = 2)
abline(h=seq(0.45,0.8,0.05),v=lambda.iris, lty=3, col="gray")
abline(v=v.iris, lty=2, lwd = 3, col="gray")
mtext(side=2, line=2.2, "ARI / NMI value\n", col="black", font=2,cex=1)
mtext(side=1, line=3, "\nlambda (Penalty parameter)", 
      col="black", font=2,cex=1)
legend(1.8, 0.8, legend=c("ARI","NMI"),
       col=c("red", "blue"), lty=1, cex=1.2, lwd = 3, text.font = 17)
