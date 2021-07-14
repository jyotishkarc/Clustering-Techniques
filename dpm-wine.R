v.wine <- sum(apply(wine.data.1, 2, var))
lambda.wine <- seq(v.wine - 9*0.5,v.wine + 9*0.5, by = 0.5)

dpm.wine <- lapply(lambda.wine, function(x) {
   DP.means.proto(wine.data.1, x, gt.wine)})

C.wine <- ari.wine <- nmi.wine <- c()

for(i in 1:19){
   C.wine[i] <- dpm.wine[[i]]$`Number-of-Clusters`
   ari.wine[i] <- dpm.wine[[i]]$ARI
   nmi.wine[i] <- dpm.wine[[i]]$NMI
}

plot(lambda.wine, nmi.wine, type = "o", ylim = c(0.45,0.8), lwd=2, 
     main="winebean (Small) Dataset : Plot for lambda vs ARI and NMI value", 
     cex.main = 1.5, xaxt="none", yaxt="none", 
     xlab = " ", ylab = " ", col = "blue")
lines(lambda.wine, ari.wine, type = "o", lwd = 2, col = "red")

text(lambda.wine, nmi.wine, labels = C.wine, cex = 1, pos = 3, font = 2)
axis(1, at = seq(v.wine-4.5,v.wine+4.5, by = 1.5), font = 2)
axis(2, at=seq(0.45,0.8,0.05), las = 2, font = 2)
abline(h=seq(0.45,0.8,0.05),v=lambda.wine, lty=3, col="gray")
abline(v=v.wine, lty=2, lwd = 3, col="gray")
mtext(side=2, line=2.2, "ARI / NMI value\n", col="black", font=2,cex=1)
mtext(side=1, line=3, "\nlambda (Penalty parameter)", 
      col="black", font=2,cex=1)
legend(10.35, 0.507, legend=c("ARI","NMI"),
       col=c("red", "blue"), lty=1, cex=1.2, lwd = 3, text.font = 17)
