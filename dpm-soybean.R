v.soy <- sum(apply(soybean.small, 2, var))
lambda.soy <- seq(v.soy - 9*0.5,v.soy + 9*0.5, by = 0.5)

dpm.soy <- lapply(lambda.soy, function(x) {
   DP.means.proto(soybean.small, x, gt.soy)})

C.soy <- ari.soy <- nmi.soy <- c()

for(i in 1:19){
   C.soy[i] <- dpm.soy[[i]]$`Number-of-Clusters`
   ari.soy[i] <- dpm.soy[[i]]$ARI
   nmi.soy[i] <- dpm.soy[[i]]$NMI
}

plot(lambda.soy, nmi.soy, type = "o", ylim = c(0.45,0.8), lwd=2, 
     main="Soybean (Small) Dataset : Plot for lambda vs ARI and NMI value", 
     cex.main = 1.5, xaxt="none", yaxt="none", 
     xlab = " ", ylab = " ", col = "blue")
lines(lambda.soy, ari.soy, type = "o", lwd = 2, col = "red")

text(lambda.soy, nmi.soy, labels = C.soy, cex = 1, pos = 3, font = 2)
axis(1, at = seq(v.soy-4.5,v.soy+4.5, by = 1.5), font = 2)
axis(2, at=seq(0.45,0.8,0.05), las = 2, font = 2)
abline(h=seq(0.45,0.8,0.05),v=lambda.soy, lty=3, col="gray")
abline(v=v.soy, lty=2, lwd = 3, col="gray")
mtext(side=2, line=2.2, "ARI / NMI value\n", col="black", font=2,cex=1)
mtext(side=1, line=3, "\nlambda (Penalty parameter)", 
      col="black", font=2,cex=1)
legend(10.35, 0.507, legend=c("ARI","NMI"),
       col=c("red", "blue"), lty=1, cex=1.2, lwd = 3, text.font = 17)
