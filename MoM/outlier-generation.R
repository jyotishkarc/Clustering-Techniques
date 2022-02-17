
N <- 300

X1 <- mvtnorm::rmvnorm(N, mean = rep(0,2), sigma = diag(0.1,2))
X2 <- mvtnorm::rmvnorm(N, mean = rep(1,2), sigma = diag(0.1,2))
X3 <- mvtnorm::rmvnorm(N, mean = rep(2,2), sigma = diag(0.1,2))

out <- matrix(c(runif(3*N/4, min = min(c(X1[,1],X2[,1],X3[,1])),
                      max = max(c(X1[,1],X2[,1],X3[,1]))),
                runif(3*N/4, min = min(c(X1[,2],X2[,2],X3[,2])),
                      max = max(c(X1[,2],X2[,2],X3[,2])))),
               nrow = 3*N/4, ncol = 2)
 
# out1 <- matrix(c(runif(N/4, min = min(X1[,1]), max = max(X1[,1])),
#                  runif(N/4, min = min(X1[,2]), max = max(X1[,2]))), 
#                nrow = N/4, ncol = 2)
# 
# out2 <- matrix(c(runif(N/4, min = min(X2[,1]), max = max(X2[,1])),
#                  runif(N/4, min = min(X2[,2]), max = max(X2[,2]))), 
#                nrow = N/4, ncol = 2)
# 
# out3 <- matrix(c(runif(N/4, min = min(X3[,1]), max = max(X3[,1])),
#                  runif(N/4, min = min(X3[,2]), max = max(X3[,2]))), 
#                nrow = N/4, ncol = 2)


par(mfrow = c(2,1))

X <- rbind(X1,X2,X3)
X.out <- rbind(X1,X2,X3,out)

plot(X)
plot(X.out)


# plot(rbind(X1,X2,X3,out1,out2,out3))
