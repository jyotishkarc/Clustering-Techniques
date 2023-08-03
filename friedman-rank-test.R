
library(dplyr)

res_DPM <- read_excel("D:/All Downloads/DP-MoM Results (1).xlsx")

# rank.all.long <- rank.all %>% gather(key = "methods", value = "dataset", 
#                                      V2, V3, V4, V5, V6, V7, V8, V9)

friedman <- function(df){
   rank.mat <- df %>% apply(1, function(vals) rank(desc(vals))) %>% t() %>% as.data.frame()
   # colnames(rank.mat) <- paste0('V',1:9)
   
   N <- nrow(rank.mat)
   K <- ncol(rank.mat)
   
   Q <- 12*N/(K*(K+1)) * sum((colMeans(rank.mat) - (K+1)/2)^2)
   
   return(list(test.stat = Q, pval = 1 - pchisq(Q, df = K-1)))
}

df.temp <- res_DPM
colnames(df.temp) <- paste0('V',1:9)

# df.temp %>% 
#    select(-c(V5,V8,V9)) %>% 
#    apply(1, function(vals) rank(desc(vals))) %>% t() %>% as.data.frame() -> df.temp.1

friedman(df.temp)
friedman(df.temp %>% select(-c(V9)))
friedman(df.temp %>% select (-c(V8,V9)))
friedman(df.temp %>% select (-c(V5,V8,V9)))
