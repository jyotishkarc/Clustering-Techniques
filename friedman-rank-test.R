
library(dplyr)
library(readxl)

res_DPM <- read_excel("D:/All Downloads/DP-MoM Results (4).xlsx")

# rank.all.long <- rank.all %>% gather(key = "methods", value = "dataset", 
#                                      V2, V3, V4, V5, V6, V7, V8, V9)

friedman <- function(df){
   rank.mat <- df %>% 
                  apply(1, function(vals) rank(desc(vals))) %>% 
                  t() %>% 
                  as.data.frame()
   # colnames(rank.mat) <- paste0('V',1:9)
   
   N <- nrow(rank.mat)
   K <- ncol(rank.mat)
   
   Q <- 12*N/(K*(K+1)) * sum((colMeans(rank.mat) - (K+1)/2)^2)
   print((colMeans(rank.mat) - (K + 1)/2)^2)
   return(list(test.stat = Q, pval = 1 - pchisq(Q, df = K-1), rank.mat = rank.mat))
}

df.temp <- res_DPM[,-1]
colnames(df.temp) <- paste0('V',1:10)


# 1: K-means++ 
# 2: SKM
# 3: K-medians
# 4: PAM
# 5: RCC
# 6: DP-means
# 7: Kb-MoM
# 8: MOMPKM
# 9: OWL-KM
# 10: DP-MoM


friedman(df.temp[1:8, ])
friedman(df.temp[9:16, ])
friedman(df.temp %>% select(-c(V10)))
friedman(df.temp %>% select (-c(V9,V10)))
friedman(df.temp %>% select (-c(V8,V9,V10)))
