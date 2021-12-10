

library(dplyr)
library(stringr)
library(readxl)

path <- "E:/Jyotishka/Data/UCR-Results/UCR/"

# path <- "C:/Users/JYOTISHKA/Dropbox/Robust_Energy_based_Classifier/Gof_Project/Results/Real/UCR/"

files <- list.files(path = path)
files.df <- as.data.frame(files, "filenames" = files)

# all.UCR <- read_excel("E:/Jyotishka/UCR-datasets.xlsx", col_names = FALSE)
all.UCR <- read_excel("~/UCR-datasets.xlsx", col_names = FALSE)

colnames(all.UCR) <- c('datasets','length')
N <- nrow(all.UCR)

M.orig <- read.csv("D:/My Documents/Real Data Analysis_ Three Databases - UCR.csv")
# M.orig <- read.csv("E:/Jyotishka/Real Data Analysis_ Three Databases - UCR.csv")

M <- M.orig[ ,-c(2:8)]

df.table <- as.data.frame(matrix(NA, nrow = N, ncol = 15))

colnames(df.table) <- H <- c("Dataset", "Length", 
                             "delta0.sin","delta0.sin.comp",
                             "delta2.sin", "delta2.sin.comp", 
                             "GLMNET", "RF", "RP","SVMlin", "SVMRBF",
                             "Nnet", "1NN", "SAVG","SVMRBF.new")

mV.which <- which(str_detect(files, "majorityVoting"))
proj.which <- which(str_detect(files, "projavg"))

Z <- list()

for (k in 1:N) {
   print(k)
   
   Z[[k]] <- files.df %>% filter(str_detect(files, all.UCR$datasets[k]))
   
   if (prod(dim(Z[[k]])) != 0) {
      
      mV <- Z[[k]] %>% filter(str_detect(files, "majorityVoting"))
      popular <- Z[[k]] %>% filter(str_detect(files, "popularclassifiers"))
      SVMRBF <- Z[[k]] %>% filter(str_detect(files, "SVMRBF"))
      SAVG <- Z[[k]] %>% filter(str_detect(files, "SAVG"))
      
      
      if (prod(dim(mV)) != 0) {
         temp.mV <- read.csv(paste(path, as.character(mV), sep = ""))
         delta.mean <- colMeans(temp.mV)
         delta.se <- apply(temp.mV, 2, function(val) sciplot::se(val))
      }
      else {delta.mean <- delta.se <- rep(NA, 4)}
      
      
      if (prod(dim(popular)) != 0) {
         temp.popular <- read.csv(paste(path, as.character(popular), sep = ""))
         temp.popular <- temp.popular[,-1]
         popular.mean <- colMeans(temp.popular)
         
         RF.pos <- which.min(popular.mean[5:8]) + 4
         popular.mean <- c(popular.mean[1:4], min(popular.mean[5:8]))
         
         popular.se <- apply(temp.popular, 2, function(val) sciplot::se(val))
         popular.se <- c(popular.se[1:4], popular.se[RF.pos])
         
         popular.mean <- c(popular.mean[1], popular.mean[5], popular.mean[2],
                           popular.mean[3], popular.mean[4], NA, NA)
         popular.se <- c(popular.se[1], popular.se[5], popular.se[2],
                         popular.se[3], popular.se[4], NA, NA)
         
         bin <- str_detect(M[,1], fixed(all.UCR$datasets[k], 
                                        ignore_case = TRUE))
         if (sum(bin) > 0) {
            pos <- which(bin)
            popular.mean <- M[pos,-1]
            popular.se <- rep(NA, 7)
         }
      }
      else {bin <- str_detect(M[,1], fixed(all.UCR$datasets[k], 
                                           ignore_case = TRUE))
      if (sum(bin) > 0) {
         pos <- which(bin)
         popular.mean <- M[pos,-1]
         popular.se <- rep(NA, 7)
      }
      else {popular.mean <- popular.se <- rep(NA, 7)}
      }
      
      
      if (prod(dim(SAVG)) != 0) {
         temp.SAVG <- read.csv(paste(path, as.character(SAVG), sep = ""))
         savg.mean <- colMeans(temp.SAVG)
         savg.se <- apply(temp.SAVG, 2, function(val) sciplot::se(val))
      }
      else {savg.value <- savg.se <- NA}
      
      
      if (prod(dim(SVMRBF)) != 0) {
         temp.SVMRBF <- read.csv(paste(path, as.character(SVMRBF), sep = ""))
         temp.SVMRBF <- as.matrix(temp.SVMRBF[,2])
         svmrbf.mean <- colMeans(temp.SVMRBF)
         svmrbf.se <- apply(temp.SVMRBF, 2, function(val) sciplot::se(val))
      }
      else {svmrbf.value <- svmrbf.se <- NA}
      
      df.table[k, ] <- c(all.UCR$datasets[k], all.UCR$length[k],
                         delta.mean, popular.mean, savg.mean, svmrbf.mean)
   }
   else {
      df.table[k, ] <- c(all.UCR$datasets[k], all.UCR$length[k],
                         rep(NA, 13))
   }
}


df.table$SVMRBF <- df.table$SVMRBF.new
df.table <- df.table[,-15]

for(j in 2:14){
   df.table[,j] <- round(as.numeric(df.table[,j]), 5)
}

View(df.table)


pref.mat <- matrix(0, N, 10)
B <- H[c(3:4,6:7,9:14)]

for (i in 1:N) {
   if(is.na(sum(df.table[i,c(3:4,6:7,9:14)])) == FALSE){
      R <- df.table[i,c(3:4,6:7,9:14)] %>%
         as.numeric() %>%
         rank(ties.method = "first")
      
      v <- c()
      for(j in 1:10){
         v[j] <- B[which(R == j)]
      }
      
      pref.mat[i,] <- v
   }
   else pref.mat[i,] <- rep(NA,10)
}

pref.mat <- cbind(all.UCR$datasets, pref.mat)
pref.mat.no.NA <- pref.mat %>% na.omit() %>% as.data.frame()

View(pref.mat.no.NA)


