
library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)

# df.sim <- data.frame("outliers" = seq(0,60,20),
#                      "DP-MoM" = c(0.7442983, 0.717398, 0.7280627, 0.71878),
#                      "DP-Means" = c(0.5550988, 0.4502984, 0.4208794, 0.3041251),
#                      "PAM" = c(0.2641038, 0.1888746, 0.1771991, 0.2106863),
#                      "K-Medians" = c(0.2574187, 0.341772, 0.2353737, 0.2199716),
#                      "SKM" = c(0.247888, 0.1715387, 0.1947336, 0.1615273),
#                      "KM++" = c(0.4948044, 0.2919705, 0.2711893, 0.2264784),
#                      "Kernel-KM" = c(0.58701041, 0.5179163, 0.48834578, 0.5088454))

df.sim.df <- read_excel("D:/All Downloads/DP-MoM Results (2).xlsx") %>% na.omit()

# colnames(df.sim.df) <- c('datapoints','outliers',
#                               'Kb.MoM','RCC','DP.MoM','DP.Means','PAM','K.medians',
#                               'SKM','K.means','Kernel.K.means','OWL-KM')

df.sim.ari.long <- df.sim.df %>% 
                     # select (-c(datapoints)) %>% 
                     melt(id = "outliers")

df.sim.ari.long %>% 
   ggplot(aes(x = outliers)) +
   geom_line(aes(y = value, color = variable), linewidth = 0.9) +
   geom_point(aes(y = value, color = variable), shape = "square", size = 2) +
   labs(x = "No. of Outliers", y="ARI", color = "Algorithm  ") +
   theme_minimal() +
   ylim(c(0,1)) + 
   theme(title = element_text(face="bold"),
         # axis.title.x = element_text("No. of Outliers"),
         # axis.title.y = element_blank(),
         # legend.position = 'none',
         # legend.text = element_text('Algorithm'),
         legend.position = 'right')



