
library(dplyr)
library(ggplot2)
library(reshape2)

if(FALSE){
   df.jain <- data.frame("outliers" = seq(0,80,20),
                         "DP-MoM" = c(0.8999924, 0.9409352, 0.9002468, 0.8624698, 0.8747635),
                         "DP-Means" = c(0.324108, 0.3729507, 0.3486739, 0.3428117, 0.3784562),
                         "PAM" = c(0.260679, 0.260679, 0.2552051, 0.2771709, 0.2717709),
                         "K-Medians" = c(0.1931851, 0.4338798, 0.2444003, 0.4133365, 0.2390691),
                         "SKM" = c(0.08471628, 0.1049909, 0.08471628, 0.08939477, 0.105689),
                         "KM++" = c(0.3180938, 0.3483907, 0.324108, 0.3483907, 0.3729507),
                         "Kernel-KM" = c(0.427809, 0.300472, 0.492140, 0.441711, 0.462983))
   
   df.jain.long <- df.jain %>% melt(id = "outliers")
   
   df.jain.long %>% 
      ggplot(aes(x = outliers)) +
      geom_line(aes(y = value, color = variable), size = 0.8) +
      geom_point(aes(y = value, color = variable), shape = "diamond", size = 2) +
      theme_minimal() +
      ylim(c(0,1))
}


if(TRUE){

   Jain_Simulation_4[ ,-1] %>%
      cbind("legend" = Jain_Ground[,-1] %>% as.numeric() %>% c(rep(3,80)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,45)) + 
      ggtitle("Before clustering") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"))
   
   
   Jain_Simulation_4[ ,-1] %>%
      cbind("legend" = Jain_Sim4_Clustering[1:373,-1] %>% 
                           as.numeric() %>% c(rep(3,80)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('blue','#E69F00','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,45)) +
      ggtitle("After implementing DP-MoM Algorithm") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"))
   
   
   # Jain_Simulation_4[ ,-1] %>%
   #    # cbind("assignment" = Jain_Sim4_Clustering[,-1] %>% as.numeric() %>% as.factor()) %>%
   #    ggplot() +
   #    geom_point(aes(x = V1, y = V2)) + 
   #    # scale_color_discrete() +
   #    theme_minimal()
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
}


