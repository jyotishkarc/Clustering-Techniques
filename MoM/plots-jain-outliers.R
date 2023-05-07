
library(dplyr)
library(ggplot2)
library(reshape2)
library(multipanelfigure)

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

add.outliers <- function(dataset, no.of.outliers){
   N <- nrow(dataset)
   d <- ncol(dataset)
   
   ranges <- dataset %>% apply(2, range)
   
   out.mat <- matrix(0, nrow = no.of.outliers, ncol = d)
   
   for(j in 1:d){
      out.mat[,j] <- runif(no.of.outliers, min = ranges[1,j], max = ranges[2,j])
   }
   
   return(rbind(dataset,out.mat))
}



if(FALSE){

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



if(TRUE){
   
   jain.out.20 <- add.outliers(jain[,-3],20) %>% as.data.frame()
   jain.out.40 <- add.outliers(jain[,-3],40) %>% as.data.frame()
   jain.out.60 <- add.outliers(jain[,-3],60) %>% as.data.frame()
   jain.out.80 <- add.outliers(jain[,-3],80) %>% as.data.frame()
   
   plot.jain.00 <- jain[,-3] %>% as.data.frame() %>%
      cbind("legend" = Jain_Clusters$Jain %>% as.numeric() %>% c(rep(7,0)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,35)) + 
      ggtitle("Original Dataset") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none')
   
   
   plot.jain.20 <- jain.out.20 %>%
      cbind("legend" = Jain_Clusters$Jain..20.outliers. %>% 
                           as.numeric() %>% c(rep(7,20)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,35)) + 
      ggtitle("Outlier count = 20") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none')
   
   
   plot.jain.40 <- jain.out.40 %>%
      cbind("legend" = Jain_Clusters$Jain..40.outliers. %>% 
               as.numeric() %>% c(rep(7,40)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('#E69F00','red','blue','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,35)) + 
      ggtitle("Outlier count = 40") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none')
   
   
   plot.jain.60 <- jain.out.60 %>%
      cbind("legend" = Jain_Clusters$Jain..60.outliers. %>% 
               as.numeric() %>% c(rep(7,60)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,35)) +
      ggtitle("Outlier count = 60") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none')
   
   
   plot.jain.80 <- jain.out.80 %>%
      cbind("legend" = Jain_Clusters$Jain..80.outliers. %>% 
               as.numeric() %>% c(rep(7,80)) %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend)) + 
      scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
      theme_minimal() +
      xlim(c(0,45)) +
      ylim(c(0,35)) + 
      ggtitle("Outlier count = 80") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none')
   
   gridExtra::grid.arrange(plot.jain.20, plot.jain.40, plot.jain.60, plot.jain.80,
                           nrow = 2, ncol = 2)
   
   # figure1 <- multi_panel_figure(columns = 2, rows = 2, panel_label_type = "none")
   # 
   # figure1 %<>%
   #    fill_panel(plot.jain.20, column = 1, row = 1) %<>%
   #    fill_panel(plot.jain.40, column = 2, row = 1) %<>%
   #    fill_panel(plot.jain.60, column = 1, row = 2) %<>%
   #    fill_panel(plot.jain.80, column = 2, row = 2)
   # 
   # figure1
}






