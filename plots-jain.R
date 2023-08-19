
library(dplyr)
library(ggplot2)
library(reshape2)
library(multipanelfigure)



jain.sim.data <- read.delim("D:/All Downloads/Jain (80 outliers).csv", sep = ",") %>% na.omit()
jain.sim.methods <- read.delim("D:/All Downloads/Clus Vec Jain (80 outliers).csv", sep = ",") %>% na.omit()

for(j in c(1,4,8,10)){
   jain.sim.methods[,j] <- case_when(jain.sim.methods[,j] == 1 ~ 2, 
                                     jain.sim.methods[,j] == 2 ~ 1)
}

for(j in c(7)){
   jain.sim.methods[,j] <- case_when(jain.sim.methods[,j] == 0 ~ 2,
                                     jain.sim.methods[,j] == 1 ~ 1)
}

colnames(jain.sim.methods) <- c("K-means++","SKM","K-medians","PAM",      
                                "RCC","DP-Means","Kb-MoM","MoMPKM","OWL-KM",   
                                "DP-MoM")

if(TRUE){
   
   (plot.jain.sim.data <- jain.sim.data[,-3] %>%
      cbind("legend" = jain.sim.data[,3] %>% as.numeric() %>% as.factor()) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_point(aes(color = legend), size = 0.7) + 
      # scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
      scale_color_manual(values=c('red','blue','darkolivegreen3')) +
      # theme_minimal() +
      theme_light() +
      xlab(" ") + ylab(" ") +
      xlim(c(0,42)) +
      ylim(c(0,42)) + 
      ggtitle("Original") +
      ggeasy::easy_center_title() +
      theme(title = element_text(face="bold"),
            legend.position = 'none'))
   
   plot.jain.sim.methods <- list()
   
   for(k in 1:ncol(jain.sim.methods)){
      
      plot.jain.sim.methods[[k]] <- jain.sim.data[,-3] %>%
         cbind("legend" = jain.sim.methods[,k] %>% 
                              as.numeric() %>% 
                              c(rep(3,80)) %>% 
                              as.factor()) %>%
         ggplot(aes(x = V1, y = V2)) +
         geom_point(aes(color = legend), size = 0.7) + 
         # scale_color_manual(values=c('#E69F00','blue','darkolivegreen3')) +
         scale_color_manual(values=c('red','blue','darkolivegreen3')) +
         # theme_minimal() +
         theme_light() +
         xlab(" ") + ylab(" ") +
         xlim(c(0,42)) +
         ylim(c(0,42)) + 
         ggtitle(paste0(colnames(jain.sim.methods)[k])) +
         ggeasy::easy_center_title() +
         theme(title = element_text(face="bold"),
               legend.position = 'none')
   }
}

(ggpubr::ggarrange(plot.jain.sim.data, plot.jain.sim.methods[[1]],
                  plot.jain.sim.methods[[3]], plot.jain.sim.methods[[4]],
                  plot.jain.sim.methods[[6]], plot.jain.sim.methods[[7]],
                  plot.jain.sim.methods[[8]], plot.jain.sim.methods[[10]],
                  ncol = 4, nrow = 2) -> plot.jain.sim.comparison.0.7)

ggsave('plot-jain-sim-comparison.png')

#######################################################################################################





























#################################################################################################


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






