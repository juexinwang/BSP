## Generate Figure 2A

library(ggplot2)
library(scales)
library(ggpubr)
setwd("D:\\PhD\\Lab\\Project_3_SVG\\Figures\\")


Input_data <- read.csv(".\\Data\\Manuscript_2D_Simulation_Pattern\\sim_MOB_pattern2_fc5_tau50_count_power1.csv")
Figure1 <- ggplot(Input_data, aes(x = x, y = y,
                               color = Input_data[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradient2(low=rgb(0,1,0), 
                        high=rgb(1,0,1), 
                        mid = rgb(245,245,220, maxColorValue = 255), 
                        midpoint = median(Input_data[,"gene1"])) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


Input_data2 <- read.csv(".\\Data\\Manuscript_2D_Simulation_Pattern\\sim_MOB_pattern3_fc5_tau50_count_power1.csv")
Figure2 <- ggplot(Input_data2, aes(x = x, y = y,
                       color = Input_data2[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradient2(low=rgb(0,1,0), 
                        high=rgb(1,0,1), 
                        mid = rgb(245,245,220, maxColorValue = 255), 
                        midpoint = median(Input_data2[,"gene1"])) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())



Input_data3 <- read.csv(".\\Data\\Manuscript_2D_Simulation_Pattern\\sim_MOB_pattern4_fc5_tau50_count_power1.csv")
Figure3 <- ggplot(Input_data3, aes(x = x, y = y,
                       color = Input_data3[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradient2(low=rgb(0,1,0), 
                        high=rgb(1,0,1), 
                        mid = rgb(245,245,220, maxColorValue = 255), 
                        midpoint = median(Input_data3[,"gene1"])) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

Final_figure <- ggarrange(Figure1, Figure2, Figure3, nrow = 1, ncol = 3)

png(file=".\\Outputs\\Manuscript_2D_Simulation_Pattern.png",width = 10,height = 3,units = "in",
    res = 600)
print(Final_figure)
dev.off()