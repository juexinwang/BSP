#############################################################
# Generate Figure 2A
#############################################################

library(ggplot2)
library(scales)
library(ggpubr)
Input_data <- read.csv(".\\..\\..\\Data\\Manuscript_2D_Simulation_Pattern\\sim_MOB_pattern2_fc5_tau50_count_power1.csv")
Figure1 <- ggplot(Input_data, aes(x = x, y = y,
                               color = Input_data[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradient2(low = rgb(230, 230, 250, maxColorValue = 255), 
                        high = rgb(255, 0, 0, maxColorValue = 255), 
                        mid = rgb(255, 240, 245, maxColorValue = 255),
                        midpoint = median(Input_data[,"gene1"], 0.25),
                        breaks = unlist(lapply(list(min, median, max), function(f) f(Input_data[,"gene1"]))),
                        labels = c("Low","", "High")) +
  theme_bw() + 
  labs(color = "Relative expression") +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


Input_data2 <- read.csv(".\\..\\..\\Data\\Manuscript_2D_Simulation_Pattern\\sim_MOB_pattern3_fc5_tau50_count_power1.csv")
Figure2 <- ggplot(Input_data2, aes(x = x, y = y,
                       color = Input_data2[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradient2(low = rgb(230, 230, 250, maxColorValue = 255), 
                        high = rgb(255, 0, 0, maxColorValue = 255), 
                        mid = rgb(255, 240, 245, maxColorValue = 255),
                        midpoint = median(Input_data[,"gene1"], 0.25),
                        breaks = unlist(lapply(list(min, median, max), function(f) f(Input_data2[,"gene1"]))),
                        labels = c("Low","", "High"))+
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())



Input_data3 <- read.csv(".\\..\\..\\Data\\Manuscript_2D_Simulation_Pattern\\sim_MOB_pattern4_fc5_tau50_count_power1.csv")
Figure3 <- ggplot(Input_data3, aes(x = x, y = y,
                       color = Input_data3[,"gene1"])) +
  geom_point(size = 3) + 
  scale_color_gradient2(low = rgb(230, 230, 250, maxColorValue = 255), 
                        high = rgb(255, 0, 0, maxColorValue = 255), 
                        mid = rgb(255, 240, 245, maxColorValue = 255),
                        midpoint = median(Input_data[,"gene1"], 0.25),
                        breaks = unlist(lapply(list(min, median, max), function(f) f(Input_data3[,"gene1"]))),
                        labels = c("Low","", "High"))+
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

Final_figure <- ggarrange(Figure1, Figure2, Figure3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")


png(file=".\\..\\..\\Outputs\\Manuscript\\Manuscript_2D_Simulation_Pattern.png",width = 10,height = 3.5, units = "in",
    res = 600)
print(Final_figure)
dev.off()