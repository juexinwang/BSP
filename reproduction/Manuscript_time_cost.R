#############################################################
# Generate Figure 2C,2D, and Supplimentary Figure 3
#############################################################

library(ggplot2)
library(openxlsx)
library(ggpubr)
library(reshape2)
Col_plate <- c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF", 
               "#8491B4FF","#7E6148FF", "#F39B7FFF",  "#DC0000FF")

# Data for Figure 1 ----------------------------------------------------------------


InputData <- read.xlsx("..\\Data\\Time_Memory\\Time_Memory.xlsx", sheet = 1)
colnames(InputData)[1] <- "Var"
InputData_Long <- melt(InputData, id = "Var")
InputData_Long$variable <- factor(InputData_Long$variable, levels = c("BSP",
                                                                      "spark",
                                                                      "sparkx",
                                                                      "spatialDE",
                                                                      "Moran's I",
                                                                      "nnSVG"))
levels(InputData_Long$variable) <- c("BSP",
                                     "SPARK",
                                     "SPARK-X",
                                     "SpatialDE",
                                     "Moran's I",
                                     "nnSVG")
colnames(InputData_Long) <- c("NGenes", "Methods", "Value")
InputData_Long$Value <- log10(InputData_Long$Value)
Fig1 <- ggplot(InputData_Long, aes(x=NGenes, y=Value, group = Methods, color = Methods)) +
  geom_line(size = 1.0) +
  xlab("Number of Genes") + 
  ylab("Log10 Computational Time (Sec)") +
  ggtitle("Computational Time against Number of Genes") +
  scale_color_manual(values=Col_plate[1:6], drop = FALSE)



# Data for Figure 2 ----------------------------------------------------------------


InputData <- read.xlsx("..\\Data\\Time_Memory\\Time_Memory.xlsx", sheet = 2)
colnames(InputData)[1] <- "Var"
InputData_Long <- melt(InputData, id = "Var")
InputData_Long$variable <- factor(InputData_Long$variable, levels = c("BSP",
                                                                      "spark",
                                                                      "sparkx",
                                                                      "spatialDE",
                                                                      "Moran's I",
                                                                      "nnSVG"))
levels(InputData_Long$variable) <- c("BSP",
                                     "SPARK",
                                     "SPARK-X",
                                     "SpatialDE",
                                     "Moran's I",
                                     "nnSVG")
colnames(InputData_Long) <- c("NGenes", "Methods", "Value")
InputData_Long$Value <- log10(InputData_Long$Value)
Fig2 <- ggplot(InputData_Long, aes(x=NGenes, y=Value, group = Methods, color = Methods)) +
  geom_line(size = 1.0) +
  xlab("Number of Spots") + 
  ylab("Log10 Computational Time (Sec)") +
  ggtitle("Computational Time against Number of Spots") +
  scale_color_manual(values=Col_plate[1:6], drop = FALSE)





# Data for Figure 3 ----------------------------------------------------------------


InputData <- read.xlsx("..\\Data\\Time_Memory\\Time_Memory.xlsx", sheet = 3)
colnames(InputData)[1] <- "Var"
InputData_Long <- melt(InputData, id = "Var")
InputData_Long$variable <- factor(InputData_Long$variable, levels = c("BSP",
                                                                      "spark",
                                                                      "sparkx",
                                                                      "spatialDE",
                                                                      "Moran's I",
                                                                      "nnSVG"))
levels(InputData_Long$variable) <- c("BSP",
                                    "SPARK",
                                    "SPARK-X",
                                    "SpatialDE",
                                    "Moran's I",
                                    "nnSVG")
colnames(InputData_Long) <- c("NGenes", "Methods", "Value")
InputData_Long$Value <- log10(InputData_Long$Value)
Fig3 <- ggplot(InputData_Long, aes(x=NGenes, y=Value, group = Methods, color = Methods)) +
  geom_line(size = 1.0) +
  xlab("Number of Genes") + 
  ylab("Log10 Memory Usage (MB)") +
  ggtitle("Memory Usage against Number of Genes") +
  scale_color_manual(values=Col_plate[1:6], drop = FALSE)





# Data for Figure 4 ----------------------------------------------------------------


InputData <- read.xlsx("..\\Data\\Time_Memory\\Time_Memory.xlsx", sheet = 4)
colnames(InputData)[1] <- "Var"
InputData_Long <- melt(InputData, id = "Var")
InputData_Long$variable <- factor(InputData_Long$variable, levels = c("BSP",
                                                                      "spark",
                                                                      "sparkx",
                                                                      "spatialDE",
                                                                      "Moran's I",
                                                                      "nnSVG"))
levels(InputData_Long$variable) <- c("BSP",
                                     "SPARK",
                                     "SPARK-X",
                                     "SpatialDE",
                                     "Moran's I",
                                     "nnSVG")
colnames(InputData_Long) <- c("NGenes", "Methods", "Value")
InputData_Long$Value <- log10(InputData_Long$Value)
Fig4 <- ggplot(InputData_Long, aes(x=NGenes, y=Value, group = Methods, color = Methods)) +
  geom_line(size = 1.0) +
  xlab("Number of Spots") + 
  ylab("Log10 Memory Usage (MB)") +
  ggtitle("Memory Usage against Number of Spots") +
  scale_color_manual(values=Col_plate[1:6], drop = FALSE)





# Print Figures -----------------------------------------------------------
# Figure 2C,2D
png(file=paste0("..\\Outputs\\Time_cost.png"),
    width=9, height=4, units = "in", res = 600)
ggarrange(Fig1, Fig2, ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
dev.off()

# Supplementary Figure 3
png(file=paste0("..\\Outputs\\Memory_usage.png"),
    width=9, height=4, units = "in", res = 600)
ggarrange(Fig3, Fig4, ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
dev.off()
