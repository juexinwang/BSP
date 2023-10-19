library(ggplot2)
library(scales)
# raw figure
set.seed(1)
Figure_Raw <- png::readPNG(".\\..\\..\\Data\\Main_Figure\\rsz_1raw.png")
#Figure_Raw <- as.data.frame(Figure_Raw[,,1])
Figure <- data.frame(x = rep(1:ncol(Figure_Raw), rep(nrow(Figure_Raw), ncol(Figure_Raw))),
                     y = rep(c((nrow(Figure_Raw) + 1) - 1:nrow(Figure_Raw)), ncol(Figure_Raw)), 
                     value = as.numeric(unlist(as.data.frame(Figure_Raw))))
Spiked_Index <- which(Figure$value>0)
Figure$value[Spiked_Index] <- rnorm(length(Spiked_Index), mean = 3, sd = 1)
Figure$value[-Spiked_Index] <- rnorm(nrow(Figure)-length(Spiked_Index), mean = 1, sd = 1)
Figure$value[Figure$value>4] <- 4
Figure$value[Figure$value<0] <- 0

# permuted figure
Figure_Perm <- Figure
Figure_Perm$value <-sample(Figure_Perm$value, replace = FALSE)

# generate outputs
write.csv(Figure, row.names = FALSE,
          ".\\..\\..\\Data\\Main_Figure\\SVG.csv")
write.csv(Figure_Perm, row.names = FALSE,
          ".\\..\\..\\Data\\Main_Figure\\NULL.csv")


# read outputs
rm(list = ls())
Figure <- read.csv(".\\..\\..\\Data\\Main_Figure\\SVG_Mean.csv")
Figure_Perm <- read.csv(".\\..\\..\\Data\\Main_Figure\\NULL_Mean.csv")

# plot
png(file=".\\..\\..\\Outputs\\Manuscript\\SVG_Raw.png",
    width=6, height=6, units = "in", res = 600)
ggplot(Figure,
       aes(x = x, y = y, color = value)) + 
  geom_point() + 
  scale_color_gradient2(low = rgb(0,1,0),
                        high = rgb(1,0,1),
                        mid = rgb(245,245,220, maxColorValue = 255),
                        midpoint = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()


png(file=".\\..\\..\\Outputs\\Manuscript\\SVG_Patch3.png",
    width=6, height=6, units = "in", res = 600)
ggplot(Figure,
       aes(x = x, y = y, color = PatchMean_3)) + 
  geom_point() + 
  scale_color_gradient2(low = rgb(0,1,0),
                        high = rgb(1,0,1),
                        mid = rgb(245,245,220, maxColorValue = 255),
                        midpoint = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

png(file=".\\..\\..\\Outputs\\Manuscript\\SVG_Patch10.png",
    width=6, height=6, units = "in", res = 600)
ggplot(Figure,
       aes(x = x, y = y, color = PatchMean_10)) + 
  geom_point() + 
  scale_color_gradient2(low = rgb(0,1,0),
                        high = rgb(1,0,1),
                        mid = rgb(245,245,220, maxColorValue = 255),
                        midpoint = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()



png(file=".\\..\\..\\Outputs\\Manuscript\\NonSVG.png",
    width=6, height=6, units = "in", res = 600)
ggplot(Figure_Perm,
       aes(x = x, y = y, color = value)) + 
  geom_point() + 
  scale_color_gradient2(low = rgb(0,1,0),
                        high = rgb(1,0,1),
                        mid = rgb(245,245,220, maxColorValue = 255),
                        midpoint = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

png(file=".\\..\\..\\Outputs\\Manuscript\\NonSVG_Patch3.png",
    width=6, height=6, units = "in", res = 600)
ggplot(Figure_Perm,
       aes(x = x, y = y, color = PatchMean_3)) + 
  geom_point() + 
  scale_color_gradient2(low = rgb(0,1,0),
                        high = rgb(1,0,1),
                        mid = rgb(245,245,220, maxColorValue = 255),
                        midpoint = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

png(file=".\\..\\..\\Outputs\\Manuscript\\NonSVG_Patch10.png",
    width=6, height=6, units = "in", res = 600)
ggplot(Figure_Perm,
       aes(x = x, y = y, color = PatchMean_10)) + 
  geom_point() + 
  scale_color_gradient2(low = rgb(0,1,0),
                        high = rgb(1,0,1),
                        mid = rgb(245,245,220, maxColorValue = 255),
                        midpoint = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()
