#############################################################
# Generate Figure 4D
# Plot 3D plots of Discrete 3D patterns
#############################################################
library(plotly)
library(scales)
library(plot3D)
InputData <- read.csv(".\\..\\..\\Data\\discrete\\Discrete_900by10_width2_qt88_Noise0_pw1.csv")

spikedIndex <- read.csv(".\\..\\..\\Data\\discrete\\Discrete_Index_900by10_width2_qt88_Noise0_pw1.csv")

Input_Data_PlotSub <- InputData[spikedIndex$Indx,]



Createcolorbar <- function(InputValues){
  Cont_part1 <- ramp.col(n = 100, c(rgb(230, 230, 250, maxColorValue = 255), 
                                    rgb(255, 0, 0, maxColorValue = 255)))
  Cont_part2 <- ramp.col(n = floor(30 / (max(quantile(InputValues)[4] - min(InputValues), 1)) * (max(InputValues) - min(InputValues))) - 100, 
                         c(rgb(255, 0, 0, maxColorValue = 255), rgb(255, 0, 0, maxColorValue = 255)), alpha = 1.0)
  Cont_part1_trans <- stringr::str_to_upper(as.hexmode(floor(seq(0.6, 1.0, length.out =100)*256)))
  Cont_part1_rev <- sapply(1:length(Cont_part1), function(i){
    paste0(substr(Cont_part1[i], 1, 7), Cont_part1_trans[i])
  })
  return(c(Cont_part1_rev[-100], Cont_part2))
}

# 3d real data 1
my_palette1 <- Createcolorbar(InputData[,"SVG_1"])
Fig1 <- plot_ly() %>% 
  add_trace(data = InputData, 
            x = ~x, y = ~y, z = ~z, 
            marker = list(size = 3),
            mode = "markers", 
            opacity = 0.3,
            type = "scatter3d", 
            color = ~SVG_1, 
            colors = my_palette1) %>%
  add_trace(
    data = Input_Data_PlotSub,
    x = ~x,
    y = ~y,
    z = ~z, 
    marker = list(size = 3),
    mode = "markers", 
    type = "scatter3d", 
    opacity = 1.0,
    color = ~SVG_1,
    colors = my_palette1,
    showlegend = FALSE) %>%
  layout(showlegend = FALSE,
         margin = list(
           l = 20,
           r = 1,
           b = 20,
           t = 1,
           pad = 1
         ),
         scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
                                    eye = list(x = 0.75, y = -2.0, z = 1.0)),
                      aspectratio = list(x = 1, y = 1, z = 0.5),
                      xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                      yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                      zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>% 
  hide_colorbar()

htmlwidgets::saveWidget(as_widget(Fig1), ".\\..\\..\\Outputs\\Manuscript\\Manuscript_3D_DiscretePattern.html")
