#############################################################
# Generate Figure 4A
#############################################################

library(plotly)
library(scales)
setwd("D:\\PhD\\Lab\\Project_3_SVG\\Figures\\")

for(Pattern_index in 1:3){
  Input_Folder <- ".\\Data\\Manuscript_3D_Simulation_Pattern\\Simulation_3D_Example\\"
  Input_File <- paste0("Pattern_",Pattern_index,".csv")
  Input_Data <- read.csv(paste0(Input_Folder, Input_File))
  
  SEED = Pattern_index 
  N = 225
  numSlices = 10
  width = 3.5
  SubPattern = Pattern_index
  set.seed(SEED)
  SimData <- data.frame()
  for(i in 1:numSlices){
    Coords <- spatstat.random::rpoispp(N, win = spatstat.geom::owin(c(0, 1), c(0, 1)))
    SimData <- rbind(SimData, data.frame(x = Coords$x * (sqrt(N) - 1), y = Coords$y * (sqrt(N) - 1), z = i))
  }
  coordDF <- SimData
  
  set.seed(SEED)
  if(SubPattern == "1"){
    Phi_const <- pi * round(runif(1, 0, 2)) - pi/2
    Phi <- pi*runif(20) + Phi_const
    Theta_const <- 0.5 * pi * round(runif(1, 0, 3))
    Theta <- 0.5 * pi*runif(20) + Theta_const
  }
  if(SubPattern == "2"){
    Phi_const <- pi * round(runif(1, 0, 2)) - pi/2
    Phi <- pi*runif(20) + Phi_const
    Theta_const <- pi * round(runif(1, 0, 2))
    Theta <- pi * runif(20) + Theta_const
  }                
  if(SubPattern == "3"){
    Phi <- 2 * pi * runif(20)
    Theta <- 2 * pi * runif(20)
  }
  Center_pts <- cbind(cumsum(2 * cos(Phi) * sin(Theta)) + sqrt(N)/2,
                      cumsum(2 * cos(Phi) * cos(Theta)) + sqrt(N)/2,
                      cumsum(2 * sin(Phi)) + numSlices/2)
  
  #plot(cumsum(2*sin(Theta)), cumsum(2*cos(Theta)),type = "b", xlim = c(-10,10), ylim =c(-10,10) )
  
  spikedCoords <- sapply(1:length(Phi), function(Center_pt){
    # get distance from every cell to the center of the coordinate set; assuming hot spot is centered here
    distFromCenter <- sp::spDists(
      x = coordDF %>% dplyr::select(x, y, z) %>% as.matrix(),
      y = matrix(Center_pts[Center_pt, ], ncol = 3),
      longlat = FALSE
    ) %>% as.vector()
    
    # identify coordinates in hot spot
    spikedCoord <- coordDF[distFromCenter <= (width-1),] %>% row.names() %>% as.numeric()
  })
  spikedCoords <- unique(unlist(spikedCoords))
  
  
  
  
  Input_Data_Plot <- Input_Data[,1:4]
  # trim extreme values
  Input_Data_Plot[Input_Data_Plot[,4]>quantile(Input_Data_Plot[,4],0.95),4] <- quantile(Input_Data_Plot[,4],0.95)
  Input_Data_PlotSub <- Input_Data_Plot[spikedCoords,]
  Center_pts <- as.data.frame(Center_pts)
  colnames(Center_pts) <- c("x", "y", "z")
  Center_pts <- Center_pts[which(apply(Center_pts, 1, max)<sqrt(N) & apply(Center_pts, 1, min)>0),]
  my_palette <- c(rgb(0,1,0),
                  rgb(245,245,220, maxColorValue = 255),
                  rgb(1,0,1))
  
  Fig <- plot_ly() %>% 
    add_trace(data = Input_Data_Plot, 
              x = ~x, y = ~y, z = ~z, 
              marker = list(size = 5),
              mode = "markers", 
              type = "scatter3d", 
              opacity = 0.6,
              color = ~SVG_1, 
              colors = my_palette) %>%
    add_trace(
      data = Input_Data_PlotSub,
      x = ~x,
      y = ~y,
      z = ~z, 
      marker = list(size = 5),
      mode = "markers", 
      type = "scatter3d", 
      opacity = 1.0,
      color = ~SVG_1,
      colors = my_palette,
      showlegend = FALSE)%>%
    add_trace(
      data = Center_pts,
      x = ~x,
      y = ~y,
      z = ~z, 
      type="scatter3d", 
      opacity = 0.3,
      mode="lines",
      line = list(color = rgb(1,0,1), width = 5),
      showlegend = FALSE)%>%
    layout(showlegend = FALSE,
           margin = list(
             l = 20,
             r = 1,
             b = 20,
             t = 1,
             pad = 1
           ),
           scene = list(camera = list(eye = list(x = 1.25, y = -1.25, z = 1.25)),
                        xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                        yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
                        zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>% 
    hide_colorbar()
  
  
  
  htmlwidgets::saveWidget(as_widget(Fig), 
                          paste0(".\\Outputs\\Manuscript_3D_Pattern_",Pattern_index,".html"))
}

