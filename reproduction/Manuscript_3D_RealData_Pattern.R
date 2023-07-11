## Generate Figure 5C,5D,5E and Supplementary Figure 25
## Plot 3D plots of genes in RA

library(plot3D)
library(plotly)
setwd("../")
InputFile <- ".\\Data\\Manuscript_3D_RealData_Pattern\\" # Figure 5C,5D,5E
# InputFile <- "/Users/wangjuex/projects/svgGranularity/3dst/" # Supplementary Figure 25


InputData_Loc <- read.csv(paste0(InputFile, "3dst_RA1_raw_loc.csv"))
InputData_Exp <- read.csv(paste0(InputFile, "3dst_RA1_raw_exp.csv"))
InputData <- as.data.frame(t(InputData_Exp[,-1]))
InputData <- cbind(InputData_Loc, InputData)
colnames(InputData) <- c("Cell", "X", "Y", "Z", InputData_Exp[,1])


Createcolorbar <- function(InputValues){
  Cont_part1 <- ramp.col(n = 100, c(rgb(0,1,0), rgb(1,0,1)))
  Cont_part2 <- ramp.col(n = floor(100 / (max(median(InputValues) - min(InputValues), 1)) * (max(InputValues) - min(InputValues))) - 100, c(rgb(1,0,1), rgb(1,0,1)), alpha = 1.0)
  Cont_part1_trans <- stringr::str_to_upper(as.hexmode(floor(seq(0.6, 1.0, length.out =100)*256)))
  Cont_part1_rev <- sapply(1:length(Cont_part1), function(i){
    paste0(substr(Cont_part1[i], 1, 7), Cont_part1_trans[i])
  })
  return(c(Cont_part1_rev[-100],Cont_part2))
}

# 3d real data 1
# my_palette1 <- Createcolorbar(InputData[,"MAN1A2"])
# Fig1 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~MAN1A2,
#             colors = my_palette1) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig1), ".\\Outputs\\Manuscript_3D_RealData_1.html")


# # 3d real data 2
# my_palette2 <- Createcolorbar(InputData[,"SEMA4D"])
# Fig2 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~SEMA4D,
#             colors = my_palette2) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig2), ".\\Output\\Manuscript_3D_RealData_2.html")


# # 3d real data 3
# my_palette3 <- Createcolorbar(InputData[,"RAC2"])
# Fig3 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~RAC2,
#             colors = my_palette3) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig3), ".\\Output\\Manuscript_3D_RealData_3.html")

##################################################
## Plot representative genes of four RA patterns
##################################################
# # RA cluster 1
# my_palette4 <- Createcolorbar(InputData[,"ACTB"])
# Fig4 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~ACTB,
#             colors = my_palette4) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig4), "/Users/wangjuex/projects/svgGranularity/Manuscript_3D_RealData_RA_0.html")

# # RA cluster 2
# my_palette5 <- Createcolorbar(InputData[,"ABCC3"])
# Fig5 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~ABCC3,
#             colors = my_palette5) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig5), "/Users/wangjuex/projects/svgGranularity/Manuscript_3D_RealData_RA_1.html")


# # RA cluster 3
# my_palette6 <- Createcolorbar(InputData[,"A2M"])
# Fig6 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~A2M,
#             colors = my_palette6) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig6), "/Users/wangjuex/projects/svgGranularity/Manuscript_3D_RealData_RA_2.html")


# # RA cluster 4
# my_palette7 <- Createcolorbar(InputData[,"ACP5"])
# Fig7 <- plot_ly() %>%
#   add_trace(data = InputData,
#             x = ~X, y = ~Y, z = ~Z,
#             marker = list(size = 5),
#             mode = "markers",
#             opacity = 0.8,
#             type = "scatter3d",
#             color = ~ACP5,
#             colors = my_palette6) %>%
#   layout(showlegend = FALSE,
#          margin = list(
#            l = 20,
#            r = 1,
#            b = 20,
#            t = 1,
#            pad = 1
#          ),
#          scene = list(camera = list(center = list(x = 0, y = 0, z = 0),
#                                     eye = list(x = 0.75, y = -2.0, z = 1.0)),
#                       aspectratio = list(x = 1, y = 1, z = 1.5),
#                       xaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       yaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE),
#                       zaxis = list(autotick = TRUE, ticks = '', showticklabels = FALSE)) ) %>%
#   hide_colorbar()

# htmlwidgets::saveWidget(as_widget(Fig7), "/Users/wangjuex/projects/svgGranularity/Manuscript_3D_RealData_RA_3.html")
