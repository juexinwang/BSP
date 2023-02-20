# Generate Supplementary Figure 4, Supplementary Figure 5, Supplementary Figure 6
setwd("D:\\PhD\\Lab\\Project_3_SVG\\Figures\\")

library(ggplot2)
library(gridExtra)
library(ggpubr)
InputDataPath <- ".\\Data\\Manuscript_3D_Simulation_Power\\"
InputFolders <- list.files(InputDataPath)
source(".\\Program\\Supp_3D_Simulation_PowerFunction.R")


for(InputFolder in InputFolders){
  InputDataFiles <- list.files(paste0(InputDataPath, InputFolder))
  for(InputDataFile in InputDataFiles){
    png(file=paste0(".\\Outputs\\", InputFolder, "_", InputDataFile, "supp_3D_power.png"),
        width=10, height=10, units = "in", res = 600)
    powerplot(paste0(InputDataPath, InputFolder, "\\", InputDataFile, "\\"))
    dev.off()
  }
}


