#Ref: https://bioconductor.org/packages/release/bioc/html/nnSVG.html
#https://bioconductor.org/packages/release/bioc/vignettes/nnSVG/inst/doc/nnSVG.html
#
# Usage:
# Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene2000_Cell2000.csv

library(SpatialExperiment)
library(STexampleData)
library(scran)
library(nnSVG)
gc()
args = commandArgs(trailingOnly=TRUE)
# passing from input
# da<-read.csv('~/Time_Mem_Cost_Test_v2/usage/Gene2000_Cell2000.csv')
da<-read.csv(args[1])
X <- da[, c("x", "y")]
expr <- t(da[,3:dim(da)[2]])
spe1 <- SpatialExperiment(
    assay = expr, 
    colData = X, 
    spatialCoordsNames = c("x", "y"))
assayNames(spe1)<-"counts"
spe1 <- computeLibraryFactors(spe1)
spe1 <- logNormCounts(spe1)
#debug
start_time <- Sys.time()
spe1 <- nnSVG(spe1)
## debug
end_time <- Sys.time()
print(mem_used())
print(end_time - start_time)