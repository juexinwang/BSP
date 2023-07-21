##########################################################################################################################
# Analysis using nnSVG in the comparison
#
# Ref: https://bioconductor.org/packages/release/bioc/html/nnSVG.html
# https://bioconductor.org/packages/release/bioc/vignettes/nnSVG/inst/doc/nnSVG.html
#
# Usage:
# Rscript nnSVG.r '/home/wangjue/fig1/20220415_DEFAULT-changeNoise_sim/sim_MOB_pattern4_fc4_tau80_count_power9.csv'
# Rscript nnSVG.r /home/wangjue/Time_Mem_Cost_Test_v2/usage/Gene2000_Cell2000.csv
##########################################################################################################################

library(SpatialExperiment)
library(STexampleData)
library(scran)
library(nnSVG)
library(pryr)
gc()
args = commandArgs(trailingOnly=TRUE)
# passing from input
# da<-read.csv('/home/wangjue/fig1/20220415_DEFAULT-changeNoise_sim/sim_MOB_pattern4_fc4_tau80_count_power9.csv')
infile<-paste('~/fig1/20220415_DEFAULT_sim/',args[1],sep='')
da<-read.csv(infile)
outfile<-paste('~/fig1_result/',args[1],sep='')
# outfile<-paste('~/fig1_result_default/',args[1],sep='')
X <- da[, c("x", "y")]
expr <- t(da[,3:dim(da)[2]])
print(dim(expr))
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
write.table(rowData(spe1)['padj'],file=outfile,quote=FALSE,sep=",",col.names=FALSE)
## debug
end_time <- Sys.time()
print(mem_used())
print(end_time - start_time)

