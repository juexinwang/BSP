##########################################################################################################################
# Analysis using SpatialDE in the comparison
#
# http://bioconductor.org/packages/release/bioc/vignettes/spatialDE/inst/doc/spatialDE.html
# Changed for only performance usage and record
# Ref: http://bioconductor.org/packages/release/bioc/vignettes/spatialDE/inst/doc/spatialDE.html
#
# Usage:
# Rscript spatialDE.r /home/wangjue/svgGranularity/By_Cell_Numbers/Sim_data_Cells_500.csv
##########################################################################################################################

library(spatialDE)
library(pryr)
gc()
args = commandArgs(trailingOnly=TRUE)
# da<-read.csv('../By_Cell_Numbers/Sim_data_Cells_500.csv')
# passing from input
da<-read.csv(args[1])
X <- da[, c("x", "y")]
expr <- t(da[,3:dim(da)[2]])
total_counts <- colSums(expr)
da_info<-data.frame(X$x, X$y, total_counts)
## debug
start_time <- Sys.time()
norm_expr <- stabilize(data.matrix(expr))
resid_expr <- regress_out(norm_expr, sample_info = da_info)
results <- spatialDE::run(resid_expr, coordinates = X)
## debug
end_time <- Sys.time()
print(mem_used())
print(end_time - start_time)

