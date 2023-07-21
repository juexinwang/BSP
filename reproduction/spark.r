##########################################################################################################################
# Analysis using SPARK in the comparison
#
# Ref: https://xzhoulab.github.io/SPARK/
# https://xzhoulab.github.io/SPARK/02_SPARK_Example/
#
# Usage:
# Rscript spark.r /home/wangjue/svgGranularity/By_Cell_Numbers/Sim_data_Cells_500.csv
##########################################################################################################################

library('SPARK')
library(pryr)
gc()
args = commandArgs(trailingOnly=TRUE)
# da<-read.csv('../By_Cell_Numbers/Sim_data_Cells_500.csv')
# passing from input
da<-read.csv(args[1])
expr <- t(da[,3:dim(da)[2]])
info<-data.frame(x=da[,"x"],
                 y=da[,"y"],
                 total_counts=apply(expr,2,sum))
rownames(info)<-rownames(da)
colnames(expr)<-rownames(da)
## debug
start_time <- Sys.time()
spark <- CreateSPARKObject(counts=expr, 
                             location=info[,1:2],
                             percentage = 0.1, 
                             min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
## Estimating Parameter Under Null
spark <- spark.vc(spark, 
                   covariates = NULL, 
                   lib_size = spark@lib_size, 
                   num_core = 1,
                   verbose = F)
## Calculating pval
spark <- spark.test(spark, 
                     check_positive = T, 
                     verbose = F)
## debug
end_time <- Sys.time()
print(mem_used())
print(end_time - start_time)
## Results
#head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])