#############################################################
# Generate ablation test for single granlarity
# Note: BSP uses paired granularities
#############################################################

InputPath <- "../Data/Single_Granular_Test/Var/"
InputFiles <- list.files(InputPath)
OutputPath <- "../Data/Single_Granular_Test/P_value/"
for(InputFile in InputFiles){
  InputData <- read.csv(paste0(InputPath, InputFile))
  OutputData <- InputData[,1]
  for(i in 2:ncol(InputData)){
    Fitted_K <- suppressWarnings(MASS::fitdistr(InputData[1001:10000, i], "chi-squared",
                                                start=list(df=mean(InputData[1001:10000, i]))))
    P_values <- 1 - pchisq(InputData[, i], df = Fitted_K$estimate)
    OutputData <- cbind(OutputData, P_values)
  }
  colnames(OutputData) <- colnames(InputData)
  OutputData <- as.data.frame(OutputData)
  write.csv(OutputData, row.names = FALSE,
            paste0(OutputPath, gsub("\\.csv", "_mergedPvalue\\.csv", InputFile)))
}
