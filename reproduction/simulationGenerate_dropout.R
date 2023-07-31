#############################################################
# Generate 3D simulations with dropout issues
#############################################################

InputDir <- "..\\Data\\3D_Standard_Sim\\"
OutputDir <- "..\\Data\\3D_Standard_Sim_Dropout\\"

InputFiles <- list.files(InputDir)
for(InputFile in InputFiles){
  InputData <- read.csv(paste0(InputDir, InputFile))
  for(i in 1:3){
    InputExp_Dropout <- as.numeric(as.matrix(InputData[,-c(1:3)]))
    InputExp_Exc <- sample(1:length(InputExp_Dropout), size = length(InputExp_Dropout) * i / 10)
    InputExp_Dropout[InputExp_Exc] <- 0
    OutputExp <- as.data.frame(matrix(InputExp_Dropout, nrow = nrow(InputData)))
    OutputData <- cbind(InputData[,1:3], OutputExp)
    colnames(OutputData) <- colnames(InputData)
    write.csv(OutputData, row.names = FALSE, 
              paste0(OutputDir, paste0(c(strsplit(InputFile, "_")[[1]][1:3], paste0("Dropout", i*10), strsplit(InputFile, "_")[[1]][4]), collapse = "_")))
  }
}

