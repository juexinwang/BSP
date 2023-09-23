#############################################################
# Calculate spatial stratified heterogeneity
#
# Ref: Wang, Jin-Feng, Tong-Lin Zhang, and Bo-Jie Fu. "A measure of spatial stratified heterogeneity." Ecological indicators 67 (2016): 250-256
#############################################################

# Loading data ------------------------------------------------------------

InputFile <- "..\\.\\Data\\Manuscript_3D_RealData_Pattern\\"


InputData_Loc <- read.csv(paste0(InputFile, "3dst_RA1_raw_loc.csv"))
InputData_Exp <- read.csv(paste0(InputFile, "3dst_RA1_raw_exp.csv"))
InputData <- as.data.frame(t(InputData_Exp[,-1]))
InputData <- cbind(InputData_Loc, InputData)
colnames(InputData) <- c("Cell", "X", "Y", "Z", InputData_Exp[,1])


# Calculating SSH ---------------------------------------------------------

SSH_Calculation <- function(Interested_Gene){
  SSH_Fml <- as.formula(paste0(Interested_Gene, "~ X + Y"))
  SSH <- GD::gdm(SSH_Fml,
                 continuous_variable = c("X", "Y"),
                 data = InputData,
                 discmethod = c("equal","natural","quantile"), discitv = c(1:10))
  return(SSH$Factor.detector$Factor$sig)
}

SSH_Calculation("MAN1A2")
SSH_Calculation("SEMA4D")
