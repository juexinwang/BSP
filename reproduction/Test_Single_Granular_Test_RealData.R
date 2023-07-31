#############################################################
# Test for single granlarity
# Note: BSP uses paired granularities
#############################################################

library(EnvStats)
library(fitdistrplus)
library(ggplot2)
library(ggpubr)
InputPath <- "../Data/Real_Data/Test Score/"
OutputPath <- "../Data/Real_Data/GOF Test/" 

InputFiles <- list.files(InputPath)

for(InputFile in InputFiles){
  InputData <- read.csv(paste0(InputPath, InputFile))
  binwidth <- quantile(InputData$X3.0, 0.95)/50
  N <- nrow(InputData)
  Beta_par <- fitdist(InputData$X3.0, "beta")
  Beta_hist <- ggplot(InputData, aes(x=X3.0)) +
    geom_histogram(fill="#69b3a2", color="#69b3a2", alpha=0.6, position = 'identity', binwidth = binwidth) +
    labs(fill="") +
    xlim(c(0, quantile(InputData$X3.0, 0.95) ))+
    xlab("Test Score") +
    ylab("Frequency") +
    stat_function(fun = function(x) dbeta(x, Beta_par$estimate[1], Beta_par$estimate[2]) * N * binwidth,
                  color = "darkred")+
    ggtitle("fitted beta distribution")
  
  Beta_QQ <- ggplot(InputData, aes(sample = X3.0))+
    stat_qq(distribution = qbeta, dparams = c(Beta_par$estimate[1], Beta_par$estimate[2]))+ 
    geom_abline(intercept = 0, slope = 1, color="red", 
                linetype="dashed", size=1.5)+
    ylab("Observed") +
    xlim(c(0, quantile(InputData$X3.0, 0.95))) +
    ylim(c(0, quantile(InputData$X3.0, 0.95))) +
    ggtitle("Q-Q plot (beta)")
  png(file=paste0("../Outputs/", gsub(".csv", "", InputFile), "_Test_Score_Dist_Beta.png"),
      width=6, height=3, units = "in", res = 600)
  ggarrange(Beta_hist, Beta_QQ, nrow = 1, ncol = 2)
  dev.off()
  
  binwidth <-   binwidth <- quantile(InputData$X3.0, 0.95)/50
  N <- nrow(InputData)
  Lnorm_par <- fitdist(InputData$X3.0, "lnorm")
  LNorm_hist <- ggplot(InputData, aes(x=X3.0)) +
    geom_histogram(fill="#69b3a2", color="#69b3a2", alpha=0.6, position = 'identity', binwidth = binwidth) +
    labs(fill="") +
    xlim(c(0, quantile(InputData$X3.0, 0.95) ))+
    xlab("Test Score") +
    ylab("Frequency") +
    stat_function(fun = function(x) dlnorm(x, Lnorm_par$estimate[1], Lnorm_par$estimate[2]) * N * binwidth,
                  color = "darkred")+
    ggtitle("fitted lognormal distribution")
  
  LNorm_QQ <- ggplot(InputData, aes(sample = X3.0))+
    stat_qq(distribution = qlnorm, dparams = c(Lnorm_par$estimate[1], Lnorm_par$estimate[2]))+ 
    geom_abline(intercept = 0, slope = 1, color="red", 
                linetype="dashed", size=1.5) +
    xlim(c(0, quantile(InputData$X3.0, 0.95))) +
    ylim(c(0, quantile(InputData$X3.0, 0.95))) +
    ylab("Observed") +
    ggtitle("Q-Q plot (lognormal)")
  png(file=paste0("../Outputs/", gsub(".csv", "", InputFile), "_Test_Score_Dist_Lognorm.png"),
      width=6, height=3, units = "in", res = 600)
  ggarrange(LNorm_hist, LNorm_QQ, nrow = 1, ncol = 2)
  dev.off()
}



for(InputFile in InputFiles){
  InputData <- read.csv(paste0(InputPath, InputFile))
  P_values <- sapply(1:9, function(i){
    P_value = sapply(2:ncol(InputData), function(j){
      return(round(as.numeric(gofTest(InputData[1:100 + (i - 1)*100, j], 
                                      distribution = "beta", 
                                      test = "cvm")$p.value), 2))
    })
  })
  P_values <- as.data.frame(t(P_values))
  colnames(P_values) <- colnames(InputData)[2:ncol(InputData)]
  write.csv(P_values, row.names = FALSE,
            paste0(OutputPath, InputFile))
}


OutputPath <- "../Data/Real_Data/GOF Test2/" 

InputFiles <- list.files(InputPath)
for(InputFile in InputFiles){
  InputData <- read.csv(paste0(InputPath, InputFile))
  P_values <- sapply(1:9, function(i){
    P_value = sapply(2:ncol(InputData), function(j){
      return(round(as.numeric(gofTest(InputData[1:100 + (i - 1)*100, j], 
                                      distribution = "lnorm", 
                                      test = "cvm")$p.value), 2))
    })
  })
  P_values <- as.data.frame(t(P_values))
  colnames(P_values) <- colnames(InputData)[2:ncol(InputData)]
  write.csv(P_values, row.names = FALSE,
            paste0(OutputPath, InputFile))
}
