library(dplyr)
InputExp <- openxlsx::read.xlsx("../Data/mmc6.xlsx", colNames = FALSE)[,-1] %>% unlist() %>% as.numeric()
OutputDir <- "../Data/Simulation_3D_Test/"


sim_discrete <- function(InputExp, OutputDir,
                         NGenes = 1000, 
                         NX = 15, NY = 15, numSlices = 10, 
                         Rad = 2, NCenter = 100,
                         Quant = 0.80, Noise_lv = 1,
                         SEED = 1){
  set.seed(SEED)
  SimData <- data.frame()
  N <- NX * NY
  for(i in 1:numSlices){
    Coords <- spatstat.random::rpoispp(N, win = spatstat.geom::owin(c(0, 1), c(0, 1)))
    SimData <- rbind(SimData, data.frame(x = Coords$x * (sqrt(N) - 1), y = Coords$y * (sqrt(N) - 1), z = i * 10))
  }
  N_Obs <- nrow(SimData)
  
  set.seed(SEED)
  Sim_SVGs <- sapply(1:NGenes, function(Gene){
    return(sample(InputExp, N_Obs, replace = TRUE))
  })
  

  spikedcells <- function(coordDF, Rad, NCenter){
    Phi <- runif(NCenter, pi / 4, pi * 3 / 4)
    Theta <- 2 * pi * runif(NCenter)
    
    Center_pts <- cbind((NX + 1)/2 + cumsum(2 * cos(Phi) * sin(Theta)),
                        (NY + 1)/2 + cumsum(2 * cos(Phi) * cos(Theta)),
                        cumsum(2 * sin(Phi)) + Rad)
    #plot(cumsum(2*sin(Theta)), cumsum(2*cos(Theta)),type = "b", xlim = c(-10,10), ylim =c(-10,10) )
    
    spikedCoords <- sapply(1:length(Phi), function(Center_pt){
      # get distance from every cell to the center of the coordinate set; assuming hot spot is centered here
      distFromCenter <- sp::spDists(
        x = coordDF %>% dplyr::select(x, y, z) %>% as.matrix(),
        y = matrix(Center_pts[Center_pt, ], ncol = 3),
        longlat = FALSE
      ) %>% as.vector()
      
      # identify coordinates in hot spot
      spikedCoord <- coordDF[distFromCenter <= Rad,] %>% row.names() %>% as.numeric()
    })
    spikedCoords <- unique(unlist(spikedCoords))
    
    # return coordinates
    return(spikedCoords)
  }

  SpikedCellIndex <- spikedcells(SimData, Rad, NCenter)
  SpikedCells <- length(SpikedCellIndex)
  SpikedValues_upper <- InputExp[InputExp>=quantile(InputExp, Quant)]
  SpikedCellValues_upper <- sapply(1:NGenes, function(Gene){
    return(sample(SpikedValues_upper, SpikedCells, replace = TRUE))
  })
  SpikedCellValues <- SpikedCellValues_upper
  SpikedCellValues <- as.matrix(SpikedCellValues)
  Sim_SVGs[SpikedCellIndex, ] <- SpikedCellValues
  
  colnames(Sim_SVGs) <- paste0("SVG_", 1:NGenes)
  SimData$z <- SimData$z / 10
  SimData <- cbind(SimData, Sim_SVGs)
  
  perm_func <- function(Input_Data, Rep = 9){
    NumGenes <- ncol(Input_Data) - 3
    NumCells <- nrow(Input_Data)
    set.seed(SEED)
    Output_Data <- sapply(1:NumGenes, function(Gene_Index){
      return(as.data.frame(sapply(1:Rep, 
                                  function(j) 
                                    return(Input_Data[sample(1:NumCells, 
                                                             replace = FALSE), 
                                                      Gene_Index + 3])
      )
      )
      )
    })
    Output_Data <- do.call("cbind", Output_Data)
    return(Output_Data)
  }
  Sim_NullGenes <- perm_func(SimData)
  colnames(Sim_NullGenes) <- paste0("NULL_", 1:(NGenes * 9))
  
  SimData <- cbind(SimData, Sim_NullGenes)
  
  SD_mean <- mean(apply(SimData[,4:ncol(SimData)], 2, sd))
  Noise_Matrix<- matrix(round(Noise_lv*rnorm((nrow(SimData) * (ncol(SimData)-3)), mean = 0, sd = SD_mean)),
                        nrow = nrow(SimData), ncol = (ncol(SimData)-3))
  SimData_AllExp <- SimData[,4:ncol(SimData)]
  SimData_AllExp <- SimData_AllExp + Noise_Matrix
  SimData_AllExp[SimData_AllExp<0] <- 0
  SimData[,4:ncol(SimData)] <- SimData_AllExp
  
  
  write.csv(SimData, 
            paste0(OutputDir,
                   "3DStack_",
                   N,
                   "by",
                   numSlices,
                   "_width",
                   Rad,
                   "_qt",
                   Quant*100,
                   "_Noise",
                   Noise_lv,
                   "_pw",
                   SEED,
                   ".csv"), row.names = FALSE, quote = FALSE)
  
}






