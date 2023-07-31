library(dplyr)
InputExp <- openxlsx::read.xlsx("../Data/mmc6.xlsx", colNames = FALSE)[,-1] %>% unlist() %>% as.numeric()
OutputDir <- "../Data/Simulation_3D_Test/"


sim_discrete <- function(InputExp, OutputDir,
                         NGenes = 1000, 
                         NX = 30, NY = 30, numSlices = 10, 
                         Rad = 2, Breaks = 8,
                         Quant = 0.80, Noise_lv = 2,
                         SEED = 1){
  set.seed(SEED)
  SimData <- data.frame()
  N <- NX * NY
  for(i in 1:numSlices){
    Coords <- spatstat.random::rpoispp(N, win = spatstat.geom::owin(c(0, 1), c(0, 1)))
    SimData <- rbind(SimData, data.frame(x = Coords$x * (sqrt(N) - 1), y = Coords$y * (sqrt(N) - 1), z = i))
  }
  N_Obs <- nrow(SimData)
  
  set.seed(SEED)
  Sim_SVGs <- sapply(1:NGenes, function(Gene){
    return(sample(InputExp, N_Obs, replace = TRUE))
  })
  
  # local patterns / discrete patterns
  spikedcells <- function(coordDF, Rad, Breaks){
    Center_pts <- as.data.frame(expand.grid(seq(from = 3, to = NX, by = Breaks), 
                                           seq(from = 3, to = NY, by = Breaks)))
    Center_pts$Z <- (numSlices + 1) / 2
    Center_pts <- Center_pts + runif(length(Center_pts), -2, 2)
    
    spikedCoords <- sapply(1:nrow(Center_pts), function(Center_pt){
      # get distance from every cell to the center of the coordinate set; assuming hot spot is centered here
      distFromCenter <- sp::spDists(
        x = coordDF %>% dplyr::select(x, y, z) %>% as.matrix(),
        y = matrix(as.numeric(Center_pts[Center_pt, ]), ncol = 3),
        longlat = FALSE
      ) %>% as.vector()
      
      # identify coordinates in hot spot
      spikedCoord <- coordDF[distFromCenter <= (Rad),] %>% row.names() %>% as.numeric()
      return(spikedCoord)
    })
    spikedCoords <- unique(unlist(spikedCoords))
    return(spikedCoords)
  }
  
  # global patterns / continuous patterns
  spikedcells2 <- function(coordDF, Rad, Breaks){
    NCenter <- length(seq(from = 3, to = NX, by = Breaks)) * length(seq(from = 3, to = NY, by = Breaks))
    
    Phi <- 2 * pi * runif(NCenter)
    Theta <- 2 * pi * runif(NCenter)
    
    Center_pts <- cbind((NX + 1)/2 + cumsum(2 * cos(Phi) * sin(Theta)),
                        (NY + 1)/2 + cumsum(2 * cos(Phi) * cos(Theta)),
                        cumsum(2 * sin(Phi)) + numSlices/2)
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
  
  
  Sim_SVGs1 <- Sim_SVGs[,1:floor(NGenes/2)]
  Sim_SVGs2 <- Sim_SVGs[,(1:ceiling(NGenes/2) + floor(NGenes/2))]
  
  # local pattern
  SpikedCellIndex <- spikedcells(SimData, Rad, Breaks)
  SpikedCells <- length(SpikedCellIndex)
  SpikedValues_upper <- InputExp[InputExp>=quantile(InputExp, Quant)]
  SpikedCellValues_upper <- sapply(1:floor(NGenes/2), function(Gene){
    return(sample(SpikedValues_upper, SpikedCells, replace = TRUE))
  })
  SpikedCellValues <- SpikedCellValues_upper
  SpikedCellValues <- as.matrix(SpikedCellValues)
  Sim_SVGs1[SpikedCellIndex, ] <- SpikedCellValues
  
  # global pattern
  SpikedCellIndex <- spikedcells2(SimData, Rad, Breaks)
  SpikedCells <- length(SpikedCellIndex)
  SpikedValues_upper <- InputExp[InputExp>=quantile(InputExp, Quant)]
  SpikedCellValues_upper <- sapply(1:ceiling(NGenes/2), function(Gene){
    return(sample(SpikedValues_upper, SpikedCells, replace = TRUE))
  })
  SpikedCellValues <- SpikedCellValues_upper
  SpikedCellValues <- as.matrix(SpikedCellValues)
  Sim_SVGs2[SpikedCellIndex, ] <- SpikedCellValues
  
  Sim_SVGs <- cbind(Sim_SVGs1, Sim_SVGs2)
  colnames(Sim_SVGs) <- paste0("SVG_", 1:NGenes)
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
  Noise_Matrix<- matrix(round(1*rnorm((nrow(SimData) * (ncol(SimData)-3)), mean = 0, sd = SD_mean)),
                        nrow = nrow(SimData), ncol = (ncol(SimData)-3))
  SimData_AllExp <- SimData[,4:ncol(SimData)]
  SimData_AllExp <- SimData_AllExp + Noise_Matrix
  SimData_AllExp[SimData_AllExp<0] <- 0
  SimData[,4:ncol(SimData)] <- SimData_AllExp
  
  
  
  write.csv(SimData, 
            paste0(OutputDir,
                   "Discrete_",
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






