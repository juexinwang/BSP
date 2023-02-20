## Generate simulation data

InputData <- openxlsx::read.xlsx("mmc6.xlsx",
                                 colNames =  FALSE)
InputData_exp <- as.numeric(unlist(InputData[,-1]))
 
Num_genes <- c(2000)
Num_cells <- c(1000000)
for(Num_cell in Num_cells){
  for(Num_gene in Num_genes){
    SimData <- data.frame()
    Coords <- spatstat.random::rpoispp(Num_cell, win = spatstat.geom::owin(c(0, 1), c(0, 1)))
    N <- Coords$n
    SimData <- rbind(SimData, data.frame(x = Coords$x * (sqrt(N) - 1), y = Coords$y * (sqrt(N) - 1)))
    Sim_exp_mat <- sapply(1:Num_gene, function(i){
      Sim_exp <- sample(as.numeric(InputData_exp), N, replace＝TRUE)
      return(Sim_exp)
    })
    colnames(Sim_exp_mat) <- paste0("Gene_", 1:Num_gene)
    SimData <- cbind(SimData, Sim_exp_mat)
    write.csv(SimData,row.names = FALSE,
              paste0("Gene",
                     Num_gene,
                     "_Cell",
                     as.character(as.numeric(Num_cell)),
                     ".csv"))
  }
}