#############################################################
# Generate Figure 3A, 3B, Supplementary Figure 4,6 
#############################################################

library(ggplot2)
library(scales)
library(ggpubr)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)
library(ComplexHeatmap)
library(ComplexUpset)
library(ggvenn)
library(viridis)
#library(tidyverse)
#library(CellChat)

# kidney ------------------------------------------------------------------
wb <- createWorkbook()
InputData <- read.csv("..\\Data\\Manuscript_2D_RealData\\kidney085_XY01_20-0038_P_values_Fast.csv")
ego <- enrichGO(gene = InputData$X[which(InputData$p_values<0.05)],
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)
GeneEnrich <- dotplot(ego, showCategory=20, color="qvalue", label_format = 50)
png(file="..\\Outputs\\Manuscript_2D_RealData_Kidney_Enrich.png",width = 10,height = 10,units = "in",
    res = 600)
print(GeneEnrich)
dev.off()

Gene_IDs <- mapIds(org.Hs.eg.db,
                   keys=InputData$X[which(InputData$p_values<0.05)], 
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
egoKEGG <- enrichPathway(Gene_IDs)
ReacEnrich <- dotplot(egoKEGG, showCategory=20, color="qvalue", label_format = 50)
png(file="..\\Outputs\\Manuscript_2D_RealData_Kidney_Pathway.png",width = 10,height = 10,units = "in",
    res = 600)
print(ReacEnrich)
dev.off()

addWorksheet(wb, "kidney - GO enrichment")
writeData(wb, "kidney - GO enrichment",
          as.data.frame(ego))
addWorksheet(wb, "kidney - pathway")
writeData(wb, "kidney - pathway",
          as.data.frame(egoKEGG))

saveWorkbook(wb, paste0("..\\Data\\kidney_GO_JL_",
                        gsub("-", "", Sys.Date()),
                        ".xlsx"))

rm(list=ls())

# MOB ---------------------------------------------------------------------
Col_plate <- c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF",  
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF")
# InputData_MOB <- read.csv("..\\Data\\Manuscript_2D_RealData\\MOB_Pvalues.csv")
InputData_MOB <- read.csv("/Users/wangjuex/Library/CloudStorage/OneDrive-IndianaUniversity/BSP/Figshare/Archive/Manuscript_2D_RealData/MOB_Pvalues.csv")
# Methods <- colnames(InputData_MOB)[2:6]
Methods<-c("BSP","nnSVG","SPARK","SpatialDE") # SPARKX is 0, so we ignore it
for(Method in Methods){
  InputData_MOB[,Method] <- as.numeric(InputData_MOB[,Method]<0.05)
}

barintersection <- function(InputData){
  OutputData <- cbind(names(colSums(InputData)), colSums(InputData))
  for(Combn_num in 2:ncol(InputData)){
    for(Combn_index in 1:choose(ncol(InputData), Combn_num)){
      Intersection_Benchmarks <- sum(apply(InputData[,combn(1:ncol(InputData), Combn_num)[,Combn_index]], 1, min)==1)
      OutputData <- rbind(OutputData, 
                          c(paste0(colnames(InputData)[combn(1:ncol(InputData), Combn_num)[,Combn_index]],collapse = " "), Intersection_Benchmarks))
    }
  }
  OutputData <- as.data.frame(OutputData)
  rownames(OutputData) <- NULL
  colnames(OutputData) <- c("Methods", "Benchmarks")
  OutputData$Methods <- stringr::str_wrap(OutputData$Methods, width = 5)
  OutputData$Methods <- factor(OutputData$Methods, levels = unique(OutputData$Methods))
  OutputData$Benchmarks <- as.numeric(OutputData$Benchmarks)
  return(OutputData)
}

## Figure 3A
InputData_MOB_Inter <- barintersection(InputData_MOB[1:10,2:5])
Upset_fig1 <- ggplot(data=InputData_MOB_Inter, aes(x=Methods, y=Benchmarks, fill = Methods)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab(NULL)+ ylab("Number of Benchmarks") +
  scale_y_continuous(breaks= pretty_breaks()) + 
  scale_fill_viridis_d()
  # scale_fill_manual(values=c(Col_plate[1:4], terrain.colors(12)[1:11]))
# png(file="..\\Outputs\\Manuscript_2D_RealData_MOB.png",width = 10,height = 6,units = "in",
#     res = 600)
png(file="~/Downloads/Manuscript_2D_RealData_MOB.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig1)
dev.off()

## Supplementary Figure 4
InputData_MOB_venn <- lapply(Methods, function(Method){
  InputData_MOB$Gene[which(InputData_MOB[1:10,Method]==1)]
})
names(InputData_MOB_venn) <- Methods
Upset_fig1_2 <- ggvenn(
  InputData_MOB_venn, 
  fill_color = viridis(4),
  # fill_color = Col_plate[1:4],
  stroke_size = 0.5, set_name_size = 4
)
# png(file="..\\Outputs\\Manuscript_2D_RealData_MOB_venn.png",width = 10,height = 6,units = "in",
#     res = 600)
png(file="~/Downloads/Manuscript_2D_RealData_MOB_venn.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig1_2)
dev.off()

rm(list=ls())






# HBC ---------------------------------------------------------------------
Col_plate <- c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF",  
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF")
# InputData_HBC <- read.csv("..\\Data\\Manuscript_2D_RealData\\HBC_Pvalues.csv")
InputData_HBC <- read.csv("/Users/wangjuex/Library/CloudStorage/OneDrive-IndianaUniversity/BSP/Figshare/Archive/Manuscript_2D_RealData/HBC_Pvalues.csv")
Methods <- colnames(InputData_HBC)[2:5]
for(Method in Methods){
  InputData_HBC[,Method] <- as.numeric(InputData_HBC[,Method]<0.05)
}

barintersection <- function(InputData){
  OutputData <- cbind(names(colSums(InputData)), colSums(InputData))
  for(Combn_num in 2:ncol(InputData)){
    for(Combn_index in 1:choose(ncol(InputData), Combn_num)){
      Intersection_Benchmarks <- sum(apply(InputData[,combn(1:ncol(InputData), Combn_num)[,Combn_index]], 1, min)==1)
      OutputData <- rbind(OutputData, 
                          c(paste0(colnames(InputData)[combn(1:ncol(InputData), Combn_num)[,Combn_index]],collapse = " "), Intersection_Benchmarks))
    }
  }
  OutputData <- as.data.frame(OutputData)
  rownames(OutputData) <- NULL
  colnames(OutputData) <- c("Methods", "Benchmarks")
  OutputData$Methods <- stringr::str_wrap(OutputData$Methods, width = 5)
  OutputData$Methods <- factor(OutputData$Methods, levels = unique(OutputData$Methods))
  OutputData$Benchmarks <- as.numeric(OutputData$Benchmarks)
  return(OutputData)
}

## Figure 3B
InputData_HBC_Inter <- barintersection(InputData_HBC[1:14,2:5])
Upset_fig2 <- ggplot(data=InputData_HBC_Inter, aes(x=Methods, y=Benchmarks, fill = Methods)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab(NULL) + ylab("Number of Benchmarks") +
  scale_y_continuous(breaks= pretty_breaks()) + 
  scale_fill_viridis_d()
  # scale_fill_manual(values=c(Col_plate[1:4], terrain.colors(12)[1:11]))
# png(file="..\\Outputs\\Manuscript_2D_RealData_HBC.png",width = 10,height = 6,units = "in",
#     res = 600)
png(file="~/Downloads/Manuscript_2D_RealData_HBC.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig2)
dev.off()

## Supplementary Figure 6
InputData_HBC_venn <- lapply(Methods, function(Method){
  InputData_HBC$Gene[which(InputData_HBC[1:14,Method]==1)]
})
names(InputData_HBC_venn) <- Methods
Upset_fig2_2 <- ggvenn(
  InputData_HBC_venn, 
  # fill_color = Col_plate[1:4],
  fill_color = viridis(4),
  stroke_size = 0.5, set_name_size = 4
)
# png(file="..\\Outputs\\Manuscript_2D_RealData_HBC_venn.png",width = 10,height = 6,units = "in",
#     res = 600)
png(file="~/Downloads/Manuscript_2D_RealData_HBC_venn.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig2_2)
dev.off()


rm(list=ls())


# RA ---------------------------------------------------------
InputData_RA <- read.csv("..\\Data\\Manuscript_2D_RealData\\RA_SVGs.csv")

VennPlot <- ggvenn(list(`SVGs from 2D` = InputData_RA$Gene[which(InputData_RA$Source%in%c("Intersection", "2D"))],
                        `SVGs from 3D` = InputData_RA$Gene[which(InputData_RA$Source%in%c("Intersection", "3D"))]), 
                   c("SVGs from 2D", "SVGs from 3D"),
                   text_size = 5,
                   fill_color = c("#F8766D", "#00BA38"),
                   fill_alpha = 0.3)
png(file="..\\Outputs\\Manuscript_3D_RealData_RA_Venn.png",width = 10,height = 10,units = "in",
    res = 600)
print(VennPlot)
dev.off()

ego_inter <- enrichGO(gene = InputData_RA$Gene[which(InputData_RA$Source=="Intersection")],
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05)
GeneEnrich_Inter <- dotplot(ego_inter, showCategory=20, color="qvalue", label_format = 50)
png(file="..\\Outputs\\Manuscript_3D_RealData_RA_Inter.png",width = 10,height = 10,units = "in",
    res = 600)
print(GeneEnrich_Inter)
dev.off()

ego_3d <- enrichGO(gene = InputData_RA$Gene[which(InputData_RA$Source=="3D")],
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
GeneEnrich_3d <- dotplot(ego_3d, showCategory=20, color="qvalue", label_format = 50)
png(file="..\\Outputs\\Manuscript_3D_RealData_RA_3d.png",width = 10,height = 10,units = "in",
    res = 600)
print(GeneEnrich_3d)
dev.off()

ego_2d <- enrichGO(gene = InputData_RA$Gene[which(InputData_RA$Source=="2D")],
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

ego_3d_all <- enrichGO(gene = InputData_RA$Gene[which(InputData_RA$Source%in%c("3D", "Intersection"))],
                       OrgDb = org.Hs.eg.db,
                       keyType = 'SYMBOL',
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05)
GeneEnrich_3d_all <- dotplot(ego_3d_all, showCategory=20, color="qvalue", label_format = 50)
png(file="..\\Outputs\\Manuscript_3D_RealData_RA_3d_all.png",width = 10,height = 10,units = "in",
    res = 600)
print(GeneEnrich_3d_all)
dev.off()


wb <- createWorkbook()
addWorksheet(wb, "RA - 3D unique")
writeData(wb, "RA - 3D unique",
          as.data.frame(ego_3d))
addWorksheet(wb, "RA - 3D all")
writeData(wb, "RA - 3D all",
          as.data.frame(ego_3d_all))
addWorksheet(wb, "RA - intersection")
writeData(wb, "RA - intersection",
          as.data.frame(ego_inter))

saveWorkbook(wb, paste0("..\\Data\\Human_Rheumatoid_Arthritis_GO_JL_",
                        gsub("-", "", Sys.Date()),
                        ".xlsx"))

rm(list=ls())


# Large scale -------------------------------------------------------------
InputData_large <- read.table("..\\Data\\Manuscript_2D_RealData\\slideseq2_genes.txt",
                              header = TRUE)

ego_large <- enrichGO(gene = InputData_large$Gene,
                      OrgDb = org.Mm.eg.db,
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05)
GeneEnrich_large <- dotplot(ego_large, showCategory=20, color="qvalue", label_format = 50)
png(file="..\\Outputs\\Manuscript_2D_RealData_large.png",width = 10,height = 10,units = "in",
    res = 600)
print(GeneEnrich_large)
dev.off()


wb <- createWorkbook()
addWorksheet(wb, "Mouse brain")
writeData(wb, "Mouse brain",
          as.data.frame(ego_large))

saveWorkbook(wb, paste0("..\\Data\\Mouse_brain_GO_JL_",
                        gsub("-", "", Sys.Date()),
                        ".xlsx"))
