# Plot Figure 3, Supplementary Figure 4,6,7

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
InputData_MOB <- read.csv("..\\Data\\Manuscript_2D_RealData\\MOB_Pvalues.csv")
Methods <- colnames(InputData_MOB)[2:5]
for(Method in Methods){
  InputData_MOB[,Method] <- as.numeric(InputData_MOB[,Method]<0.05)
}
InputData_MOB$Benchmarks <- 0
InputData_MOB$Benchmarks[1:10] <- 1
InputData_MOB$Benchmarks <- as.factor(InputData_MOB$Benchmarks)
InputData_MOB_Benchmarks <- reshape::melt(InputData_MOB[1:10,], id = "Gene")
InputData_MOB_Benchmarks$value <- as.factor(InputData_MOB_Benchmarks$value)
Upset_fig1 <- upset(InputData_MOB[1:10,], 
      sort_intersections = FALSE,intersect = c("BSP", "SPARK", "SPARKX", "SpatialDE"),
      intersections = list(c('BSP'),
                           c('SPARK'),
                           c('BSP', 'SPARK'),
                           c('BSP', 'SPARK', 'SpatialDE')),
      queries=list(upset_query( intersect=c('BSP'), color=Col_plate[1],fill=Col_plate[1]),
        upset_query( intersect=c('SPARK'), color=Col_plate[2],fill=Col_plate[2])))
png(file="..\\Outputs\\Manuscript_2D_RealData_MOB.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig1)
dev.off()



Upset_fig1_2 <- upset(InputData_MOB, 
                    sort_intersections = FALSE,
                    intersect = c("BSP", "SPARK", "SPARKX", "SpatialDE"),
                    intersections = c(list("BSP", "SPARK", "SPARKX", "SpatialDE"),
                                         lapply(1:ncol(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),2)), function(i) return(unlist(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),2)[,i]))),
                                         lapply(1:ncol(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),3)), function(i) return(unlist(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),3)[,i]))),
                                         list(c("BSP", "SPARK", "SPARKX", "SpatialDE"))),
                    queries=list(upset_query( intersect=c('BSP'), color=Col_plate[1],fill=Col_plate[1]),
                                 upset_query( intersect=c('SPARK'), color=Col_plate[2],fill=Col_plate[2]),
                                 upset_query( intersect=c('SPARKX'), color=Col_plate[3],fill=Col_plate[3]),
                                 upset_query( intersect=c('SpatialDE'), color=Col_plate[4],fill=Col_plate[4])))
png(file="..\\Outputs\\Manuscript_2D_RealData_MOB_all.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig1_2)
dev.off()
rm(list=ls())






# HBC ---------------------------------------------------------------------
Col_plate <- c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF",  
               "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#DC0000FF")
InputData_HBC <- read.csv("..\\Data\\Manuscript_2D_RealData\\HBC_Pvalues.csv")
Methods <- colnames(InputData_HBC)[2:5]
for(Method in Methods){
  InputData_HBC[,Method] <- as.numeric(InputData_HBC[,Method]<0.05)
}
InputData_HBC$Benchmarks <- 0
InputData_HBC$Benchmarks[1:14] <- 1
InputData_HBC$Benchmarks <- as.factor(InputData_HBC$Benchmarks)
InputData_HBC_Benchmarks <- reshape::melt(InputData_HBC[1:14,], id = "Gene")
InputData_HBC_Benchmarks$value <- as.factor(InputData_HBC_Benchmarks$value)
Upset_fig2 <- upset(InputData_HBC[1:14,], 
                    sort_intersections = FALSE,intersect = c("BSP", "SPARK", "SPARKX", "SpatialDE"),
                    intersections = list(c('BSP'),
                                         c('BSP', 'SPARK', 'SPARKX', 'SpatialDE'),
                                         c('BSP', 'SPARK'),
                                         c('BSP',  'SPARKX'),
                                         c('SPARK', 'SPARKX', 'SpatialDE'),
                                         c('BSP', 'SPARK', 'SpatialDE'),
                                         c('BSP', 'SPARK', 'SPARKX')),
                    queries=list(
                      upset_query( intersect=c('BSP'), color=Col_plate[1],fill=Col_plate[1])))
png(file="..\\Outputs\\Manuscript_2D_RealData_HBC.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig2)
dev.off()


Upset_fig2_2 <- upset(InputData_HBC, 
                      sort_intersections = FALSE,
                      intersect = c("BSP", "SPARK", "SPARKX", "SpatialDE"),
                      intersections = c(list("BSP", "SPARK", "SPARKX", "SpatialDE"),
                                        lapply(1:ncol(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),2)), function(i) return(unlist(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),2)[,i]))),
                                        lapply(1:ncol(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),3)), function(i) return(unlist(combn(c("BSP", "SPARK", "SPARKX", "SpatialDE"),3)[,i]))),
                                        list(c("BSP", "SPARK", "SPARKX", "SpatialDE"))),
                      queries=list(upset_query( intersect=c('BSP'), color=Col_plate[1],fill=Col_plate[1]),
                                   upset_query( intersect=c('SPARK'), color=Col_plate[2],fill=Col_plate[2]),
                                   upset_query( intersect=c('SPARKX'), color=Col_plate[3],fill=Col_plate[3]),
                                   upset_query( intersect=c('SpatialDE'), color=Col_plate[4],fill=Col_plate[4])))
png(file="..\\Outputs\\Manuscript_2D_RealData_HBC_all.png",width = 10,height = 6,units = "in",
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



# RA 3d ---------------------------------------------------------
InputData_RA3d <- read.table("..\\Data\\Manuscript_2D_RealData\\RA_3d_pvalues.txt", header = TRUE)

GO_Results <- list()
for(Subject in 1:6){
  Sig_genes <- InputData_RA3d$Genes[which(InputData_RA3d[,paste0("RA", Subject)]<0.05)]
  ego_subject <- enrichGO(gene = Sig_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = 'SYMBOL',
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05)
  GeneEnrich_3dsubject <- dotplot(ego_subject, showCategory=30, color="qvalue", label_format = 50)
  GO_Results[[Subject]] <- ego_subject[,1]
  
  png(file=paste0("..\\Outputs\\Manuscript_3D_RealData_RA_", Subject ,".png"),
      width = 10,height = 10,units = "in",
      res = 600)
  print(GeneEnrich_3dsubject)
  dev.off()
}
names(GO_Results) <- paste0("RA", 1:6)

GO_Results_table <- lapply(GO_Results, function(RA_GO) return(data.frame(GO_term = unlist(RA_GO),
                                                                         Value = 1)))

GO_Results_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GO_term", all = TRUE),
                           GO_Results_table)
colnames(GO_Results_table)[2:7] <- paste0("RA", 1:6)
GO_Results_table[is.na(GO_Results_table)] <- 0
Upset_fig_all <- upset(GO_Results_table, 
      width_ratio = 0.1,   
      base_annotations=list(
        'Intersection size'=intersection_size(counts=FALSE)
      ),
      intersect = paste0("RA", 1:6))


png(file="..\\Outputs\\Manuscript_2D_RealData_RA_All.png",width = 10,height = 6,units = "in",
    res = 600)
print(Upset_fig_all)
dev.off()



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


InputData_Celltypes <- read.table("..\\Data\\Manuscript_2D_RealData\\slideseq2_celltypes.txt",
                              header = FALSE, sep = "\t")
colnames(InputData_Celltypes)[3] <- "Cell_Type"
InputData_Celltypes <- InputData_Celltypes[which(!is.na(InputData_Celltypes$Cell_Type)),]
CellType_Table <- InputData_Celltypes %>% 
  group_by(Cell_Type) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(per=`n`/sum(`n`)) %>% 
  arrange(desc(Cell_Type))
# aggregate cell types with less than 20 genes
CellType_Table <- CellType_Table[order(CellType_Table$n, decreasing = TRUE),]
CellType_Table <- rbind(CellType_Table[which(CellType_Table$n>=20),],
                        c("Others", as.numeric(colSums(CellType_Table[which(CellType_Table$n<20),2:3]))))
CellType_Table$per <- as.numeric(CellType_Table$per)
CellType_Table$label <- paste0(round(CellType_Table$per * 100), "%")


PieChart <- ggplot(CellType_Table, aes("", per, fill = Cell_Type)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5)) +
  guides(fill = guide_legend(reverse = TRUE, title = "Cell Types")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position="bottom")

png(file="..\\Outputs\\Manuscript_2D_RealData_CellType.png",width = 10,height = 6,units = "in",
    res = 600)
print(PieChart)
dev.off()
