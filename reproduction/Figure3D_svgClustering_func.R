################################################################
# Generate Figure 3D, Supplementary Figure 25, and
# Corresponding Supplementary Figures 26-33
################################################################

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
# InputData <- read.csv("/Users/wangjuex/projects/svgGranularity/kidney_0_genes.txt",header=FALSE)
InputData <- read.csv("/Users/wangjuex/projects/svgGranularity/kidney_1_genes.txt",header=FALSE)

ego <- enrichGO(gene = InputData$V1,
                 OrgDb = org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
GeneEnrich <- dotplot(ego, showCategory=20, color="qvalue", label_format = 50)
# png(file="/Users/wangjuex/projects/svgGranularity/kidney_0_Kidney_Enrich.png",width = 10,height = 10,units = "in",
#     res = 600)
png(file="/Users/wangjuex/projects/svgGranularity/kidney_1_Kidney_Enrich.png",width = 10,height = 10,units = "in",
    res = 600)
print(GeneEnrich)
dev.off()

Gene_IDs <- mapIds(org.Hs.eg.db,
                   keys=InputData$V1, 
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
egoKEGG <- enrichPathway(Gene_IDs)
ReacEnrich <- dotplot(egoKEGG, showCategory=20, color="qvalue", label_format = 50)
# png(file="/Users/wangjuex/projects/svgGranularity/kidney_0_Kidney_Pathway.png",width = 10,height = 10,units = "in",
#     res = 600)
png(file="/Users/wangjuex/projects/svgGranularity/kidney_1_Kidney_Pathway.png",width = 10,height = 10,units = "in",
    res = 600)
print(ReacEnrich)
dev.off()

addWorksheet(wb, "kidney - GO enrichment")
writeData(wb, "kidney - GO enrichment",
          as.data.frame(ego))
addWorksheet(wb, "kidney - pathway")
writeData(wb, "kidney - pathway",
          as.data.frame(egoKEGG))

# saveWorkbook(wb, paste0("/Users/wangjuex/projects/svgGranularity/kidney_0_GO_",
#                         gsub("-", "", Sys.Date()),
#                         ".xlsx"))

saveWorkbook(wb, paste0("/Users/wangjuex/projects/svgGranularity/kidney_1_GO_",
                        gsub("-", "", Sys.Date()),
                        ".xlsx"))

rm(list=ls())


# RA ------------------------------------------------------------------
wb <- createWorkbook()
InputData <- read.csv("/Users/wangjuex/projects/svgGranularity/RA_0_genes.txt",header=FALSE)
# InputData <- read.csv("/Users/wangjuex/projects/svgGranularity/RA_1_genes.txt",header=FALSE)
# InputData <- read.csv("/Users/wangjuex/projects/svgGranularity/RA_2_genes.txt",header=FALSE)
# InputData <- read.csv("/Users/wangjuex/projects/svgGranularity/RA_3_genes.txt",header=FALSE)

ego <- enrichGO(gene = InputData$V1,
                 OrgDb = org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
GeneEnrich <- dotplot(ego, showCategory=20, color="qvalue", label_format = 50)
png(file="/Users/wangjuex/projects/svgGranularity/RA_0_RA_Enrich.png",width = 10,height = 10,units = "in",
    res = 600)
# png(file="/Users/wangjuex/projects/svgGranularity/RA_1_RA_Enrich.png",width = 10,height = 10,units = "in",
#     res = 600)
# png(file="/Users/wangjuex/projects/svgGranularity/RA_2_RA_Enrich.png",width = 10,height = 10,units = "in",
#     res = 600)
# png(file="/Users/wangjuex/projects/svgGranularity/RA_3_RA_Enrich.png",width = 10,height = 10,units = "in",
#     res = 600)
print(GeneEnrich)
dev.off()

Gene_IDs <- mapIds(org.Hs.eg.db,
                   keys=InputData$V1, 
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
egoKEGG <- enrichPathway(Gene_IDs)
ReacEnrich <- dotplot(egoKEGG, showCategory=20, color="qvalue", label_format = 50)
png(file="/Users/wangjuex/projects/svgGranularity/RA_0_RA_Pathway.png",width = 10,height = 10,units = "in",
    res = 600)
# png(file="/Users/wangjuex/projects/svgGranularity/RA_1_RA_Pathway.png",width = 10,height = 10,units = "in",
#     res = 600)
# png(file="/Users/wangjuex/projects/svgGranularity/RA_2_RA_Pathway.png",width = 10,height = 10,units = "in",
#     res = 600)
# png(file="/Users/wangjuex/projects/svgGranularity/RA_3_RA_Pathway.png",width = 10,height = 10,units = "in",
#     res = 600)
print(ReacEnrich)
dev.off()

addWorksheet(wb, "RA - GO enrichment")
writeData(wb, "RA - GO enrichment",
          as.data.frame(ego))
addWorksheet(wb, "RA - pathway")
writeData(wb, "RA - pathway",
          as.data.frame(egoKEGG))

saveWorkbook(wb, paste0("/Users/wangjuex/projects/svgGranularity/RA_0_GO_",
                        gsub("-", "", Sys.Date()),
                        ".xlsx"))

# saveWorkbook(wb, paste0("/Users/wangjuex/projects/svgGranularity/RA_1_GO_",
#                         gsub("-", "", Sys.Date()),
#                         ".xlsx"))

# saveWorkbook(wb, paste0("/Users/wangjuex/projects/svgGranularity/RA_2_GO_",
#                         gsub("-", "", Sys.Date()),
#                         ".xlsx"))

# saveWorkbook(wb, paste0("/Users/wangjuex/projects/svgGranularity/RA_3_GO_",
#                         gsub("-", "", Sys.Date()),
#                         ".xlsx"))

rm(list=ls())