library(openxlsx)
library(data.table)
library(vegan)
library(ggplot2)
library(ggpubr)
library(colorhcplot)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
library (DESeq2)
mygenere <- read.xlsx ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/conteggio_generi1.xlsx")
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")

#T1

conditionlist <- c("Ster_soil", "Ster_root", "Soil")
mygenere$nome_ID <- gsub("_[^_]*$", "", mygenere$nome_ID)
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
biggenere <- mygenere
bigmetadata <- mymetadata
for (cond in 1:length(conditionlist)) 
{
  mygenere <- biggenere
  mymetadata <- bigmetadata
  condition <- conditionlist[cond]
  mymetadata <- mymetadata[,c("Sample_name", condition)]
  names(mymetadata)[2]<-"condition"
  rownames(mymetadata) <- mymetadata$Sample_name
  mymetadata$Sample_name <- NULL
  if(sum(rownames(mymetadata)==names(mygenere)) != ncol(mygenere)) stop("Names in countdata do not match names in metadata")
  #genero matrice per deseq che deve avere generi, metadata e la condizione da confrontare
  dds <- DESeqDataSetFromMatrix(countData = mygenere,
                                colData = mymetadata,
                                design = ~ condition)
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  newcounts <- as.data.frame(counts(dds, normalized = TRUE))
  write.xlsx(newcounts, "C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/normcount_T1.xlsx", row.names = TRUE)
  soglia <- order(rowMeans(newcounts), decreasing = TRUE)[1:25]
  newcounts$media <- rowMeans(newcounts)
  pcond <- data.frame(condition=mymetadata$condition,row.names=rownames(mymetadata))
  if(condition=="Soil"){
    anncolors <- list(condition =c("Manure" = "black", "Peat" = "red", "Sand"="gold"))
  }else {
    anncolors <- list(condition =c(Yes = "green",No = "red"))
  }
  pheatmap(assay(vsd)[soglia,], annotation_col = pcond, annotation_colors = anncolors, fontsize=4, cellwidth=6, cellheight=4, filename = paste0("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/immagini/deseq/", "heatmap_T1_", condition, ".pdf"))
}

