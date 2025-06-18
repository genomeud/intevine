library(DESeq2)
library(data.table)
library(openxlsx)
library(gtools)
genes <- fread("/projects/marroni/intevine/analysis/transcriptomic/rna_table.txt",data.table=F)
rownames(genes)<- genes$Gene
genes$Gene<-NULL
genes[is.na(genes)]=0
genes<- genes[5:nrow(genes),]
 colnames(genes)<-gsub("Radice","",colnames(genes))
metadata<-read.xlsx("/projects/marroni/intevine/docs/metadata.xlsx")
metadata$Sample_label <- NULL
names(metadata)[7]<-"condition"
rownames(metadata)<- metadata$Sample_name
ddsHTSeqq <- DESeqDataSetFromMatrix(countData  = genes,
                                        colData = metadata,
                                        design = ~ condition)
ddsHTSeqq <- DESeq(ddsHTSeqq)		
p <- counts(ddsHTSeqq, normalized=TRUE)
p <- as.data.frame(p)
p <- tibble::rownames_to_column(p, "Gene")			
write.table(p, "norm_counts.txt", quote=F, sep="\t")