#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT0"), type="character", default="/home/massimo.guazzini/0.001_raw_counts_edited_16s_T0.xlsx",
              help="Input directory for T0 [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1"), type="character", default="/home/massimo.guazzini/0.001_raw_counts_edited_16s_T1.xlsx",
              help="Input directory for T1 [default= %default]", metavar="character"),
  make_option(c("-C", "--abundanceITST0"), type="character", default="/home/massimo.guazzini/0.001_raw_counts_edited_ITS_T0.xlsx",
              help="Input directory for ITS T0 [default= %default]", metavar="character"),
  make_option(c("-D", "--abundanceITST1"), type="character", default="/home/massimo.guazzini/0.001_raw_counts_edited_ITS_T1.xlsx",
              help="Input directory for ITS T1 [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/home/massimo.guazzini/map_grapevine_final.xlsx", 
              help="Metafile directory [default= %default]", metavar="character"),
  make_option(c("-O", "--outputT0"), type="character", default="/home/massimo.guazzini/deseq2/Species_aggregated/16sT0/", 
              help="Output directory for T0 [default= %default]", metavar="character"),
  make_option(c("-P", "--outputT1"), type="character", default="/home/massimo.guazzini/deseq2/Species_aggregated/16sT1/", 
              help="Output directory for T1 [default= %default]", metavar="character"),
  make_option(c("-Q", "--outputITST0"), type="character", default="/home/massimo.guazzini/deseq2/Species_aggregated/ITST0/", 
              help="Output directory for ITS T0 [default= %default]", metavar="character"),
  make_option(c("-R", "--outputITST1"), type="character", default="/home/massimo.guazzini/deseq2/Species_aggregated/ITST1/", 
              help="Output directory for ITS T1 [default= %default]", metavar="character"),
  make_option(c("-S", "--outtableT0"), type="character", default="/home/massimo.guazzini/deseq2/T0/06_tables_silva/DESeq2/Species_aggregated/", 
              help="Output table directory for T0 [default= %default]", metavar="character"),
  make_option(c("-T", "--outtableT1"), type="character", default="/home/massimo.guazzini/deseq2/T1/06_tables_silva/DESeq2/Species_aggregated/", 
              help="Output table directory for T1 [default= %default]", metavar="character"),
  make_option(c("-U", "--outtableITST0"), type="character", default="/home/massimo.guazzini/deseq2/T0_ITS/06_tables_unite/DESeq2/Species_aggregated/", 
              help="Output table directory for ITS T0 [default= %default]", metavar="character"),
  make_option(c("-V", "--outtableITST1"), type="character", default="/home/massimo.guazzini/deseq2/T1_ITS/06_tables_unite/DESeq2/Species_aggregated/", 
              help="Output table directory for ITS T1 [default= %default]", metavar="character"),
  make_option(c("-W", "--conditionT0"), type="character", default="Soil,Ster_soil", 
              help="Condition for T0 [default= %default]", metavar="character"),
  make_option(c("-X", "--condition"), type="character", default="Soil,Ster_soil,Ster_root", 
              help="Condition [default= %default]", metavar="character"),
  make_option(c("-Y", "--listrelevant"), type="character", default="/home/massimo.guazzini/Supplementary_file_8_list_of_agricultural_relevant_microorganisms.xlsx", 
              help="List of beneficial bacteria [default= %default]", metavar="character"),
  make_option(c("-Y", "--bacben"), type="character", default="Bacteria_beneficial", 
              help="List of beneficial bacteria [default= %default]", metavar="character"),
  make_option(c("-Z", "--bacdam"), type="character", default="Damaging_bacteria", 
              help="List of damaging bacteria [default= %default]", metavar="character"),
  make_option(c("-E", "--funben"), type="character", default="Fungal_beneficial", 
              help="List of beneficial fungi [default= %default]", metavar="character"),
  make_option(c("-F", "--fundam"), type="character", default="Damaging_fungi", 
              help="List of damaging fungi [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser) 

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundanceT0)) {
  stop("WARNING: No abundanceT0 specified with '-A' flag.")
} else {  
  cat("abundanceT0 is ", opt$abundanceT0, "\n")
  aT0 <- opt$abundanceT0  
}

if (is.null(opt$abundanceT1)) {
  stop("WARNING: No abundanceT1 specified with '-B' flag.")
} else {  
  cat("abundanceT1 is ", opt$abundanceT1, "\n")
  aT1 <- opt$abundanceT1  
}

if (is.null(opt$abundanceITST0)) {
  stop("WARNING: No abundanceITST0 specified with '-C' flag.")
} else {  
  cat("abundanceITST0 is ", opt$abundanceITST0, "\n")
  aITST0 <- opt$abundanceITST0  
}

if (is.null(opt$abundanceITST1)) {
  stop("WARNING: No abundanceITST1 specified with '-D' flag.")
} else {  
  cat("abundanceITST1 is ", opt$abundanceITST1, "\n")
  aITST1 <- opt$abundanceITST1  
}

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  
  cat("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
}

if (is.null(opt$outputT0)) {
  stop("WARNING: No outputT0 specified with '-O' flag.")
} else {  
  cat("outputT0 dir is ", opt$outputT0, "\n")
  oT0 <- opt$outputT0  
}

if (is.null(opt$outputT1)) {
  stop("WARNING: No outputT1 specified with '-P' flag.")
} else {  
  cat("outputT1 dir is ", opt$outputT1, "\n")
  oT1 <- opt$outputT1  
}

if (is.null(opt$outputITST0)) {
  stop("WARNING: No outputITST0 directory specified with '-Q' flag.")
} else {  
  cat("outputITST0 dir is ", opt$outputITST0, "\n")
  oITST0 <- opt$outputITST0 
}

if (is.null(opt$outputITST1)) {
  stop("WARNING: No outputITST1 directory specified with '-R' flag.")
} else {  
  cat("outputITST1 is ", opt$outputITST1, "\n")
  oITST1 <- opt$outputITST1
}

if (is.null(opt$outtableT0)) {
  stop("WARNING: No outtableT0 specified with '-S' flag.")
} else {  
  cat("outtableT0 dir is ", opt$outtableT0, "\n")
  otableT0 <- opt$outtableT0  
}

if (is.null(opt$outtableT1)) {
  stop("WARNING: No outtableT1 specified with '-T' flag.")
} else {  
  cat("outtableT1 dir is ", opt$outtableT1, "\n")
  otableT1 <- opt$outtableT1  
}

if (is.null(opt$outtableITST0)) {
  stop("WARNING: No outtableITST0 specified with '-U' flag.")
} else {  
  cat("outtableITST0 dir is ", opt$outtableITST0, "\n")
  otableITST0 <- opt$outtableITST0  
}

if (is.null(opt$outtableITST1)) {
  stop("WARNING: No outtableITST1 specified with '-V' flag.")
} else {  
  cat("outtableITST1 dir is ", opt$outtableITST1, "\n")
  otableITST1 <- opt$outtableITST1 
}

if (is.null(opt$conditionT0)) {
  stop("WARNING: No conditionT0 specified with '-W' flag.")
} else {  
  cat("conditionT0 is ", opt$conditionT0, "\n")
  conditionT0  <- opt$conditionT0 
}

if (is.null(opt$condition)) {
  stop("WARNING: No condition specified with '-X' flag.")
} else {  
  cat("condition is ", opt$condition, "\n")
  condition <- opt$condition 
}

if (is.null(opt$listrelevant)) {
  stop("WARNING: No listrelevant specified with '-X' flag.")
} else {  
  cat("listrelevant is ", opt$listrelevant, "\n")
  listrelevant <- opt$listrelevant 
}

if (is.null(opt$bacben)) {
  stop("WARNING: No bacben specified with '-Y' flag.")
} else {  
  cat("bacben is ", opt$bacben, "\n")
  bacben <- opt$bacben 
}

if (is.null(opt$bacdam)) {
  stop("WARNING: No bacdam specified with '-Z' flag.")
} else {  
  cat("bacdam is ", opt$bacdam, "\n")
  bacdam <- opt$bacdam 
}

if (is.null(opt$funben)) {
  stop("WARNING: No funben specified with '-E' flag.")
} else {  
  cat("funben is ", opt$funben, "\n")
  funben <- opt$funben 
}

if (is.null(opt$fundam)) {
  stop("WARNING: No fundam specified with '-F' flag.")
} else {  
  cat("fundam is ", opt$fundam, "\n")
  fundam <- opt$fundam 
}

runDEseq<-function() 
{
    library(DESeq2)
    library(data.table)
	library(ape)
	library(openxlsx)
	library("RColorBrewer")
	library("pheatmap")
	library("ggplot2")
	library("ggrepel")
  library("stringr")
  library("dplyr")
  library("gtools")
  library("vegan")
  inputlist <- c(aT0,aT1,aITST0,aITST1)
  damag <- c(bacdam,bacdam,fundam,fundam)
  benef <- c(bacben,bacben,funben,funben)
  outputlist <- c(oT0,oT1,oITST0,oITST1)
     tablelist <- c(otableT0, otableT1,otableITST0,otableITST1)
	  metadatalist <- c("T0", "Campionate_edited", "T0", "Campionate_edited")
  condlist <- c(conditionT0,condition,conditionT0,condition)
     for (myinput in 1:length(inputlist)) {	 
	 damagio <- data.frame(class = "formerging", match = "formerging", stringsAsFactors = FALSE)
	 print("RESET")
	 	 damagiogrape <- data.frame(class = "formerging", match = "formerging", stringsAsFactors = FALSE)
benefo <- data.frame(class = "formerging", match = "formerging", stringsAsFactors = FALSE)
	 damagiogrape <- data.frame(class = "formerging", match = "formerging", stringsAsFactors = FALSE)
	 print("RESET")
benefogrape <- data.frame(class = "formerging", match = "formerging", stringsAsFactors = FALSE)
	allcond<-unlist(strsplit(condlist[myinput],","))
 metadata<-read.xlsx(metafile, sheet = metadatalist[myinput])
 countdata <- read.xlsx(inputlist[myinput])
 rownames(countdata)<- countdata$Taxonomy
countdata$Taxonomy<- countdata$seq <- NULL
countdata[is.na(countdata)]=0
  #rownames(countdata)<- paste(countdata$Name, countdata$Taxon_ID, sep= "_")
	#countdata$Name<- countdata$Taxon_ID<- countdata$abundance<-NULL
  #countdata<- countdata[2:nrow(countdata),]
	#Condition based on soil type: LC_simpl_2018
	#No covariates
	#We store bigmeta in memory for later re-use

  #row_remove<-unlist(strsplit(removeme,","))
  #for(rown in 1:length(row_remove))
  #{
  #metadata <- metadata[!grepl(row_remove[rown], metadata$BARCODE),]
  #}
  countdata$Taxonomy <- rownames(countdata)
	  browser()
  countdata$Taxonomy <- sub("^[^;]*;", "", countdata$Taxonomy)
  countdata<- aggregate(. ~ Taxonomy, data = countdata, FUN = sum)
  semicolon_counts <- nchar(gsub("[^;]", "", countdata$Taxonomy))
  countdata <- countdata[semicolon_counts >= 5, ]
  rownames(countdata)<- countdata$Taxonomy
  countdata$Taxonomy <- NULL
  bigmeta<-metadata
  options(scipen = 999)
    bigcount <- countdata

	for(bbb in 1:length(allcond))
	{
		# Initialize an empty list to store individual datax data frames

	print(allcond[bbb])
	#if(bbb==2) 
	metadata<-bigmeta
 countdata <- bigcount
  mycondition<-allcond[bbb]
	dir.create(paste0(outputlist[myinput],mycondition), recursive = TRUE) 
	mydir<-paste0(outputlist[myinput],mycondition,"/")
 padj="/padj"
  pdir<-paste0(outputlist[myinput],mycondition,padj,"/")
  dir.create(pdir)
  	metadata<-metadata[,c("Sample_name",mycondition)]
	names(metadata)[2]<-"condition"
	#readcount<-fread(readfile,data.table=F)
 #keep only the samples present in metadata, this allow to input, for instance, only the samples tagged as "italy" in the metadata and perform the analysis only on those samples
  #readcount<-readcount[readcount$V1%in%metadata$Sample_name,]
	countdata<-countdata[,names(countdata)%in%metadata$Sample_name]
	metadata<-metadata[metadata$Sample_name%in%metadata$Sample_name,]
 countdata <- countdata[, mixedsort(names(countdata))]
	rownames(metadata)<-metadata$Sample_name
	#Add (unmapped) read counts to the count table (needed to normalize)
	#Trick to sort the read counts in the same order of the gene counts based on sample names
	#mc<-match(names(countdata),readcount[,1])
	#sreadcount<-readcount[mc,]
	#tcount<-sreadcount[,2]
	#totmapped<-apply(countdata,2,sum)
 	#unmapped<-tcount-totmapped
	#unmapped_countdata<-data.frame(rbind(countdata,unmapped),stringsAsFactors=F)
 	#rownames(unmapped_countdata)[nrow(unmapped_countdata)]<-"Unmapped"
	#Check if metadata has NA and in case assign them to Unknown
	metadata$condition[is.na(metadata$condition)]<-"Unknown"
	#Strip leading and trailing whitespace from condition
 	metadata$condition<-gsub("/","",metadata$condition)
	metadata$condition<-gsub(" ","",metadata$condition)
	#Check that samples are exactly the same and in the same order between counts and metadata
  if(sum(rownames(metadata)==names(countdata))!=ncol(countdata)) stop("Names in countdata do not match names in metadata")
  #countdata <- as.matrix(countdata)
  
  ddsHTSeq <- DESeqDataSetFromMatrix(countData  = countdata,
                                        colData = metadata,
                                        design= ~ condition)
	#Remove low abundance drug classes
      #cat("starting estimateSizeFactors \n")
    #ddsHTSeq<-estimateSizeFactors(ddsHTSeq, type = "poscounts")
    ddsHTSeq<-DESeq(ddsHTSeq)
     write.table(counts(ddsHTSeq,normalized=TRUE),paste(mydir,"norm_counts.txt",sep=""),sep="\t",quote=F,col.names=NA)
    write.table(counts(ddsHTSeq,normalized=FALSE),paste(mydir,"raw_counts.txt",sep=""),sep="\t",quote=F,col.names=NA)
    # #ppm<-counts(ddsHTSeq,normalized=FALSE)
	# #ppm<-ppm*1000000
	# #ppm<-t(apply(ppm,1,"/",tcount))
	# #write.table(ppm,paste(mydir,"ppm.txt",sep=""),sep="\t",quote=F,col.names=NA)
	    #ppm<-counts(ddsHTSeq,normalized=FALSE)
	#ppm<-ppm*1000000
	#ppm<-t(apply(ppm,1,"/",tcount))
	#write.table(ppm,paste(mydir,"ppm.txt",sep=""),sep="\t",quote=F,col.names=NA)
 #unmapped_ddsHTSeq <- DESeqDataSetFromMatrix(countData  = unmapped_countdata,
                                        #colData = metadata,
                                        #design= ~ condition)
    #percentage<-counts(unmapped_ddsHTSeq,normalized=FALSE)
	#percentage<-percentage*100
	#percentage<-t(apply(percentage,1,"/",tcount))
	#write.table(percentage,paste(mydir,"ppm.txt",sep=""),sep="\t",quote=F,col.names=NA)
	vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
	#Perform my own PCA analysis to retrieve the proportion of variance and write it in the axis
	mydat<-data.frame(assay(vsd),check.names=F)
	myvar<-apply(mydat,1,var)
	mydat<-mydat[order(myvar,decreasing=T),]
	#mydat<-mydat[1:500,]
	mypc<-prcomp(t(mydat))
	percentVar <- round(100*mypc$sdev^2/sum(mypc$sdev^2),0)
	myxlab<-paste("PC1: ",percentVar[1],"% variance",sep="")
	myylab<-paste("PC2: ",percentVar[2],"% variance",sep="")
 
 ColorsDT <-  data.table(group=c("Manure","Peat","Sand"), Color=c("#333BFF", "#CC6600", "#9633FF"), key="group")
	pdf(paste0(mydir,"VST_PCCA.pdf"))
	#plotPCA(vsd, intgroup=c("condition"))
	plca<-plotPCA(vsd, intgroup=c("condition"),returnData=T)
	  if(allcond[bbb] == "Soil")
  {
    ann_colors <- c(Sand = "darkgoldenrod1",Manure = "brown",Peat = "red")
 group_colors <- c("Sand" = "Gold3", "Manure" = "Black", "Peat" ="Red")
 } else {
 ann_colors <- c(Yes = "Green",No = "orange")
 group_colors <- c("Yes" = "Green", "No" = "orange")
 }
  graph<-ggplot(plca, aes(x=PC1, y=PC2, color=group),size=3)+geom_point(size=3) + geom_text_repel(aes(label=name),size=4) + coord_fixed() + scale_color_manual(values=group_colors) + theme(legend.text=element_text(size=10)) +xlab(myxlab) + ylab(myylab) 
	print(graph)
	dev.off()
	# pdf("dispersion_fit.pdf")
	# plotDispEsts(dd_1)
	# dev.off()
 	#Build heatmap 
	newcounts<-counts(ddsHTSeq,normalized=TRUE)
	#We don't want to use unmapped, and we set them to zero
	#newcounts=newcounts[rownames(newcounts)!="Unmapped",]
	keepme <- order(rowMeans(newcounts), decreasing=TRUE)[1:min(nrow(newcounts),500)]
	metadata_for_annotations <- read.xlsx(metafile, sheet = "Campionate_edited")
	pcond<-data.frame(Soil=metadata_for_annotations$Soil, Autoclave=metadata_for_annotations$Ster_soil, Heat_root=metadata_for_annotations$Ster_root ,row.names=rownames(metadata_for_annotations))
	#pheatmap(newcounts[keepme,], cluster_rows=FALSE, annotation_col = pcond, fontsize=4, cellwidth=6, cellheight=4, filename=paste0(graphdir,"50-most-abundant-genes_clust.pdf"))
  #Build heatmap on vsd-corrected data
  somma=rowSums(assay(vsd))
  mapvsd=assay(vsd)
  mapvsd=mapvsd[order(somma, decreasing = T),][1:min(nrow(mapvsd),500),]
    pcond_color=pcond
    if(allcond[bbb] == "Soil")
  {
    ann_colors <- list(
  Soil = c(Sand = "Gold3", Manure = "Black", Peat = "Red"),
  Autoclave = c(Yes = "DarkGreen", No = "LightGreen"),
  Heat_root = c(Yes = "DarkBlue", No = "LightBlue")
)
  } else {
  ann_colors <- list(condition =c(Yes = "Green",No = "orange"))
  }
  forwilcox <- head(mapvsd, 35)
  rownames(forwilcox) <- sub(".*;", "", rownames(forwilcox))
perform_wilcoxon <- function(forwilcox, metadata, condition_col, output_file) {
  
  common_samples <- intersect(colnames(forwilcox), metadata$Sample_name)
  
  forwilcox <- forwilcox[, common_samples, drop = FALSE]
  metadata <- metadata %>% filter(Sample_name %in% common_samples)
  
  results_list <- list()
  
  for (taxon in rownames(forwilcox)) {
    taxon_data <- data.frame(
      Abundance = as.numeric(forwilcox[taxon, ]),
      Condition = metadata[[condition_col]]
    )
    
    unique_conditions <- unique(taxon_data$Condition)
    
    if (length(unique_conditions) > 1) {
      p_value <- NA
      posthoc_pvals <- NA
      comparisons <- paste(unique_conditions, collapse = " vs ")
      
      # Wilcoxon test for two groups
      if (length(unique_conditions) == 2) {
        wilcox_res <- wilcox.test(Abundance ~ Condition, data = taxon_data, exact = FALSE)
        p_value <- wilcox_res$p.value
      }
      
      # Kruskal-Wallis + Pairwise Wilcoxon for 3+ groups
      if (length(unique_conditions) > 2) {
        kruskal_res <- kruskal.test(Abundance ~ Condition, data = taxon_data)
        p_value <- kruskal_res$p.value
        
        # Pairwise Wilcoxon test with p-value adjustment
        posthoc_res <- pairwise.wilcox.test(
          taxon_data$Abundance, 
          taxon_data$Condition, 
          p.adjust.method = "BH"  # Adjust p-values using Benjamini-Hochberg (FDR)
        )
        
        # Extract pairwise comparisons and format correctly
        posthoc_pvals <- as.data.frame(as.table(posthoc_res$p.value))
        posthoc_pvals <- posthoc_pvals %>%
          filter(!is.na(Freq)) %>% 
          mutate(Comparison = paste(Var1, "vs", Var2, "=", round(Freq, 5))) %>% 
          pull(Comparison) %>% 
          paste(collapse = "; ")
      }
      
      results_list[[taxon]] <- data.frame(
        Taxon = taxon,
        Comparisons = comparisons,
        P_Value = p_value,
        Posthoc_P_Values = ifelse(length(unique_conditions) > 2, posthoc_pvals, NA)  # Show posthoc only if 3+ groups
      )
    }
  }
  
  final_results <- do.call(rbind, results_list)
  write.table(final_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(final_results)
}

output_file <- paste0(mydir, "/", allcond[bbb], "_35_wilcoxon.txt")
perform_wilcoxon(forwilcox, metadata, "condition", output_file)
pheatmap(mapvsd, cluster_rows=FALSE, annotation_col = pcond, annotation_color = ann_colors, fontsize=4, cellwidth=6, annotation_legend = TRUE, cellheight=4, filename=paste0(mydir,"vsd_50-most-abundant-genes_clust.pdf"))

 # #unmapped_ddsHTSeq <- DESeqDataSetFromMatrix(countData  = unmapped_countdata,
                                        # #colData = metadata,
                                        # #design= ~ condition)
    # #percentage<-counts(unmapped_ddsHTSeq,normalized=FALSE)
	# #percentage<-percentage*100
	# #percentage<-t(apply(percentage,1,"/",tcount))
	# #write.table(percentage,paste(mydir,"ppm.txt",sep=""),sep="\t",quote=F,col.names=NA)
	# vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
    # # write.table(assay(vsd),paste(outputlist[myinput],"VST.txt",sep=""),sep="\t",quote=F,col.names=NA)
	# png(paste(mydir,"PCA.png",sep=""),height=10,width=10,units="cm",res=600,type="cairo")
    # myplot<-plotPCA(vsd,intgroup="condition")
    # print(myplot)
	# dev.off()
	
	# #Perform my own PCA analysis to retrieve the proportion of variance and write it in the axis
	# mydat<-data.frame(assay(vsd),check.names=F)
	# myvar<-apply(mydat,1,var)
	# mydat<-mydat[order(myvar,decreasing=T),]
	# #mydat<-mydat[1:500,]
	# mypc<-prcomp(t(mydat))
	# percentVar <- round(100*mypc$sdev^2/sum(mypc$sdev^2),0)
	# myxlab<-paste("PC1: ",percentVar[1],"% variance",sep="")
	# myylab<-paste("PC2: ",percentVar[2],"% variance",sep="")
 
 # ColorsDT <-  data.table(group=c("Manure","Peat","Sand"), Color=c("#333BFF", "#CC6600", "#9633FF"), key="group")
	# png(paste0(mydir,"VST_PCCA.png"),height=20,width=20,units="cm",res=600,type="cairo")
	# #plotPCA(vsd, intgroup=c("condition"))
	# plca<-plotPCA(vsd, intgroup=c("condition"),returnData=T)
	  # if(allcond[bbb] == "Soil")
  # {
    # ann_colors <- c(Sand = "darkgoldenrod1",Manure = "brown",Peat = "red")
 # group_colors <- c("Sand" = "Gold3", "Manure" = "Black", "Peat" ="Red")
 # } else {
 # ann_colors <- c(Yes = "Green",No = "orange")
 # group_colors <- c("Yes" = "Green", "No" = "orange")
 # }
  # graph<-ggplot(plca, aes(x=PC1, y=PC2, color=group),size=3)+geom_point(size=3) + geom_text_repel(aes(label=name),size=4) + coord_fixed() + scale_color_manual(values=group_colors) + theme(legend.text=element_text(size=10)) +xlab(myxlab) + ylab(myylab) 
	# print(graph)
	# dev.off()
	# # pdf("dispersion_fit.pdf")
	# # plotDispEsts(dd_1)
	# # dev.off()
 	# #Build heatmap 
	# newcounts<-counts(ddsHTSeq,normalized=TRUE)
	# #We don't want to use unmapped, and we set them to zero
	# #newcounts=newcounts[rownames(newcounts)!="Unmapped",]
	# keepme <- order(rowMeans(newcounts), decreasing=TRUE)[1:min(nrow(newcounts),500)]
	# metadata_for_annotations <- read.xlsx(metafile, sheet = "Campionate_edited")
	# pcond<-data.frame(Soil=metadata_for_annotations$Soil, Autoclave=metadata_for_annotations$Ster_soil, Heat_root=metadata_for_annotations$Ster_root ,row.names=rownames(metadata_for_annotations))
	# #pheatmap(newcounts[keepme,], cluster_rows=FALSE, annotation_col = pcond, fontsize=4, cellwidth=6, cellheight=4, filename=paste0(graphdir,"50-most-abundant-genes_clust.pdf"))
  # #Build heatmap on vsd-corrected data
  # somma=rowSums(assay(vsd))
  # mapvsd=assay(vsd)
  # mapvsd=mapvsd[order(somma, decreasing = T),][1:min(nrow(mapvsd),500),]
    # pcond_color=pcond
    # if(allcond[bbb] == "Soil")
  # {
    # ann_colors <- list(
  # Soil = c(Sand = "Gold3", Manure = "Black", Peat = "Red"),
  # Autoclave = c(Yes = "DarkGreen", No = "LightGreen"),
  # Heat_root = c(Yes = "DarkBlue", No = "LightBlue")
# )
  # } else {
  # ann_colors <- list(condition =c(Yes = "Green",No = "orange"))
  # }
# # pheatmap(mapvsd, cluster_rows=FALSE, annotation_col = pcond, annotation_color = ann_colors, fontsize=4, cellwidth=6, annotation_legend = FALSE, cellheight=4, filename=paste0(mydir,"vsd_50-most-abundant-genes_clust.pdf"))

    # #dev.off()
	# #Plot distance matrix as heatmap
	# sampleDist <- dist(t(assay(vsd)))
	# sampleDistMatrix <- as.matrix(sampleDist)
	# #rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
	# rownames(sampleDistMatrix) <- colnames(vsd)
	# colnames(sampleDistMatrix) <- colnames(vsd)
	# colors <- colorRampPalette ( rev( brewer.pal(9 , "Blues")) )( 255)
  # pheatmap(sampleDistMatrix, clustering_distance_rows =sampleDist, annotation_col = pcond, annotation_row = pcond, clustering_distance_cols =sampleDist, filename=paste0(mydir,"VST_samples-distances.pdf"))
	# #dev.off()

	#We are now comparing all vs all
	to.contrast<-unique(metadata$condition)
	compare<-combn(to.contrast,2)
    for(aaa in 1:ncol(compare))
    {
 		res <- results(ddsHTSeq,contrast=c("condition",compare[2,aaa],compare[1,aaa]))
        res<-res[order(res$padj),]
        pres<-subset(res, padj < 0.05)
        res=data.frame(class=row.names(res),res)
        pres=data.frame(class=row.names(pres),pres)
        resfile<-paste(mydir,compare[2,aaa],"_vs_",compare[1,aaa],".txt",sep="")
        presfile<-paste(pdir,compare[2,aaa],"_vs_",compare[1,aaa],".txt",sep="")
		
		# allorno <- c("all","grapevine")
		# for (grapes in allorno)
		# {
				# listben<-read.xlsx(listrelevant, sheet = benef[myinput])
		# listben$References <- listben$Remarks <- NULL
		# if (grapes == "grapevine")
		# {
		# listben <- listben %>%
  # filter(grepl("Vitis|vitis|Grapevine|grapevine", Food.Crops, ignore.case = TRUE))
		# }
  # listben$References <- listben$Remarks <- listben$Food.Crops <- NULL
		# cycleben <- strsplit(as.character(listben$PGPR), ",")
		# cycleben <- unlist(cycleben)
		# cycleben <- sub("\\(.*", "", cycleben)
		
		# listdam<-read.xlsx(listrelevant, sheet = damag[myinput])
		# listdam$References <- listdam$Remarks <- NULL
		# listdam <- listdam %>%
  # filter(grepl("Vitis|vitis|Grapevine|grapevine", Food.Crops, ignore.case = TRUE))
  # listdam$References <- listdam$Remarks <- listdam$Food.Crops <- NULL
		# cycledam <- strsplit(as.character(listdam$PGPR), ",")
		# cycledam <- unlist(cycledam)
		# cycledam <- sub("\\(.*", "", cycledam)		
# # Step 1: Clean up cycleben and cycledam by removing 'spp.' and 'sp.'
# cycleben <- gsub("(_spp\\.|_sp\\.)$", "", cycleben)
# cycleben <- unique(cycleben)
# cycledam <- gsub("(_spp\\.|_sp\\.)$", "", cycledam)
# cycledam <- unique(cycledam)

# # Initialize empty data frames to store the matched rows for both cycleben and cycledam
# matched_ben_df <- data.frame()
# matched_dam_df <- data.frame()
# # Step 2: Find matches for cycleben in pres$class
# for (i in 1:length(cycleben)) {
  # entry <- cycleben[i]
  # print(entry)
  # # Check if the entry contains an underscore
  # if (grepl("_", entry)) {
    # # Replace underscore with semicolon
    # entry_modified <- gsub("_", ";", entry)
    # cat("Modified entry (ben):", entry_modified, "\n")
    
    # # Find matches in pres$class
    # matching_indices <- grep(entry_modified, res$class)
    # print(paste("Matching indices for entry (ben):", entry_modified, "=>", matching_indices))
    
  # } else {
# # Define the search pattern with word boundaries
# pattern <- paste0(";\\b", entry, "\\b")
# matching_indices <- grep(pattern, res$class, fixed = FALSE)  # Set fixed = FALSE to interpret \\b as regex
# print(paste("Matching indices for exact match (dam):", entry, "=>", matching_indices))
  # }
  
  # # Append the matched rows to matched_ben_df if matches are found
  # if (length(matching_indices) > 0) {
    # matched_rows <- res[matching_indices, ]
    # matched_rows$match <- entry  # Add a 'match' column with the entry name
    # matched_ben_df <- rbind(matched_ben_df, matched_rows)
  # } else {
    # print(paste("No match found for entry (ben):", entry))
  # }
# }

# # Step 3: Find matches for cycledam in pres$class
# for (i in 1:length(cycledam)) {
  # entry <- cycledam[i]
  
  # # Check if the entry contains an underscore
  # if (grepl("_", entry)) {
    # # Replace underscore with semicolon
    # entry_modified <- gsub("_", ";", entry)
    # cat("Modified entry (dam):", entry_modified, "\n")
    
    # # Find matches in pres$class
    # matching_indices <- grep(entry_modified, res$class, fixed = TRUE)
    # print(paste("Matching indices for entry (dam):", entry_modified, "=>", matching_indices))
    
  # } else {
    # # No underscore, search for exact match
# # Define the search pattern with word boundaries
# pattern <- paste0(";\\b", entry, "\\b")
# matching_indices <- grep(pattern, res$class, fixed = FALSE)  # Set fixed = FALSE to interpret \\b as regex
# print(paste("Matching indices for exact match (dam):", entry, "=>", matching_indices))
  # }
  
  # # Append the matched rows to matched_dam_df if matches are found
  # if (length(matching_indices) > 0) {
    # matched_rows <- res[matching_indices, ]
    # matched_rows$match <- entry  # Add a 'match' column with the entry name
    # matched_dam_df <- rbind(matched_dam_df, matched_rows)
  # } else {
    # print(paste("No match found for entry (dam):", entry))
  # }
# }


# # Final check: print both resulting data frames
# cat("Matched rows for cycleben:\n")
# print(matched_ben_df)

# if (compare[2, aaa] == "Peat") {
  # colorfirst <- "Red"
# }
# if (compare[1, aaa] == "Peat") {
  # colorsecond <- "Red"
# }

# if (compare[2, aaa] == "Manure") {
  # colorfirst <- "Black"
# }
# if (compare[1, aaa] == "Manure") {
  # colorsecond <- "Black"
# }

# if (compare[2, aaa] == "Sand") {
  # colorfirst <- "Gold3"
# }
# if (compare[1, aaa] == "Sand") {
  # colorsecond <- "Gold3"
# }

# if (compare[2, aaa] == "Yes") {
  # if (allcond[bbb] == "Ster_soil") {
    # colorfirst <- "DarkGreen"
  # } else {
    # colorfirst <- "Darkblue"
  # }
# }

# if (compare[1, aaa] == "No") {
  # if (allcond[bbb] == "Ster_soil") {
    # colorsecond <- "DarkGreen"
  # } else {
    # colorsecond <- "Darkblue"
  # }
# }

# if (compare[2, aaa] == "No") {
# if (allcond[bbb] == "Ster_soil") {
  # colorfirst <- "Lightgreen"
# } else {
  # colorfirst <- "Lightblue"
# }
# }


# if (compare[1, aaa] == "No") {
# if (allcond[bbb] == "Ster_soil") {
  # colorsecond <- "Lightgreen"
# } else {
  # colorsecond <- "Lightblue"
# }
# }
# cat("Matched rows for cycledam:\n")
# print(matched_dam_df)
# print(inputlist[myinput])
# if (nrow(matched_ben_df) != 0 & ncol(matched_ben_df) != 0) {
  # matched_ben_df$color <- with(matched_ben_df, ifelse(padj < 0.05 & log2FoldChange < -1.5, colorsecond,
                              # ifelse(padj < 0.05 & log2FoldChange > 1.5, colorfirst, "#4F4F4F")))
							  # print("zalgo comes")
   # # Calculate the number of differentially abundant ASVs based on padj and log2FoldChange thresholds
  # diff_abundant_asvs <- sum(matched_ben_df$padj < 0.05 & (matched_ben_df$log2FoldChange > 1.5 | matched_ben_df$log2FoldChange < -1.5), na.rm = TRUE)
   # # Calculate the numbers of positive and negative log2FoldChange ASVs
  # pos_diff_asvs <- sum(matched_ben_df$padj < 0.05 & matched_ben_df$log2FoldChange > 1.5, na.rm = TRUE)
  # neg_diff_asvs <- sum(matched_ben_df$padj < 0.05 & matched_ben_df$log2FoldChange < -1.5, na.rm = TRUE)
   # # Create the volcano plot with updated title
  # pdf(paste0(pdir, compare[2,aaa], "_vs_", compare[1,aaa], "_", grapes, "_ben_volcanoplot.pdf"))
  # total_asvs <- nrow(res)
   # volcan <- ggplot(matched_ben_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
    # geom_point(size = 2) +
    # scale_color_identity() +
    # theme_classic() +
    # labs(x = "Log2 Fold Change", y = "-log10(padj)", 
         # title = paste("Volcano Plot (Negative:",compare[1,aaa], " ", neg_diff_asvs, "Positive:",compare[2,aaa], " ", pos_diff_asvs, "Total:", total_asvs, ")")) +
    # geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "grey50") + # Dotted vertical lines for log2FoldChange threshold
    # geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") + # Dotted horizontal line for padj threshold
    # theme(legend.position = "none") # Remove legend if not needed
 # print(volcan)
  # dev.off()
# }
# if (nrow(matched_dam_df) != 0 & ncol(matched_dam_df) != 0) {
  # matched_dam_df$color <- with(matched_dam_df, ifelse(padj < 0.05 & log2FoldChange < -1.5, "#B0BEC5",
                              # ifelse(padj < 0.05 & log2FoldChange > 1.5, "#6A0D91", "#4F4F4F")))
  # # Calculate the number of differentially abundant ASVs based on padj and log2FoldChange thresholds
  # diff_abundant_asvs <- sum(matched_dam_df$padj < 0.05 & (matched_dam_df$log2FoldChange > 1.5 | matched_dam_df$log2FoldChange < -1.5), na.rm = TRUE)
  # # Calculate the numbers of positive and negative log2FoldChange ASVs
  # pos_diff_asvs <- sum(matched_dam_df$padj < 0.05 & matched_dam_df$log2FoldChange > 1.5, na.rm = TRUE)
  # neg_diff_asvs <- sum(matched_dam_df$padj < 0.05 & matched_dam_df$log2FoldChange < -1.5, na.rm = TRUE)
  # # Create the volcano plot with updated title
  # pdf(paste0(pdir, compare[2,aaa], "_vs_", compare[1,aaa],"_", grapes, "_dam_volcanoplot.pdf"))
  # total_asvs <- nrow(res)
  # volcan <- ggplot(matched_dam_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
    # geom_point(size = 2) +
    # scale_color_identity() +
    # theme_classic() +
    # labs(x = "Log2 Fold Change", y = "-log10(padj)", 
         # title = paste("Volcano Plot (Negative:",compare[1,aaa], " ", neg_diff_asvs, "Positive:",compare[2,aaa], " ", pos_diff_asvs, "Total:", total_asvs, ")")) +
    # geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "grey50") + # Dotted vertical lines for log2FoldChange threshold
    # geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") + # Dotted horizontal line for padj threshold
    # theme(legend.position = "none") # Remove legend if not needed
  # print(volcan)
  # dev.off()
		# }
# # Show the matches for both beneficial and damaging lists
		# write.table(matched_dam_df,paste(pdir,compare[2,aaa],"_vs_",compare[1,aaa],"_", grapes,"damaging.txt",sep=""),sep="\t",quote=F,row.names=F)
		# write.table(matched_ben_df,paste(pdir,compare[2,aaa],"_vs_",compare[1,aaa],"_",grapes,"beneficial.txt",sep=""),sep="\t",quote=F,row.names=F)
		
		 write.table(res,resfile,sep="\t",quote=F,row.names=F)
		 write.table(res,resfile,sep="\t",quote=F,row.names=F)
    write.table(pres,presfile,sep="\t",quote=F,row.names=F)
		 mysig<-sum(res$padj<=0.05,na.rm=T)
		 cat(compare[2,aaa],"vs",compare[1,aaa],",",mysig,"significant genes out of",nrow(res),"\n")
	 }
	     for(aaa in 1:ncol(compare))
    {
 		res <- results(ddsHTSeq,contrast=c("condition",compare[2,aaa],compare[1,aaa]))
        res<-res[order(res$padj),]
		rownames(res) <- sub(".*;", "", rownames(res))
		# Ensure row names are set properly
		res <- res[rownames(res) %in% rownames(forwilcox), , drop = FALSE]
        pres<-subset(res, padj < 0.05)
        res=data.frame(class=row.names(res),res)
        pres=data.frame(class=row.names(pres),pres)
        resfile<-paste(mydir,compare[2,aaa],"_vs_",compare[1,aaa],"_35_top_abundant.txt",sep="")
        presfile<-paste(pdir,compare[2,aaa],"_vs_",compare[1,aaa],"_35_top_abundat.txt",sep="")
		write.table(res,resfile,sep="\t",quote=F,row.names=F)
		write.table(pres,presfile,sep="\t",quote=F,row.names=F)
		fortable <- pres
		fortable$baseMean <- fortable$lfcSE	<- fortable$stat	<- fortable$pvalue	<- fortable$padj <- NULL
		fortable$log2FoldChange <- ifelse(fortable$log2FoldChange > 0, 1, 
                                  ifelse(fortable$log2FoldChange < 0, -1, 
                                         fortable$log2FoldChange))
		names(fortable)[2] <- paste0(compare[2,aaa]," vs ",compare[1,aaa])
		mysig<-sum(res$padj<=0.05,na.rm=T)
		 cat(compare[2,aaa],"vs",compare[1,aaa],",",mysig,"significant genes out of",nrow(res),"\n")
		 fortablefile<-paste(pdir,compare[2,aaa],"_vs_",compare[1,aaa],"_35_top_abundant_table.txt",sep="")
		 write.table(fortable,fortablefile,sep="\t",quote=F,row.names=F)
	 }
	 }
	 
# # List all files in the directory
# files <- list.files(pdir, full.names = TRUE)

# # Initialize a list to store names
# names_list <- list()
# # Loop through the files
# for (file in files) {
  # # Check if "grapevine" is in the file name and categorize
   # if (grepl("_vs_", basename(file)) && 
      # (grepl("damaging", basename(file)) || grepl("beneficial", basename(file)))) {

    # # Determine if it's damaging or beneficial
    # if (grepl("damaging", basename(file))) {
      # category <- "damaging"
    # } else {
      # category <- "beneficial"
    # }
    # # Extract names from the file name
    # nameA <- sub("_vs_.*", "", basename(file))
# nameB <- sub(".*_vs_([^_]+).*", "\\1", basename(file))

    # # Store the names in the list
    # names_list[[basename(file)]] <- list(nameA = nameA, nameB = nameB)
# if (grepl("Ster_soil", file))
# {
	# if (nameA == "Yes") {
	# nameA = "Autoclave" 
	# } else {nameA = "No_Autoclave"}
	# if (nameB == "Yes") {
	# nameB = "Autoclave"
	# } else {nameB = "No_Autoclave"}
# }else if (grepl("Ster_root", file))
	# {if (nameA == "Yes") {
	# nameA = "Heat_root"
	# } else {nameA = "No_Heat_root"}
	# if (nameB == "Yes") {
	# nameB = "Heat_root"
	# } else {nameB = "No_Heat_root"}
# }
    # # Load the data from the file (assuming it's a tab-delimited text file)
	# file_contents <- readLines(file)
# if (length(file_contents) == 0 || all(trimws(file_contents) == "")) {
  # print("The file is empty or contains only whitespace.")
# } else {
    # datax <- fread(file, header = TRUE, data.table = F)

# # Add the new column with the dynamic name
# datax <- datax %>%
  # mutate(
    # !!paste0(nameA, "_vs_", nameB) := case_when(
      # padj > 0.05 ~ "ns",
      # log2FoldChange > 0 ~ "1",
      # log2FoldChange < 0 ~ "-1",
      # abs(log2FoldChange) < 0.5 ~ "0"
    # )
  # )

    # # Remove unnecessary columns
    # datax$baseMean <- NULL
    # datax$log2FoldChange <- NULL
    # datax$lfcSE <- NULL
    # datax$stat <- NULL
    # datax$pvalue <- NULL
    # datax$padj <- NULL
    # datax$color <- NULL

# if (category == "damaging") {
    # # Store the processed data in the list
	# if (grepl("all", file)) {
    # if (!grepl(paste0(nameA, "_vs_", nameB), damagio)) {
        # damagio <- merge(damagio, datax, by = c("class", "match"), all = TRUE)
    
# }	else if (grepl("grapevine", file)) {
	    # if (!grepl(paste0(nameA, "_vs_", nameB), damagiogrape)) {
        # damagiogrape <- merge(damagiogrape, datax, by = c("class", "match"), all = TRUE)
		# }
		# }
		# }
# } else if (category == "beneficial") {
 	# if (grepl("all", file)) {
    # if (!grepl(paste0(nameA, "_vs_", nameB), damagio)) {
        # benefo <- merge(benefo, datax, by = c("class", "match"), all = TRUE)
    # }
# }	else if (grepl("grapevine", file)) {
	    # if (!grepl(paste0(nameA, "_vs_", nameB), benefogrape)) {
        # benefogrape <- merge(benefogrape, datax, by = c("class", "match"), all = TRUE)
		# }
		# }
		  
# }
   # }
# }
  # }
  # }
  # # #write the final tables
	# # cycleprint <- c(TRUE,FALSE)
	# # for (mycycle in cycleprint)
	# # {
	# # print("Final tables")
	# # mycycle <- as.logical(mycycle)
		# # printxlsx <- as.data.frame(counts(ddsHTSeq,normalized=mycycle))
		# # printxlsx <- printxlsx %>%
		# # mutate(Gene = rownames(.)) %>%
		# # select(Gene, everything())
		# # if (mycycle == TRUE)
			# # {
			# # write.xlsx(printxlsx, paste0(mydir, "norm_counts.xlsx"), data.table = F)
			# # } else {write.xlsx(printxlsx, paste0(mydir, "raw_counts.xlsx"), data.table = F)
			# # }
	# # # }
	# print("PRINTING")
  # write.table(benefo,paste0(outputlist[myinput],"beneficial.txt"),sep="\t",quote=F,)
  # write.table(damagio,paste0(outputlist[myinput],"damaging.txt"),sep="\t",quote=F,)
   # write.xlsx(benefo,paste0(outputlist[myinput],"beneficial.xlsx"),rownames = F)
  # write.xlsx(damagio,paste0(outputlist[myinput],"damaging.xlsx"),rownames = F)
    # write.table(benefogrape,paste0(outputlist[myinput],"beneficial_grapevine.txt"),sep="\t",quote=F,)
  # write.table(damagiogrape,paste0(outputlist[myinput],"damaging_grapevine.txt"),sep="\t",quote=F,)
   # write.xlsx(benefogrape,paste0(outputlist[myinput],"beneficial_grapevine.xlsx"),rownames = F)
  # write.xlsx(damagiogrape,paste0(outputlist[myinput],"damaging_grapevine.xlsx"),rownames = F)

  }
  }
runDEseq()
