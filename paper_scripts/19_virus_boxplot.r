#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT1"), type="character", default="/projects/marroni/intevine/alignment_vv_virus/tables/RNAseq_VirDet.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/projects/marroni/intevine/docs/map_grapevine_final.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-P", "--outputT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/07_plot_SILVA/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-C", "--conditionlistT1"), type="character", default="Soil,Ster_soil,Ster_root", 
              help="conditionlist [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")
 
if (is.null(opt$abundanceT1)) {
  stop("WARNING: No abundanceT1 specified with '-B' flag.")
} else {  cat ("abundanceT1 is ", opt$abundanceT1, "\n")
  aT1 <- opt$abundanceT1  
  }  
  
if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-V' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }
  
	if (is.null(opt$conditionlistT1)) {
  stop("WARNING: No conditionlistT1 specified with '-C' flag.")
} else {  cat ("conditionlist is ", opt$conditionlistT1, "\n")
  conT1 <- opt$conditionlistT1  
  }   
  
     if (is.null(opt$outputT1)) {
  stop("WARNING: No outputT1 specified with '-I' flag.")
} else {  cat ("outputT1 dir is ", opt$outputT1, "\n")
  oT1 <- opt$outputT1  
  }


indices <- function() {
  library(corrplot)
  library(data.table)
  library(openxlsx)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(dplyr)
  library(janitor)
  library(vegan)
  library(ggpubr)
  
  inputlist <- c(aT1)
  outputlist <- c(oT1)
  namelist <- c("species")
  
  for (myinput in 1:length(inputlist)) {
    outdir <- outputlist[myinput]
    countdata <- read.xlsx(inputlist[myinput])
    row.names(countdata) <- countdata$Name
    countdata$Name <- NULL
    metadata <- read.xlsx(metafile, sheet = "Campionate_edited")
    allcond <- unlist(strsplit(conT1, ","))
    bigmeta <- metadata
    bigcount <- countdata
    
    for (bbb in 1:length(allcond)) {
      metadata <- bigmeta
      countdata <- bigcount
      condition <- allcond[bbb]
      metadata <- metadata[, c("Sample_name", condition)]
      names(metadata)[2] <- "condition"
      metadata <- metadata[metadata$Sample_name %in% names(countdata), ]
      rownames(metadata) <- metadata$Sample_name
      countdata <- data.frame(t(countdata))
      column_sums <- colSums(countdata)
      columns_to_keep <- column_sums != 0
      countdata <- countdata[, columns_to_keep]
      # Soil settings based on 'allcond' and 'inputlist'
      if (allcond[bbb] == "Ster_soil")
	  {
        my_comparisons <- list(c("Yes", "No"))
        FillTrattamento <- c("DarkGreen", "LightGreen")
		pingas <- "Autoclave"
		} else if ( allcond[bbb] == "Ster_root") {
		 my_comparisons <- list(c("Yes", "No"))
        FillTrattamento <- c("DarkBlue", "LightBlue")
		pingas <- "Heat root"
	      } else if (allcond[bbb] == "Soil" ) {
        my_comparisons <- list(
          c("Peat", "Manure"),
          c("Peat", "Sand"),
          c("Manure", "Sand")
        )
        FillTrattamento <- c("#f5d902", "#f50202", "#070808")
		pingas <- "Soil"
      }
      
      # Merge data frames
      countmeta <- merge(countdata, metadata, by = "row.names")
	  if (allcond[bbb] == "Soil" ) {
      countmeta$condition <- factor(countmeta$condition, levels = c('Sand', 'Peat', 'Manure'), ordered = TRUE)
	  }
      
      # Output directory and file names
      outfile_base <- paste0(outdir, condition, "_", namelist, "_virus_boxplot")
      
      # Function to format y-axis labels
      format_y_label <- function(label) {
        formatted_label <- gsub("\\.", " ", label)
        return(paste0(formatted_label, " (FPKM)"))
      }
      
      # Plot without Kruskal
      pdf(paste0(outfile_base, "_no_kruskal.pdf"))
      for (i in 2:7) {
        y_label <- format_y_label(colnames(countmeta)[i])
        graphi1 <- ggboxplot(countmeta, x = "condition", y = colnames(countmeta)[i],
                             color = "condition", palette = FillTrattamento) +
          stat_compare_means(comparisons = my_comparisons) +
          ylab(y_label) +
  xlab(pingas)
        print(graphi1)
      }
      dev.off()
      
      # Plot with Kruskal
      pdf(paste0(outfile_base, ".pdf"))
      for (i in 2:7) {
        y_label <- format_y_label(colnames(countmeta)[i])
        max_y <- max(countmeta[[colnames(countmeta)[i]]], na.rm = TRUE)
        graphi2 <- ggboxplot(countmeta, x = "condition", y = colnames(countmeta)[i],
                             color = "condition", palette = FillTrattamento) +
          stat_compare_means(comparisons = my_comparisons) +
          stat_compare_means(label.y = max_y) +
          ylab(y_label) +
  xlab(pingas)
        print(graphi2)
      }
      dev.off()
    }
  }
}
	
indices()