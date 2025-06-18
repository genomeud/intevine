#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/genus_aggregated.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/genus_aggregated.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1species"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/0.001_filtered_percentage_counts.xlsx",
				help="Input directory [default= %default]", metavar="character"),
  make_option(c("-O", "--outputT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/", 
              help="output dir [default= %default]", metavar="character"),
    make_option(c("-P", "--outputT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-G", "--genescorr"), type="numeric", default=25, 
              help="Comma separate list of samples to be removed, e.g. because of low yield [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundanceT0)) {
  stop("WARNING: No abundanceT0 specified with '-A' flag.")
} else {  cat ("abundanceT0 is ", opt$abundanceT0, "\n")
  aT0 <- opt$abundanceT0  
  }

if (is.null(opt$abundanceT1)) {
  stop("WARNING: No abundanceT1 specified with '-B' flag.")
} else {  cat ("abundanceT1 is ", opt$abundanceT1, "\n")
  aT1 <- opt$abundanceT1  
  }
  
 if (is.null(opt$abundanceT1species)) {
  stop("WARNING: No abundanceT1species specified with '-B' flag.")
} else {  cat ("abundanceT1species is ", opt$abundanceT1species, "\n")
  aT1species <- opt$abundanceT1species  
  }
  
if (is.null(opt$genescorr)) {
  stop("WARNING: No genescorr specified with '-V' flag.")
} else {  cat ("genescorr is ", opt$genescorr, "\n")
  genescorr <- opt$genescorr  
  }

  if (is.null(opt$outputT0)) {
  stop("WARNING: No outputT0 specified with '-I' flag.")
} else {  cat ("outputT0 dir is ", opt$outputT0, "\n")
  oT0 <- opt$outputT0  
  }
  
    if (is.null(opt$outputT1)) {
  stop("WARNING: No outputT1 specified with '-I' flag.")
} else {  cat ("outputT1 dir is ", opt$outputT1, "\n")
  oT1 <- opt$outputT1  
  }


stacked<-function() 
{
	library(corrplot)
    library(data.table)
	library(openxlsx)
	library("RColorBrewer")
	library("ggplot2")
		library("ggrepel")
 library(tidyr)
 library(dplyr)
 library(janitor)
	library(vegan)
  inputlist <- c(aT0,aT1,aT1species)
  outputlist <- c(oT0,oT1,oT1)
  for(myinput in 1:length(inputlist))
  {
  	dir.create(outputlist[myinput], recursive = TRUE)
	if (inputlist[myinput] == aT1species)
		{
		countdata<-read.xlsx(inputlist[myinput])
		countdata <- countdata[complete.cases(countdata$Species), ]
		countdata$Species <- paste0(countdata$Genus,"_",countdata$Species)
		countdata$seq <- countdata$Kingdom <- countdata$Phylum <- countdata$Class <- countdata$Order <- countdata$Family <- countdata$Genus <- NULL
		write.table(countdata,paste0(outputlist[myinput], "001_only_species.txt"), row.names = F, quote = F, sep = "\t")
		} else {
		countdata<-read.xlsx(inputlist[myinput])
	countdata$Taxonomy <- sub("^.*;([^;]+)$", "\\1", countdata$Taxonomy)
	write.table(countdata,paste0(outputlist[myinput], "only_genus_aggregated.txt"), row.names = F, quote = F, sep = "\t")
    }
	}
	}
    stacked()

