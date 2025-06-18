#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/0.001_average_filtered_percentage_counts_edited_with_seq.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/0.001_average_filtered_percentage_counts_edited_with_seq.xlsx",
              help="Input directory [default= %default]", metavar="character"),
	make_option(c("-A", "--abundanceITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/0.001_average_filtered_percentage_counts_edited_with_seq.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/0.001_average_filtered_percentage_counts_edited_with_seq.xlsx",
              help="Input directory [default= %default]", metavar="character"),
   make_option(c("-A", "--raw_allT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/_Species_merged_count.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--raw_allT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/_Species_merged_count.txt",
              help="Input directory [default= %default]", metavar="character"),
	make_option(c("-A", "--raw_allITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/_Species_merged_count.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--raw_allITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/_Species_merged_count.txt",
              help="Input directory [default= %default]", metavar="character"),
	 make_option(c("-O", "--outtableT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/", 
              help="output dir [default= %default]", metavar="character"),
    make_option(c("-O", "--outtableT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/", 
              help="output dir [default= %default]", metavar="character"),
	make_option(c("-O", "--outtableITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/", 
              help="output dir [default= %default]", metavar="character"),
    make_option(c("-O", "--outtableITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/", 
              help="output dir [default= %default]", metavar="character")
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

if (is.null(opt$abundanceITST0)) {
  stop("WARNING: No abundanceITST0 specified with '-C' flag.")
} else {  cat ("abundanceITST0 is ", opt$abundanceITST0, "\n")
  aITST0 <- opt$abundanceITST0  
  }
  
if (is.null(opt$abundanceITST1)) {
  stop("WARNING: No abundanceITST1 specified with '-D' flag.")
} else {  cat ("abundanceITST1 is ", opt$abundanceITST1, "\n")
  aITST1 <- opt$abundanceITST1  
  }

if (is.null(opt$raw_allT0)) {
  stop("WARNING: No raw_allT0 specified with '-A' flag.")
} else {  cat ("raw_allT0 is ", opt$raw_allT0, "\n")
  rawT0 <- opt$raw_allT0  
  }

if (is.null(opt$raw_allT1)) {
  stop("WARNING: No raw_allT1 specified with '-B' flag.")
} else {  cat ("raw_allT1 is ", opt$raw_allT1, "\n")
  rawT1 <- opt$raw_allT1  
  }

if (is.null(opt$raw_allITST0)) {
  stop("WARNING: No raw_allITST0 specified with '-C' flag.")
} else {  cat ("raw_allITST0 is ", opt$raw_allITST0, "\n")
  rawITST0 <- opt$raw_allITST0  
  }
  
if (is.null(opt$raw_allITST1)) {
  stop("WARNING: No raw_allITST1 specified with '-D' flag.")
} else {  cat ("raw_allITST1 is ", opt$raw_allITST1, "\n")
  rawITST1 <- opt$raw_allITST1  
  }

    if (is.null(opt$outtableT0)) {
  stop("WARNING: No outtableT0 specified with '-I' flag.")
} else {  cat ("outtableT0 dir is ", opt$outtableT0, "\n")
  otableT0 <- opt$outtableT0  
  }
  
    if (is.null(opt$outtableT1)) {
  stop("WARNING: No outtableT1 specified with '-I' flag.")
} else {  cat ("outtableT1 dir is ", opt$outtableT1, "\n")
  otableT1 <- opt$outtableT1  
  }
  
  if (is.null(opt$outtableITST0)) {
  stop("WARNING: No outtableITST0 specified with '-I' flag.")
} else {  cat ("outtableITST0 dir is ", opt$outtableITST0, "\n")
  otableITST0 <- opt$outtableITST0  
  }
  
    if (is.null(opt$outtableITST1)) {
  stop("WARNING: No outtableITST1 specified with '-I' flag.")
} else {  cat ("outtableITST1 dir is ", opt$outtableITST1, "\n")
  otableITST1 <- opt$outtableITST1 
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
  inputlist <- c(aT0,aT1,aITST0,aITST1)
  rawlist <- c(rawT0,rawT1,rawITST0,rawITST1)
     tablelist <- c(otableT0, otableT1,otableITST0,otableITST1)
  for(myinput in 1:length(inputlist))
  {
	countdata<-read.xlsx(inputlist[myinput])
	rawdata<-fread(rawlist[myinput], data.table = FALSE, header = TRUE)
	numeric_columns <- sapply(rawdata, is.numeric)
	rawdata <- as.data.frame(lapply(rawdata, function(x) gsub("[a-z]__", "", x)))
rawdata <- as.data.frame(rawdata[!apply(rawdata, 1, function(row) any(grepl("Incertae", row))), ])
colnames(rawdata) <- gsub("^X", "", colnames(rawdata))
# Convert columns marked as TRUE in numeric_columns to numeric
rawdata[] <- lapply(seq_along(rawdata), function(i) {
  if (numeric_columns[i]) {
    as.numeric(as.character(rawdata[[i]]))
  } else {
    rawdata[[i]]
  }
})
	filtered_rawdata <- rawdata[rawdata$seq %in% countdata$seq, ]
	filtered_rawdata$Species <- NULL
filtered_rawdata$Taxonomy <- paste(
  paste0("ASV_", seq_len(nrow(filtered_rawdata))),   # This will create "ASV_" followed by the row number
  filtered_rawdata$Kingdom,
  filtered_rawdata$Phylum,
  filtered_rawdata$Class,
  filtered_rawdata$Order,
  filtered_rawdata$Family,
  filtered_rawdata$Genus,
  filtered_rawdata$Species,
  sep = ";"
) 
# Create a helper function to clean up the Taxonomy string
clean_taxonomy <- function(taxonomy) {
  # Remove "NA" and trailing semicolons
  cleaned <- gsub(";?NA;?", "", taxonomy)
  # Remove any trailing semicolon left after the above operation
  cleaned <- gsub(";$", "", cleaned)
  return(cleaned)
}
filtered_rawdata$Taxonomy <- sapply(filtered_rawdata$Taxonomy, clean_taxonomy)
filtered_rawdata$Kingdom<- filtered_rawdata$Phylum<- filtered_rawdata$Class<- filtered_rawdata$Order<- filtered_rawdata$Family<- filtered_rawdata$Genus<- filtered_rawdata$Species <- NULL
 filtered_rawdata <- filtered_rawdata[, mixedsort(names(filtered_rawdata))]
# Move the Taxonomy column to the beginning of the dataframe
filtered_rawdata <- filtered_rawdata[, c("Taxonomy", setdiff(names(filtered_rawdata), "Taxonomy"))]
	write.xlsx(filtered_rawdata, paste0(tablelist[myinput], "0.001_raw_counts_edited.xlsx"), row.names = F)
    }
	}     
	stacked()

