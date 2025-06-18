#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
    make_option(c("-A", "--abundanceT1"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1_new/06_tables_silva/_Species_merged_count.txt",
                help="Path to abundance 16s T1 file [default= %default]", metavar="character"),
	    make_option(c("-B", "--abundanceT0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/_Species_merged_count.txt",
                help="Path to abundance 16s T1 file [default= %default]", metavar="character"),
    make_option(c("-C", "--abundanceITST1"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/_Species_merged_count.txt",
                help="Path to abundance ITS T0 file [default= %default]", metavar="character"),
	    make_option(c("-D", "--abundanceITST0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/_Species_merged_count.txt",
                help="Path to abundance ITS T0 file [default= %default]", metavar="character"),
	     make_option(c("-E", "--output"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1_new/06_tables_silva/",
                help="Output directory for 16s T0 results [default= %default]", metavar="character"),
     make_option(c("-F", "--outputT0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/",
                help="Output directory for 16s T0 results [default= %default]", metavar="character"),
    make_option(c("-G", "--outputITS"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/",
                help="Output directory for 16s T1 results [default= %default]", metavar="character"),
	    make_option(c("-H", "--outputITST0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/",
                help="Output directory for 16s T1 results [default= %default]", metavar="character"),
	     make_option(c("-I", "--countfragments"), type="character",
                default="/projects/marroni/intevine/docs/count_reads_16s_output_T1.txt",
                help="Path to count fragments T0 file [default= %default]", metavar="character"),
     make_option(c("-L", "--countfragmentsT0"), type="character",
                default="/projects/marroni/intevine/docs/count_reads_16s_T0_output.txt",
                help="Path to count fragments T0 file [default= %default]", metavar="character"),
    make_option(c("-M", "--countfragmentsITS"), type="character",
                default="/projects/marroni/intevine/docs/count_reads_ITS_output_T1.txt",
                help="Path to count fragments T1 file [default= %default]", metavar="character"),
	    make_option(c("-N", "--countfragmentsITST0"), type="character",
                default="/projects/marroni/intevine/docs/count_reads_ITS_T0_output.txt",
                help="Path to count fragments T1 file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundanceT1)) {
  stop("WARNING: No abundanceT1 specified with '-A' flag.")
} else {  
  cat("abundanceT1 is ", opt$abundanceT1, "\n")
  a16sT1 <- opt$abundanceT1
}

if (is.null(opt$abundanceT0)) {
  stop("WARNING: No abundanceT0 specified with '-B' flag.")
} else {  
  cat("abundanceT0 is ", opt$abundanceT0, "\n")
  a16sT0 <- opt$abundanceT0
} 

if (is.null(opt$abundanceITST1)) {
  stop("WARNING: No abundanceITST1 specified with '-C' flag.")
} else {  
  cat("abundanceITS is ", opt$abundanceITST1, "\n")
  aITST1 <- opt$abundanceITST1  
}

if (is.null(opt$abundanceITST0)) {
  stop("WARNING: No abundanceITST0 specified with '-D' flag.")
} else {  
  cat("abundanceITST0 is ", opt$abundanceITST0, "\n")
  aITST0 <- opt$abundanceITST0  
}

if (is.null(opt$output)) {
  stop("WARNING: No output specified with '-E' flag.")
} else {  
  cat("output is ", opt$output, "\n")
  o16sT1 <- opt$output  
}

if (is.null(opt$outputT0)) {
  stop("WARNING: No outputT0 specified with '-F' flag.")
} else {  
  cat("outputT0 is ", opt$outputT0, "\n")
  o16sT0 <- opt$outputT0  
}

if (is.null(opt$outputITS)) {
  stop("WARNING: No outputITS directory specified with '-G' flag.")
} else {  
  cat("outputITS dir is ", opt$outputITS, "\n")
  oITST1 <- opt$outputITS 
}

if (is.null(opt$outputITST0)) {
  stop("WARNING: No outputITST0 directory specified with '-H' flag.")
} else {  
  cat("outputITST0 dir is ", opt$outputITST0, "\n")
  oITST0 <- opt$outputITST0 
}

if (is.null(opt$countfragments)) {
  stop("WARNING: No countfragments specified with '-I' flag.")
} else {  
  cat("countfragments is ", opt$countfragments, "\n")
  c16sT1 <- opt$countfragments 
}

if (is.null(opt$countfragmentsT0)) {
  stop("WARNING: No countfragmentsT0 specified with '-L' flag.")
} else {  
  cat("countfragmentsT0 is ", opt$countfragmentsT0, "\n")
  c16sT0 <- opt$countfragmentsT0 
}

if (is.null(opt$countfragmentsITS)) {
  stop("WARNING: No countfragmentsITS specified with '-M' flag.")
} else {  
  cat("countfragmentsITS is ", opt$countfragmentsITS, "\n")
  cITST1 <- opt$countfragmentsITS 
}

if (is.null(opt$countfragmentsITST0)) {
  stop("WARNING: No countfragmentsITST0 specified with '-N' flag.")
} else {  
  cat("countfragmentsITST0 is ", opt$countfragmentsITST0, "\n")
  cITST0 <- opt$countfragmentsITST0 
}

percentage<-function() 
{
	library(corrplot)
    library(data.table)
	library(openxlsx)
	library("RColorBrewer")
	library("ggplot2")
	library(viridis)
	library(tidyr)
	library(tibble)
	options(scipen=999)
  listinput <- c(a16sT0,a16sT1,aITST0,aITST1)
  listoutput <- c(o16sT0,o16sT1,oITST0,oITST1)
  listreads <- c(c16sT0,c16sT1,cITST0,cITST1)
  for(myinput in 1:length(listinput))
  {
  browser()
  dir.create(listoutput[myinput], recursive = TRUE)
	countdata<-fread(listinput[myinput], data.table = F, header=TRUE)
	print(listoutput[myinput])
	counts <- fread(listreads[myinput], data.table = F)
  numeric_col_names <- grep("^[0-9]+$", names(countdata), value = TRUE)
  numeric_values <- as.numeric(numeric_col_names)
  samplesnumber <- max(numeric_values)
  # Extract numeric and non-numeric column names
  numeric_names <- names(countdata)[which(names(countdata) %in% as.character(1:36))]
  non_numeric_names <- names(countdata)[!names(countdata) %in% as.character(1:36)]
  # Sort numeric names numerically
  numeric_names_sorted <- sort(as.numeric(numeric_names))
  # Combine sorted numeric names with non-numeric names
  new_names <- c("seq", as.character(numeric_names_sorted), non_numeric_names)
  # Reorder columns in dataframe based on new_names
  countdata <- countdata[, new_names]
  countdata$seq.1 <- NULL
# Check the new order of names
print(names(countdata))
	are_values_in_same_order <- all(counts$V1 == names(countdata[2:(samplesnumber+1)]))
	if (are_values_in_same_order) {
  print("Values in 'counts' column are in the same order as 'countdata' names.")
	} else {
	print("Values in 'counts' column are not in the same order as 'countdata' names.")
	stop("Execution stopped: counts are not in the same order as countdata names.")
	}
  numeric_col_names <- grep("^[0-9]+$", names(countdata), value = TRUE)
  numeric_values <- as.numeric(numeric_col_names)
  samplesnumber <- max(numeric_values)
  #for now i want to keep the column with the sequences which is in position 1, so i add +1 to the start and end of the cycle
 for (aaa in 2:(samplesnumber+1))
  {
   print (aaa)
  countdata[aaa] <- (countdata[,aaa]/counts$V2[aaa-1])*100
  }
write.xlsx(countdata, paste0(listoutput[myinput], "percentage_counts.xlsx"), row.names = F)
write.table(countdata, paste0(listoutput[myinput], "percentage_counts.txt"), row.names = F,quote = F)
countdata_editing <- countdata
countdata_editing$Species <- NULL
countdata_editing$Taxonomy <- paste(
  # This will create "ASV_" followed by the row number
  countdata_editing$Kingdom,
  countdata_editing$Phylum,
  countdata_editing$Class,
  countdata_editing$Order,
  countdata_editing$Family,
  countdata_editing$Genus,
  countdata_editing$Species,
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
countdata_editing$Taxonomy <- sapply(countdata_editing$Taxonomy, clean_taxonomy)
countdata_editing$Kingdom<- countdata_editing$Phylum<- countdata_editing$Class<- countdata_editing$Order<- countdata_editing$Family<- countdata_editing$Genus<- countdata_editing$Species <- NULL
# Move the Taxonomy column to the beginning of the dataframe
countdata_editing <- countdata_editing[, c("Taxonomy", setdiff(names(countdata_editing), "Taxonomy"))]
write.xlsx(countdata_editing, paste0(listoutput[myinput], "percentage_counts_edited.xlsx"), row.names = F)
write.table(countdata_editing, paste0(listoutput[myinput], "percentage_counts_edited.txt"), row.names = F,quote = F)
##########remove the ASV which don't have at least 0.1 abundance in at least one sample
#numeric_columns <- sapply(countdata, is.numeric)
#numeric_column_names <- names(numeric_columns)[numeric_columns]
#countdata <- countdata[apply(countdata[, numeric_columns], 1, function(row) any(row >= 0.001)), ]
# Identify the numeric columns
numeric_columns <- sapply(countdata, is.numeric)
row_means <- apply(countdata[, numeric_columns], 1, mean)
countdata <- countdata[row_means >= 0.001, ]
##########remove mitochondria, non bacteria, chloroplast
countdata <- countdata[!apply(countdata, 1, function(row) any(grepl("Mitochondria|Chloroplast", row))), ]
countdata <- countdata[grepl("Bacteria|Fungi", countdata$Kingdom), ]
countdata <- as.data.frame(lapply(countdata, function(x) gsub("[a-z]__", "", x)))
countdata <- as.data.frame(countdata[!apply(countdata, 1, function(row) any(grepl("Incertae", row))), ])
colnames(countdata) <- gsub("^X", "", colnames(countdata))
# Convert columns marked as TRUE in numeric_columns to numeric
countdata[] <- lapply(seq_along(countdata), function(i) {
  if (numeric_columns[i]) {
    as.numeric(as.character(countdata[[i]]))
  } else {
    countdata[[i]]
  }
})
write.xlsx(countdata, paste0(listoutput[myinput], "0.001_filtered_percentage_counts.xlsx"), row.names = F)
write.table(countdata, paste0(listoutput[myinput], "0.001_filtered_percentage_counts.txt"), row.names = F,quote = F)
# Add a new column 'Taxonomy' to the data frame, remove species since we don't have that resolution
countdata_species <- countdata
countdata_species$Taxonomy <- paste(
  paste0("ASV_", seq_len(nrow(countdata_species))),   # This will create "ASV_" followed by the row number
  countdata_species$Kingdom,
  countdata_species$Phylum,
  countdata_species$Class,
  countdata_species$Order,
  countdata_species$Family,
  countdata_species$Genus,
  countdata_species$Species,
  sep = ";"
)

countdata_species$Taxonomy <- sapply(countdata_species$Taxonomy, clean_taxonomy)
countdata_species$Kingdom<- countdata_species$Phylum<- countdata_species$Class<- countdata_species$Order<- countdata_species$Family<- countdata_species$Genus<- countdata_species$Species <- NULL
# Move the Taxonomy column to the beginning of the dataframe
countdata_species <- countdata_species[, c("Taxonomy", setdiff(names(countdata_species), "Taxonomy"))]
write.xlsx(countdata_species, paste0(listoutput[myinput], "0.001_average_filtered_percentage_counts_edited_with_seq_species.xlsx"), row.names = F)
write.table(countdata_species, paste0(listoutput[myinput], "0.001_average_filtered_percentage_counts_edited_with_seq_species.txt"), row.names = F,quote = F)
countdata$Species <- NULL
countdata$Taxonomy <- paste(
  paste0("ASV_", seq_len(nrow(countdata))),   # This will create "ASV_" followed by the row number
  countdata$Kingdom,
  countdata$Phylum,
  countdata$Class,
  countdata$Order,
  countdata$Family,
  countdata$Genus,
  countdata$Species,
  sep = ";"
)

countdata$Taxonomy <- sapply(countdata$Taxonomy, clean_taxonomy)
countdata$Kingdom<- countdata$Phylum<- countdata$Class<- countdata$Order<- countdata$Family<- countdata$Genus<- countdata$Species <- NULL
# Move the Taxonomy column to the beginning of the dataframe
countdata <- countdata[, c("Taxonomy", setdiff(names(countdata), "Taxonomy"))]
write.xlsx(countdata, paste0(listoutput[myinput], "0.001_average_filtered_percentage_counts_edited_with_seq.xlsx"), row.names = F)
write.table(countdata, paste0(listoutput[myinput], "0.001_average_filtered_percentage_counts_edited_with_seq.txt"), row.names = F,quote = F)
countdata$seq <- NULL
write.xlsx(countdata, paste0(listoutput[myinput], "0.001_average_filtered_percentage_counts_edited.xlsx"), row.names = F)
write.table(countdata, paste0(listoutput[myinput], "0.001_average_filtered_percentage_counts_edited.txt"), row.names = F,quote = F)
}	
}
percentage()
