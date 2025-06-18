#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
    make_option(c("-A", "--abundance"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/_Species_merged_count.txt",
                help="Path to abundance 16s T1 file [default= %default]", metavar="character"),
	    make_option(c("-B", "--abundanceT0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/_Species_merged_count.txt",
                help="Path to abundance 16s T1 file [default= %default]", metavar="character"),
    make_option(c("-C", "--abundanceITS"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/_Species_merged_count.txt",
                help="Path to abundance ITS T0 file [default= %default]", metavar="character"),
	    make_option(c("-D", "--abundanceITST0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/_Species_merged_count.txt",
                help="Path to abundance ITS T0 file [default= %default]", metavar="character"),
	     make_option(c("-E", "--output"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/07_plot_SILVA/",
                help="Output directory for 16s T0 results [default= %default]", metavar="character"),
     make_option(c("-F", "--outputT0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/07_plot_SILVA/",
                help="Output directory for 16s T0 results [default= %default]", metavar="character"),
    make_option(c("-G", "--outputITS"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/07_plot_unite/",
                help="Output directory for 16s T1 results [default= %default]", metavar="character"),
	    make_option(c("-H", "--outputITST0"), type="character",
                default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/07_plot_unite/",
                help="Output directory for 16s T1 results [default= %default]", metavar="character")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-A' flag.")
} else {  cat ("abundancephylum is ", opt$abundance, "\n")
  a16sT1 <- opt$abundance
  }
  
 if (is.null(opt$abundanceT0)) {
  stop("WARNING: No abundanceT0 specified with '-B' flag.")
} else {  cat ("abundancephylum is ", opt$abundanceT0, "\n")
  a16sT0 <- opt$abundanceT0
  } 
  
if (is.null(opt$abundanceITS)) {
  stop("WARNING: No abundanceITS specified with '-C' flag.")
} else {  cat ("abundanceITS is ", opt$abundanceITS, "\n")
  aITST1 <- opt$abundanceITS  
  }

if (is.null(opt$abundanceITST0)) {
  stop("WARNING: No abundanceITST0 specified with '-D' flag.")
} else {  cat ("abundanceITS is ", opt$abundanceITST0, "\n")
  aITST0 <- opt$abundanceITST0  
  }

if (is.null(opt$output)) {
  stop("WARNING: No 1 specified with '-E' flag.")
} else {  cat ("output is ", opt$output, "\n")
  o16sT1 <- opt$output  
  }

if (is.null(opt$outputT0)) {
  stop("WARNING: No 1 specified with '-F' flag.")
} else {  cat ("outputT0 is ", opt$outputT0, "\n")
  o16sT0 <- opt$outputT0  
  }

  if (is.null(opt$outputITS)) {
  stop("WARNING: No outputITS directory specified with '-G' flag.")
} else {  cat ("outputITS dir is ", opt$outputITS, "\n")
  oITST1 <- opt$outputITS 
  }

  if (is.null(opt$outputITST0)) {
  stop("WARNING: No outputITST0 directory specified with '-H' flag.")
} else {  cat ("outputITST0 dir is ", opt$outputITST0, "\n")
  oITST0 <- opt$outputITST0 
  }




# Load necessary libraries
rarefact <- function() 
{
  library(ggplot2)
  library(vegan) # for rarecurve function
  library(tidyr)
  library(tibble)
  library(openxlsx)
  library("RColorBrewer")
  library(viridis)
  options(scipen = 999)
  
  # Input and output lists
  listinput <- c(a16sT0, a16sT1, aITST0, aITST1)
  listoutput <- c(o16sT0, o16sT1, oITST0, oITST1)
  
  # Loop through each input/output pair
  for (i in 1:length(listinput)) {
    myinput <- listinput[i]
    myoutput <- listoutput[i]

    # Load the abundance data
    abundance_data <- read.table(myinput, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Get the column names
    col_names <- colnames(abundance_data)

    # Use a regular expression to identify columns that contain at least one digit (0-9)
    cols_with_numbers <- grep("\\d", col_names, value = TRUE)

    # Subset the data.frame to keep only columns with numbers in their names
    abundance_data <- abundance_data[, cols_with_numbers]

    # Transpose data to have species as rows and samples as columns
    abundance_data <- t(abundance_data)
# Remove "X" from rownames of abundance_data
rownames(abundance_data) <- gsub("^X", "", rownames(abundance_data))
    # Plot using ggplot2
    output_file <- paste0(myoutput, "rarefaction_curve_", basename(myinput), ".pdf")
    pdf(output_file)
aaa<- rare_curve_data <- rarecurve(abundance_data, step = 20, label = TRUE, col = "blue")
print(aaa)

    dev.off()
    
    cat("Rarefaction curve for", myinput, "saved to", output_file, "\n")
  }
}

# Call the function to generate rarefaction curves
rarefact()

# Call the function to generate rarefaction curves
rarefact()