# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni (adapted for current script)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tidyverse)
})

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default="/projects/marroni/intevine/STAR/",
              help="Directory containing 'ReadsPerGene' files [default= %default]", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="rna_table.txt",
              help="Output filename for the combined RNA table [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check if input directory exists
if (!dir.exists(opt$input_dir)) {
  stop(paste("Error: Input directory does not exist:", opt$input_dir))
}

# Get list of 'ReadsPerGene' files
outputs <- list.files(opt$input_dir, pattern="ReadsPerGene", full.names = TRUE)

if (length(outputs) == 0) {
  stop(paste("Error: No 'ReadsPerGene' files found in", opt$input_dir))
}

# Initialize the combined data table with the first file
one <- fread(outputs[1], data.table = FALSE, header = FALSE)
one$V2 <- one$V4 <- NULL
colnames(one) <- c("Gene", basename(outputs[1])) # Use basename for cleaner column names

# Loop through the rest of the files and merge them
for (f_idx in 2:length(outputs)) {
  my_cycle <- fread(outputs[f_idx], data.table = FALSE, header = FALSE)
  my_cycle$V2 <- my_cycle$V4 <- NULL
  colnames(my_cycle) <- c("Gene", basename(outputs[f_idx])) # Use basename
  one <- merge(one, my_cycle, by = "Gene", all = TRUE)
  print(paste("Processed:", outputs[f_idx]))
}

# Extract relevant parts from column names for cleaner labels
# Assuming the pattern "- something _" as in your original script
p <- str_match(colnames(one), "-\\s*(.*?)\\s*_")
# Replace NAs that might result from columns like 'Gene'
clean_colnames <- ifelse(!is.na(p[,2]), p[,2], colnames(one))
colnames(one) <- clean_colnames
names(one)[1] <- "Gene" # Ensure 'Gene' is always the first column name

# Write the combined table to the specified output file
write.table(one, opt$output_file, quote = FALSE, row.names = FALSE, sep = "\t")

print(paste("Successfully combined RNA data into:", opt$output_file))