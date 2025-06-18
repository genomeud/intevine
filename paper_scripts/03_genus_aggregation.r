#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance16sT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/0.001_raw_counts_edited.xlsx",
              help="Input directory [default= %default]", metavar="character"),
 make_option(c("-B", "--abundance16sT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/0.001_raw_counts_edited.xlsx",
              help="Input directory [default= %default]", metavar="character"),
   make_option(c("-C", "--abundanceITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/0.001_raw_counts_edited.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-D", "--abundanceITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/0.001_raw_counts_edited.xlsx",
              help="Input directory [default= %default]", metavar="character"),
make_option(c("-E", "--output16sT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/",
              help="Input directory [default= %default]", metavar="character"),
 make_option(c("-F", "--output16sT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/",
              help="Input directory [default= %default]", metavar="character"),
   make_option(c("-G", "--outputITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/",
              help="Input directory [default= %default]", metavar="character"),
      make_option(c("-H", "--outputITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV//T1/06_tables_unite/",
              help="Input directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundance16sT0)) {
  stop("WARNING: No abundance16sT0 specified with '-A' flag.")
} else {  cat ("abundance16sT0 is ", opt$abundance16sT0, "\n")
  a16sT0 <- opt$abundance16sT0  
  }

if (is.null(opt$abundance16sT1)) {
  stop("WARNING: No abundance16sT1 specified with '-B' flag.")
} else {  cat ("abundancephylum is ", opt$abundance16sT1, "\n")
  a16sT1 <- opt$abundance16sT1
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

if (is.null(opt$output16sT0)) {
  stop("WARNING: No output16sT0 specified with '-E' flag.")
} else {  cat ("output16sT0 is ", opt$output16sT0, "\n")
  o16sT0 <- opt$output16sT0  
  }

if (is.null(opt$output16sT1)) {
  stop("WARNING: No 1 specified with '-F' flag.")
} else {  cat ("output16sT1 is ", opt$output16sT1, "\n")
  o16sT1 <- opt$output16sT1  
  }

  if (is.null(opt$outputITST0)) {
  stop("WARNING: No outputITST0 directory specified with '-G' flag.")
} else {  cat ("outputITST0 dir is ", opt$outputITST0, "\n")
  oITST0 <- opt$outputITST0 
  }
  
  if (is.null(opt$outputITST1)) {
  stop("WARNING: No outputITST1 directory specified with '-H' flag.")
} else {  cat ("outputITST1 is ", opt$outputITST1, "\n")
  oITST1 <- opt$outputITST1
  }


aggregation<-function() 
{
  library(data.table)
	library(openxlsx)
	library(viridis)
	library(tidyr)
	library(tibble)
  library(dplyr)
  listinput <- c(a16sT0,a16sT1,aITST0,aITST1)
  listoutput <- c(o16sT0,o16sT1,oITST0,oITST1)
  
  for (mycycle in 1:length(listinput))
  {
  	  print( paste0("Processing ", listinput[mycycle]))
	  myinput <- read.xlsx(listinput[mycycle])
	  myinput$Taxonomy <- gsub("ASV[^;]*;", "", myinput$Taxonomy)
	  semicolons_count <- sapply(gregexpr(";", myinput$Taxonomy), function(x) sum(x > 0))
	  myinput <- myinput[semicolons_count >= 5, ]
	  myinput [] <- lapply(myinput , function(x) 
	  {
		if (is.character(x)) {
		x <- gsub(" ", "_", x)
		}
		x
		})
		# Remove incertae sedis from the column genus
	  myinput$seq  <- NULL
	  #Since not all the input have the same number of samples and each sample is numered, with this trick i can find the number of samples for the aggregation function
	  numeric_col_names <- grep("^[0-9]+$", names(myinput), value = TRUE)
	  numeric_values <- as.numeric(numeric_col_names)
	  samplesnumber <- max(numeric_values)
	 # Aggregate data by Genus and sum columns 1-samplesnumber
	  aggregated_data <- myinput %>%
	  group_by(Taxonomy, ) %>%
	  summarise(across(1:samplesnumber, sum, na.rm = TRUE), .groups = "drop")
		aggregated_data <- aggregated_data %>%
	  mutate(Sums = rowSums(select(., `1`:samplesnumber))) %>%
	  arrange(desc(Sums)) %>%
	  select(-Sums)
	  aggregated_data <- as.data.frame(aggregated_data)
	  dir.create(listoutput[mycycle])
	  if (listinput[mycycle] == aITST0 | listinput[mycycle] == aITST1)
		{
		aggregated_data$Taxonomy <- gsub("k__", "", aggregated_data$Taxonomy)
		aggregated_data$Taxonomy <- gsub("p__", "", aggregated_data$Taxonomy)
		aggregated_data$Taxonomy <- gsub("c__", "", aggregated_data$Taxonomy)
		aggregated_data$Taxonomy <- gsub("o__", "", aggregated_data$Taxonomy)
		aggregated_data$Taxonomy <- gsub("f__", "", aggregated_data$Taxonomy)
		aggregated_data$Taxonomy <- gsub("g__", "", aggregated_data$Taxonomy)
		}
	  write.table(aggregated_data,paste0(listoutput[mycycle],"genus_aggregated_raw.txt"), row.names = FALSE, quote = FALSE)
	  write.xlsx(aggregated_data,paste0(listoutput[mycycle],"genus_aggregated_raw.xlsx"), row.names = FALSE, quote = FALSE)
	  }
}
aggregation()
