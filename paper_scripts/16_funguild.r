#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/0.001_filtered_percentage_counts_funguild.guilds.reviewed.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/projects/marroni/intevine/docs/map_grapevine_final.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-P", "--outputT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/09_full_ASV/07_plot_unite/001/", 
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
  library("RColorBrewer")
  library("ggplot2")
  library("ggrepel")
  library(tidyr)
  library(dplyr)
  library(janitor)
  library(vegan)


  inputlist <- c(aT1)
  outputlist <- c(oT1)
  remove_possible <- c("Yes", "No")
  namelist <- c("genus")

  for (myinput in 1:length(inputlist)) {
    outdir <- outputlist[myinput]
        countdata <- read.xlsx(inputlist[myinput])
        countdata <- countdata[!grepl("Possible", countdata$Confidence.Ranking), ]
       unique_guilds <- unique(countdata$Guild)
      unique_guilds <- gsub("\\|", "", unique_guilds)
      split_guilds <- unlist(strsplit(unique_guilds, "-"))
      unique_guilds <- unique(split_guilds)
      unique_guilds <- unique_guilds[unique_guilds != ""]

      metadata <- read.xlsx(metafile, sheet = "Campionate_edited")
      allcond <- unlist(strsplit(conT1, ","))
      		countdata$seq <- countdata$Taxonomy <- countdata$OTU <- countdata$Taxon <- countdata$'Taxon.Level' <- countdata$'Trophic.Mode' <- countdata$Growth.Morphology <- countdata$Trait <- countdata$Confidence.Ranking <- countdata$Notes <- countdata$'Citation/Source' <- NULL

	  bigmeta <- metadata
      bigcount <- countdata

      for (bbb in 1:length(allcond)) {
        soil <- allcond[bbb]
         countdata <- bigcount
		  metadata <- bigmeta
		  forcycle<- countdata
		  forcycle <- forcycle[0,]
		  forcycle$Guild <- NULL
        # Loop through unique guilds and create a plot for each
	# soilal settings based on 'allcond' and 'inputlist'
if (allcond[bbb] == "Ster_soil" | allcond[bbb] == "Ster_root") {
  # Specify your pairwise comparisons
  my_comparisons <- list(
    c("Yes", "No")
  )
    if (allcond[bbb] == "Ster_soil") 
  {
	FillTrattamento <- c("DarkGreen", "LightGreen")
	namex <- "Autoclave"
   } else if (allcond[bbb] == "Ster_root") {
   FillTrattamento <- c("DarkBlue", "LightBlue")}
   namex <- "Heat root"
  } else if (allcond[bbb] == "Soil" && (inputlist[myinput] == aT1)) {
  my_comparisons <- list(
    c("Peat", "Manure"),
    c("Peat", "Sand"),
    c("Manure", "Sand")
	)
	namex <- "Soil"
  FillTrattamento <- c("#f5d902", "#f50202", "#070808")
}
        for (myloop in 1:length(unique_guilds)) {
		  countdatac <- bigcount
		  metadatac <- bigmeta
          print(unique_guilds[myloop])
          countdatac <- countdatac[grepl(unique_guilds[myloop], countdatac$Guild), ]
          columns_numbers <- grep("[0-9]", names(countdatac), value = TRUE)
          columns_keep <- c(columns_numbers, "Guild")
          countdatac <- countdatac[, columns_keep]
          countdatac$Guild <- NULL
          total <- colSums(countdatac, na.rm = TRUE)
          new_row <- data.frame(t(total))
          names(new_row) <- sub("^X", "", names(new_row))
          rownames(new_row) <- unique_guilds[myloop]
          countdatac <- rbind(countdatac, new_row)
          metadatac <- bigmeta
          metadatac <- metadatac[, c("Sample_name", soil)]
          names(metadatac)[2] <- "soil"
          metadatac <- metadatac[metadatac$Sample_name %in% names(countdatac), ]
          rownames(metadatac) <- metadatac$Sample_name
          countdatac <- data.frame(t(countdatac))
          columns_numbers <- grep("[0-9]", names(countdatac), value = TRUE)
          countdatac <- countdatac[, !names(countdatac) %in% columns_numbers]
				forcycle[nrow(forcycle) + 1, ] <- countdatac
			rownames(forcycle)[nrow(forcycle)] <- unique_guilds[myloop]
		}
		metadata<-metadata[,c("Sample_name",soil)]
		names(metadata)[2]<-soil
		countmeta <- as.data.frame(t(forcycle))
		countmeta$Sample_name<- rownames(countmeta)
		countmeta <- merge( countmeta, metadata, by = "Sample_name", sort = FALSE) 
		if (allcond[bbb] == "Soil") {
countmeta[[soil]] <- factor(countmeta[[soil]], levels = c('Sand', 'Peat', 'Manure'), ordered = TRUE)
}
if (allcond[bbb] == "Ster_soil" | allcond[bbb] == "Ster_root") {
countmeta[[soil]] <- factor(countmeta[[soil]], levels = c('Yes', 'No'), ordered = TRUE)
}
browser()
# Output directory and file names
outfile_base <- paste0(outdir, soil, "_", namelist, "_funguild_boxplot")
# Plot without Kruskal
abundance_threshold <- 0.01
prevalence_threshold <- 0.25
prevalent_columns <- apply(countmeta[, 2:26], 2, function(x) mean(x >= abundance_threshold) >= prevalence_threshold)
countmeta_filtered <- countmeta[, c(TRUE, prevalent_columns)]
exclude_processes <- c("Animal Parasite","Sample_name","Animal Pathogen","Algal Parasite","Nematophagous")
# Remove columns from countmeta_filtered where names match those in the exclude list
countmeta_filtered <- countmeta_filtered[, !grepl(paste(exclude_processes, collapse = "|"), colnames(countmeta_filtered))]
# Create an empty dataframe to store the results
# Create an empty dataframe to store the results
new_countmeta <- countmeta_filtered[, 1, drop = FALSE]  # Keep the first column (assuming it's 'Row.names')
# Track columns already processed or combined
processed_columns <- c()
# Loop through the columns of the dataframe
for(i in 2:(ncol(countmeta_filtered) - 1)) {  # Loop through columns (assuming column 1 is not relevant for comparison)
  if(!(colnames(countmeta_filtered)[i] %in% processed_columns)) {  # Only process columns that have not been combined
    # Initialize a vector to hold the columns that are to be combined
    columns_to_combine <- i  # Start with the current column

    # Loop through all subsequent columns to check for identical values
    for(j in (i + 1):ncol(countmeta_filtered)) {
      # Check if the columns are identical across all rows and if the second column has not been processed
      if(all(countmeta_filtered[, i] == countmeta_filtered[, j]) && !(colnames(countmeta_filtered)[j] %in% processed_columns)) {
        # Add the matching column to the list of columns to combine
        columns_to_combine <- c(columns_to_combine, j)
      }
    }
     # If multiple columns were found to be identical, combine them
    if(length(columns_to_combine) > 1) {
      # Create a new column name by combining the names of all matching columns
      new_col_name <- paste(colnames(countmeta_filtered)[columns_to_combine], collapse = ", ")
      # Assign the values from the first matching column to the new combined column
      new_countmeta[[new_col_name]] <- countmeta_filtered[, columns_to_combine[1]]
       # Mark all the columns as processed
      processed_columns <- c(processed_columns, colnames(countmeta_filtered)[columns_to_combine])
    } else {
      # If only one column was found, just add it to the new dataframe
      new_countmeta[[colnames(countmeta_filtered)[i]]] <- countmeta_filtered[, i]
    }
  }
}
new_countmeta[[soil]] <- countmeta_filtered[[soil]]
pdf(paste0(outdir, "/", soil, "_no_kruskal.pdf"))
# Maximum width for title
max_title_width <- 70

for (i in 2:(ncol(new_countmeta)-1)) {
  colname_i <- colnames(new_countmeta)[i]
  
  # Replace underscores with spaces for better readability
  colname_i_title <- gsub("_", " ", colname_i)
  
  # Use strwrap to split long titles into multiple lines
  colname_i_title <- paste(strwrap(colname_i_title, width = max_title_width), collapse = "\n")

  print(colname_i)
  
  max_y <- max(new_countmeta[[colname_i]], na.rm = TRUE)
  
  graphi2 <- ggboxplot(new_countmeta, x = allcond[bbb], y = colname_i,
                       color = allcond[bbb], palette = FillTrattamento) +
             ggtitle(colname_i_title) +
             ylab("Percentage count") +
			 xlab(namex) +
             stat_compare_means(comparisons = my_comparisons) +
             theme(legend.position = "none")
  
  print(graphi2)
}
dev.off()
# Plot with Kruskal
pdf(paste0(outfile_base, ".pdf"))
print(ncol(new_countmeta))

for (i in 2:(ncol(new_countmeta)-1)) {
  colname_i <- colnames(new_countmeta)[i]
  
  # Replace underscores with spaces for better readability
  colname_i_title <- gsub("_", " ", colname_i)
  
  # Use strwrap to split long titles into multiple lines
  colname_i_title <- paste(strwrap(colname_i_title, width = max_title_width), collapse = "\n")

  print(colname_i)
  
  max_y <- max(new_countmeta[[colname_i]], na.rm = TRUE)
  
  graphi2 <- ggboxplot(new_countmeta, x = allcond[bbb], y = colname_i,
                       color = allcond[bbb], palette = FillTrattamento) +
             ggtitle(colname_i_title) +
             ylab("Percentage count") +
			 xlab(namex) +
             stat_compare_means(comparisons = my_comparisons) +
             stat_compare_means(label.y = max_y) +
             theme(legend.position = "none")
  
  print(graphi2)
}
dev.off()
          # # Set up comparisons and colors based on condition
          # if (soil == "Ster_soil" | soil == "Ster_root") {
            # my_comparisons <- list(c("Yes", "No"))
            # FillTrattamento <- if (soil == "Ster_soil") c("DarkGreen", "LightGreen") else c("DarkBlue", "LightBlue")
          # } else if (soil == "Soil") {
            # my_comparisons <- list(c("Peat", "Manure"), c("Peat", "Sand"), c("Manure", "Sand"))
            # FillTrattamento <- c("#f5d902", "#f50202", "#070808")
          # }
          # # Merge data frames
          # countmeta <- merge(countdata, metadata, by = "row.names")
          # if (soil == "Soil") {
            # countmeta$soil <- factor(countmeta$soil, levels = c('Sand', 'Peat', 'Manure'), ordered = TRUE)
          # }
		  # browser()
          # # Plot with Kruskal
          # max_y <- max(countmeta$x, na.rm = TRUE)
          # graphi <- ggboxplot(countmeta, x = "soil", y = "x",
                              # color = "soil", palette = FillTrattamento) +
                    # ggtitle(unique_guilds[myloop]) +
                    # ylab("Total Percentage count") +
                    # stat_compare_means(comparisons = my_comparisons) +
                    # stat_compare_means(label.y = max_y) +
                    # theme(legend.position = "none")
          # print(graphi)
        }

      }
    }
  indices()