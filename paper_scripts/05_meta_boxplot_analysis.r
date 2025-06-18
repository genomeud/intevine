# Modified by Fabio Marroni on 2022/03/19

suppressPackageStartupMessages({
  library(optparse)
})

# Define command-line options
option_list = list(
  make_option(c("-M", "--metafile"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/Guazzini_Omics_analysis_grapevine_V03_Supplementary_file_1_Metadata.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/meta/new_pipeline/T1/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-C", "--conditions"), type="character", default=c("Soil", "Autoclave", "Heat root"), 
              help="conditions for analysis [default= %default]", metavar="character")    
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check for required inputs
if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  
  cat("Metafile is:", opt$metafile, "\n")
  metafile <- opt$metafile  
}

if (is.null(opt$out)) {
  stop("WARNING: No output directory specified with '-O' flag.")
} else {  
  cat("Output directory is:", opt$out, "\n")
  outdir <- opt$out  
}

if (is.null(opt$conditions)) {
  stop("WARNING: No conditions specified with '-C' flag.")
} else {  
  cat("Conditions are:", opt$conditions, "\n")
  cond <- opt$conditions  
}

correl <- function() {
  library(corrplot)
  library(data.table)
  library(openxlsx)
  library(RColorBrewer)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(rstatix)  # For post hoc tests

  # Create output directory if it does not exist
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  for (condcycle in 1:length(cond)) {
  print("test")
    metadata <- read.xlsx(metafile, sheet = "T1")
	colnames(metadata)[5:8] <- c("pH", "DOC", "ATP", "C_N")
	colnames(metadata)[13:22] <- sub("\\.\\(.*", "", colnames(metadata)[13:22])
    if (cond[condcycle] == "Soil") {  
      my_comparisons <- list(c("Manure", "Peat"), c("Peat", "Sand"), c("Sand", "Manure"))
      FillTrattamento <- c(Manure = "black", Peat = "red", Sand = "gold")
      metadata$Soil <- factor(metadata$Soil, levels = c('Sand', 'Peat', 'Manure'), ordered = TRUE)
    } else if (cond[condcycle] == "Autoclave") {
      my_comparisons <- list(c("Yes", "No"))
      FillTrattamento <- c(Yes = "DarkGreen", No = "LightGreen")
      metadata$Autoclave <- factor(metadata$Autoclave, levels = c('Yes', 'No'), ordered = TRUE)
    } else if (cond[condcycle] == "Heat root") {
      my_comparisons <- list(c("Yes", "No"))
      FillTrattamento <- c(Yes = "DarkBlue", No = "LightBlue")
      metadata$`Heat root` <- factor(metadata$`Heat root`, levels = c('Yes', 'No'), ordered = TRUE)
    }
    
    create_pdf_plots <- function(output_file, columns, include_kruskal) {
      # Ensure 'columns' is treated as a vector
      if (!is.vector(columns)) columns <- as.vector(columns)
      
      plots <- list()
      for (i in seq_along(columns)) {
	  print(colnames(metadata)[columns[i]])
        max_y <- max(metadata[[colnames(metadata)[columns[i]]]], na.rm = TRUE)
        plot <- ggboxplot(metadata, x = cond[condcycle], y = colnames(metadata)[columns[i]],
                          color = cond[condcycle], palette = FillTrattamento) +
                stat_compare_means(comparisons = my_comparisons) +
                theme(legend.position = "none",
                      axis.title = element_text(size = 16),
                      axis.text = element_text(size = 14))
        
        # Modify the y-axis label based on conditions:
        # If the column is "ATP", use its specific label.
        # Otherwise, if the plot is for "ionomics", append " (ppm)".
        if (colnames(metadata)[columns[i]] == "ATP") {
          plot <- plot + labs(y = "ATP (nmol/g soil)")
        } else if (grepl("ionomics", output_file)) {
          plot <- plot + labs(y = paste0(colnames(metadata)[columns[i]], " (ppm)"))
        }
        
        if (include_kruskal) {
          plot <- plot + stat_compare_means(label.y = max_y)
        }
        
        plots[[i]] <- plot
      }
      
      # Open a PDF device so that each plot is printed on its own page
      pdf(output_file, width = 10, height = 8)
      for (p in plots) {
        print(p)
      }
      dev.off()
    }
    
    create_pdf_plots(paste0(outdir, cond[condcycle], "_ionomics_boxplot_no_kruskal.pdf"), c(13:22), FALSE)
    create_pdf_plots(paste0(outdir, cond[condcycle], "_ionomics_boxplot.pdf"), c(13:22), TRUE)
    create_pdf_plots(paste0(outdir, cond[condcycle], "_soil_chemistry_boxplot_no_kruskal.pdf"), c(5:8), FALSE)
    create_pdf_plots(paste0(outdir, cond[condcycle], "_soil_chemistry_boxplot.pdf"), c(5:8), TRUE)
    create_pdf_plots(paste0(outdir, cond[condcycle], "_multispectral_boxplot_no_kruskal.pdf"), c(23:26), FALSE)
    create_pdf_plots(paste0(outdir, cond[condcycle], "_multispectral_boxplot.pdf"), c(23:26), TRUE)
    
    # Uncomment and modify below if needed for genetic extraction boxplots
    # create_pdf_plots(paste0(outdir, cond[condcycle], "_genetic_extraction_boxplot_no_kruskal.pdf"), c(23,24,25), FALSE)
    # create_pdf_plots(paste0(outdir, cond[condcycle], "_genetic_extraction_boxplot.pdf"), c(23,24,25), TRUE)
    
    # Additional plots (e.g., FDR-corrected) can be added here if needed.
  }
}

# Run the function
correl()
