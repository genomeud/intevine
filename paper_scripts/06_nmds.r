# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
})

# Define options for script arguments
option_list <- list(
  make_option(c("-A", "--abundanceT0"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_16sT0.xlsx",
              help="Input directory for abundanceT0 [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_16sT1.xlsx",
              help="Input directory for abundanceT1 [default= %default]", metavar="character"),
  make_option(c("-C", "--abundanceITST0"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_ITST0.xlsx",
              help="Input directory for abundanceITST0 [default= %default]", metavar="character"),
  make_option(c("-D", "--abundanceITST1"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_ITST1.xlsx",
              help="Input directory for abundanceITST1 [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/home/massimo.guazzini/map_grapevine_final.xlsx", 
              help="Metafile [default= %default]", metavar="character"),
  make_option(c("-O", "--outputT0"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/16s/ASV/SILVA/T0/07_plot_SILVA/", 
              help="Output directory for T0 [default= %default]", metavar="character"),
  make_option(c("-P", "--outputT1"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/16s/ASV/SILVA/T1/09_full_ASV/07_plot_SILVA/", 
              help="Output directory for T1 [default= %default]", metavar="character"),
  make_option(c("-Q", "--outputITST0"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/ITS/ASV/T0/07_plot_unite/", 
              help="Output directory for ITST0 [default= %default]", metavar="character"),
  make_option(c("-R", "--outputITST1"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/ITS/ASV/T1/07_plot_unite/", 
              help="Output directory for ITST1 [default= %default]", metavar="character"),
  make_option(c("-S", "--outtableT0"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/", 
              help="Output table directory for T0 [default= %default]", metavar="character"),
  make_option(c("-T", "--outtableT1"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/16s/ASV/SILVA/T1/09_full_ASV/06_tables_silva/", 
              help="Output table directory for T1 [default= %default]", metavar="character"),
  make_option(c("-U", "--outtableITST0"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/ITS/ASV/T0/06_tables_unite/", 
              help="Output table directory for ITST0 [default= %default]", metavar="character"),
  make_option(c("-V", "--outtableITST1"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/data_required/output_paper/metagenomic/ITS/ASV/T1/06_tables_unite/", 
              help="Output table directory for ITST1 [default= %default]", metavar="character"),
  make_option(c("-G", "--genescorr"), type="numeric", default=25, 
              help="Genes correlation threshold [default= %default]", metavar="numeric")
)

# Parse the options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

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
  
if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-V' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
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

# Load additional required libraries
suppressPackageStartupMessages({
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
  library(gridExtra)
})

# Define the main analysis function
run_analysis <- function(abundance_file, output_dir, metadata_sheet, table_output_dir, name_label, metadata_file) {
  dir.create(output_dir, recursive = TRUE)
  dir.create(table_output_dir, recursive = TRUE)
  metadata <- read.xlsx(metadata_file, sheet = metadata_sheet)
  names(metadata) <- sub("\\(.*", "", names(metadata))
  countdata <- read.xlsx(abundance_file)
  rownames(countdata) <- countdata$Taxonomy
  countdata$Taxonomy <- NULL
  print(abundance_file)
  # Convert data to matrix 
  count_matrix <- as.matrix(as.data.frame(t(countdata)))
  
  # Running NMDS with Bray-Curtis distance

  nmds <- metaMDS(count_matrix, distance = "bray", k = 2)
  site.scrs <- as.data.frame(scores(nmds, display = "sites"))
  if (abundance_file %in% c(aT1, aITST1)) {
  site.scrs <- cbind(site.scrs, Soil = metadata$Soil, Autoclave = metadata$Ster_soil, Heat_root <- metadata$Ster_root)
  } else site.scrs <- cbind(site.scrs, Soil = metadata$Soil, Autoclave = metadata$Ster_soil)
  stress_value <- round(nmds$stress, 3)
# Define the target column abbreviations
target_columns <- c(
  "GDVI",
  "NDVI",
  "CVI",
  "NDWI"
)
bigcount <- countdata
bigmeta <- metadata
  if (abundance_file %in% c(aT1, aITST1)) {
    metacycles <- list("", "Ionomics_", "Soilchem_","Multispectral_")
    env_columns <- list(NULL, 8:17, c(18, 20, 24,25), target_columns)
    for (i in seq_along(metacycles)) {
	countdata <- bigcount
	bigmeta <- metadata
      metacycle <- metacycles[[i]]
	  metacycle <- sub("\\(.*", "", metacycle)
      env_data <- if (!is.null(env_columns[[i]])) metadata[, env_columns[[i]], drop = FALSE] else NULL
      print(metacycle)
      # Handle Multispectral imaging special case
      if (metacycle == "Multispectral_") {
        countdata <- countdata %>% select(-`13`)
        metadata <- metadata %>% filter(Sample_name != "13")
        env_data <- if (!is.null(env_columns[[i]])) metadata[, env_columns[[i]], drop = FALSE] else NULL
        # Re-run NMDS with the modified data
        count_matrix <- as.matrix(as.data.frame(t(countdata)))
        nmds <- metaMDS(count_matrix, distance = "bray", k = 2)
        site.scrs <- as.data.frame(scores(nmds, display = "sites"))
        site.scrs <- cbind(site.scrs, Soil = metadata$Soil, Autoclave = metadata$Ster_soil, Heat_root <- metadata$Ster_root)
        stress_value <- round(nmds$stress, 3)
      } else {
	   metadata <- read.xlsx(metadata_file, sheet = metadata_sheet)
	    names(metadata) <- sub("\\(.*", "", names(metadata))
	    countdata <- read.xlsx(abundance_file)
	    rownames(countdata) <- countdata$Taxonomy
	    countdata$Taxonomy <- NULL
	   count_matrix <- as.matrix(as.data.frame(t(countdata)))
	  	    # Running NMDS with Bray-Curtis distance
		  nmds <- metaMDS(count_matrix, distance = "bray", k = 2)
		  site.scrs <- as.data.frame(scores(nmds, display = "sites"))
		  if (abundance_file %in% c(aT1, aITST1)) {
		  site.scrs <- cbind(site.scrs, Soil = metadata$Soil, Autoclave = metadata$Ster_soil, Heat_root <- metadata$Ster_root)
		  } else site.scrs <- cbind(site.scrs, Soil = metadata$Soil, Autoclave = metadata$Ster_soil)
		  stress_value <- round(nmds$stress, 3)
	  }
	   if (!is.null(env_data)) {
	        metanmds <- metaMDS(env_data,distance = "bray", k = 2)
	  metasite.scrs <- as.data.frame(scores(metanmds, display = "sites"))
	  metasite.scrs <- cbind(metasite.scrs, Soil = metadata$Soil, Autoclave = metadata$Ster_soil, Heat_root <- metadata$Ster_root)
	        # Full plot
      pdf(file.path(output_dir, paste0(metacycle, "meta_nmds_bray_curtis.pdf")))
	  metasite.scrs$Autoclave_Ster_root <- paste(metasite.scrs$Autoclave, metasite.scrs$Heat_root, sep = "_")
      pch.levels <- c("No_No" = 1, "No_Yes" = 19, "Yes_Yes" = 17, "Yes_No" = 2)
      metanmds.plot <- ggplot(metasite.scrs, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(aes(colour = factor(Soil), shape = Autoclave_Ster_root, size = 4)) + 
        scale_colour_manual(values = c("Manure" = "black", "Sand" = "darkgoldenrod2", "Peat" = "red")) + 
		scale_shape_manual(values = pch.levels) + 
        coord_fixed() +
        theme_classic() +
        theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
        labs(colour = "Soil", shape = "Autoclave") + 
        theme(legend.position = "right", legend.text = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
        ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
        print(metanmds.plot)
      dev.off()
	  }
      if (!is.null(env_data)) {
        # Fitting environmental parameters to the NMDS ordination
        envfit_res <- envfit(nmds, env_data, permutations = 999)
	capture.output(envfit_res, 
               file = file.path(output_dir, paste0(metacycle, "_envfit.txt")))
			   
        # Extracting significant vectors
        sig_arrows <- envfit_res$vectors$arrows[envfit_res$vectors$pvals < 0.05, ]
        sig_pvals <- envfit_res$vectors$pvals[envfit_res$vectors$pvals < 0.05]
        # Convert to a data frame for easy plotting
        sig_vectors <- as.data.frame(sig_arrows)
        sig_vectors$pvals <- sig_pvals
        sig_vectors$label <- rownames(sig_arrows)
		  # Scaling down the vectors
		  sig_vectors$NMDS1 <- sig_vectors$NMDS1 * 0.5
		  sig_vectors$NMDS2 <- sig_vectors$NMDS2 * 0.5

      }

      # Full plot
      pdf(file.path(output_dir, paste0(metacycle, "nmds_bray_curtis.pdf")))
	  site.scrs$Autoclave_Ster_root <- paste(site.scrs$Autoclave, site.scrs$Heat_root, sep = "_")
      pch.levels <- c("No_No" = 1, "No_Yes" = 19, "Yes_Yes" = 17, "Yes_No" = 2)
      nmds.plot <- ggplot(site.scrs, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(aes(colour = factor(Soil), shape = Autoclave_Ster_root, size = 4)) + 
        scale_colour_manual(values = c("Manure" = "black", "Sand" = "darkgoldenrod2", "Peat" = "red")) + 
		scale_shape_manual(values = pch.levels) + 
        coord_fixed() +
        theme_classic() +
        theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
        labs(colour = "Soil", shape = "Autoclave") + 
        theme(legend.position = "right", legend.text = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
        ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
      if (!is.null(env_data) && nrow(sig_vectors) > 0) {
        nmds.plot <- nmds.plot + 
          geom_segment(data = sig_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                       arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
          geom_text(data = sig_vectors, aes(x = NMDS1, y = NMDS2, label = label), 
                    hjust = 1.2, vjust = 1.2, size = 3, colour = "black")
      }
      
      print(nmds.plot)
      dev.off()
    }
  }

  # Soil only plot
  pdf(file.path(output_dir, "soil_only_nmds_bray_curtis.pdf"))
  nmds.plot <- ggplot(site.scrs, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(aes(colour = factor(Soil)), size = 4) + 
    scale_colour_manual(values = c("Manure" = "black", "Sand" = "darkgoldenrod2", "Peat" = "red")) + 
    coord_fixed() +
    theme_classic() +
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
    labs(colour = "Soil") + 
    theme(legend.position = "right", legend.text = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
    ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
  print(nmds.plot)
  dev.off()
  
    # # full Soil 
  # pdf(file.path(output_dir, "full_only_nmds_bray_curtis.pdf"))
   # pch.levels <- c("No_No" = 1, "No_Yes" = 19, "Yes_Yes" = 17, "Yes_No" = 2)
   # site.scrs$Autoclave_Ster_root <- paste(site.scrs$Autoclave, site.scrs$Heat_root, sep = "_")
  # nmds.plot <- ggplot(site.scrs, aes(x = NMDS1, y = NMDS2)) + 
    # geom_point(aes(colour = factor(Soil), shape = Autoclave_Ster_root, size = 4)) + 
    # scale_colour_manual(values = c("Manure" = "black", "Sand" = "darkgoldenrod2", "Peat" = "red")) + 
    # coord_fixed() +
    # theme_classic() +
    # theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
    # labs(colour = "Soil") + 
    # theme(legend.position = "right", legend.text = element_text(size = 10)) +
    # theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
    # ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
  # print(nmds.plot)
  # dev.off()
  
      if (abundance_file %in% c(aT0, aITST0)) {
    # Soil only plot
  pdf(file.path(output_dir, "soil_with_all_nmds_bray_curtis.pdf"))
    count_matrix <- as.matrix(as.data.frame(t(countdata)))
	  nmds <- metaMDS(count_matrix, distance = "bray", k = 2)
      site.scrs.cycle <- as.data.frame(scores(nmds, display = "sites"))
      site.scrs.cycle <- cbind(site.scrs.cycle, Soil = metadata$Soil, Autoclave = metadata$Ster_soil)
      stress_value <- round(nmds$stress, 3)
      site.scrs.cycle$Autoclave_Ster_root <- paste(site.scrs.cycle$Autoclave)
      pch.levels <- c("No" = 19, "Yes" = 17)
      nmds.plot <- ggplot(site.scrs.cycle, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(aes(colour = factor(Soil), shape = Autoclave_Ster_root), size = 4) + 
        scale_colour_manual(values = c("Manure" = "black", "Peat" = "red")) + 
        scale_shape_manual(values = pch.levels) +  # Map shapes to pch levels
        coord_fixed() +
        theme_classic() +
        theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
        labs(colour = "Soil", shape = "Autoclave") +  # Add shape legend
        theme(legend.position = "right", legend.text = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), 
              axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
        ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
      print(nmds.plot)
      dev.off()
  }
  
  # Plot per condition
  if (abundance_file %in% c(aT1, aITST1)) {
    for (condicycle in unique(metadata$Soil)) {
      cyclesamples <- metadata %>% filter(Soil == condicycle)
      cyclesample_names <- as.character(cyclesamples$Sample_name)
      cycle_countdata <- countdata %>% select(all_of(as.character(cyclesample_names)))
      count_matrix <- as.matrix(as.data.frame(t(cycle_countdata)))
	  print(condicycle)
	  nmds <- metaMDS(count_matrix, distance = "bray", k = 2)
      site.scrs.cycle <- as.data.frame(scores(nmds, display = "sites"))
      site.scrs.cycle <- cbind(site.scrs.cycle, Soil = cyclesamples$Soil, Autoclave = cyclesamples$Ster_soil, Heat_root = cyclesamples$Ster_root)
      stress_value <- round(nmds$stress, 3)
      site.scrs.cycle$Autoclave_Ster_root <- paste(site.scrs.cycle$Autoclave, site.scrs.cycle$Heat_root, sep = "_")
      pch.levels <- c("No_No" = 1, "No_Yes" = 19, "Yes_Yes" = 17, "Yes_No" = 2)
      pdf(file.path(output_dir, paste0(condicycle, "_nmds_bray_curtis.pdf")))
      nmds.plot <- ggplot(site.scrs.cycle, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(aes(colour = factor(Soil), shape = Autoclave_Ster_root), size = 4) + 
        scale_colour_manual(values = c("Manure" = "black", "Sand" = "darkgoldenrod2", "Peat" = "red")) + 
        scale_shape_manual(values = pch.levels) +  # Map shapes to pch levels
        coord_fixed() +
        theme_classic() +
        theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
        labs(colour = "Soil", shape = "Autoclave_Ster_root") +  # Add shape legend
        theme(legend.position = "right", legend.text = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), 
              axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
        ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
      print(nmds.plot)
      dev.off()
    }
  }
    if (abundance_file %in% c(aT0, aITST0)) {
    for (condicycle in unique(metadata$Soil)) {
      cyclesamples <- metadata %>% filter(Soil == condicycle)
      cyclesample_names <- as.character(cyclesamples$Sample_name)
      cycle_countdata <- countdata %>% select(all_of(as.character(cyclesample_names)))
      count_matrix <- as.matrix(as.data.frame(t(cycle_countdata)))
	  print(condicycle)
	  nmds <- metaMDS(count_matrix, distance = "bray", k = 2)
      site.scrs.cycle <- as.data.frame(scores(nmds, display = "sites"))
      site.scrs.cycle <- cbind(site.scrs.cycle, Soil = cyclesamples$Soil, Autoclave = cyclesamples$Ster_soil)
      stress_value <- round(nmds$stress, 3)
      site.scrs.cycle$Autoclave_Ster_root <- paste(site.scrs.cycle$Autoclave)
      pch.levels <- c("No" = 1, "Yes" = 17)
      pdf(file.path(output_dir, paste0(condicycle, "_nmds_bray_curtis.pdf")))
      nmds.plot <- ggplot(site.scrs.cycle, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(aes(colour = factor(Soil), shape = Autoclave_Ster_root), size = 4) + 
        scale_colour_manual(values = c("Manure" = "black", "Peat" = "red")) + 
        scale_shape_manual(values = pch.levels) +  # Map shapes to pch levels
        coord_fixed() +
        theme_classic() +
        theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
        labs(colour = "Soil", shape = "Autoclave") +  # Add shape legend
        theme(legend.position = "right", legend.text = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), 
              axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
        ggtitle(paste("NMDS Plot (Stress =", stress_value, ")"))
      print(nmds.plot)
      dev.off()
    }
  }
  
  # Writing output
  write.xlsx(site.scrs, file.path(table_output_dir, paste0("nmds_scores_", name_label, ".xlsx")), row.names = TRUE)
 

  # Return NMDS object for further analysis if needed
  return(nmds)
}

# Run analyses for each dataset
run_analysis(aT0, oT0, "T0", otableT0, "T0", metafile)
run_analysis(aT1, oT1, "Campionate_edited", otableT1, "T1", metafile)
run_analysis(aITST0, oITST0, "T0", otableITST0, "ITST0", metafile)
run_analysis(aITST1, oITST1, "Campionate_edited", otableITST1, "ITST1", metafile)



