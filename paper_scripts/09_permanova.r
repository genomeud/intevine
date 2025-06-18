	#extract samples from the ppm table, produced via DESeq2.r, by selecting a variable and a specific value. 
	#--abundance : ppm table of sample collected in "italy", produced via DESeq2.r
	#--metafile : Metadata of italian sample
	#--condition : "LC_simpl_2018", terrain type
	#--token : "Woodland", terrain that i want to select
	suppressPackageStartupMessages({
	  library(optparse)
	})

	option_list = list(
	  make_option(c("-A", "--abundanceT0"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_16sT0.xlsx",
				  help="Input directory [default= %default]", metavar="character"),
	  make_option(c("-B", "--abundanceT1"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_16sT1.xlsx",
				  help="Input directory [default= %default]", metavar="character"),
		make_option(c("-A", "--abundanceITST0"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_ITST0.xlsx",
				  help="Input directory [default= %default]", metavar="character"),
	  make_option(c("-B", "--abundanceITST1"), type="character", default="/home/massimo.guazzini/0.001_average_filtered_percentage_counts_edited_ITST1.xlsx",
				  help="Input directory [default= %default]", metavar="character"),
	  make_option(c("-M", "--metafile"), type="character", default="/home/massimo.guazzini/map_grapevine_final.xlsx", 
				  help="List of genes close to Vandal elements [default= %default]", metavar="character"),
	  make_option(c("-O", "--out"), type="character", default="/home/massimo.guazzini/permanova_new_pipeline/", 
				  help="output file name [default= %default]", metavar="character")
	)
	opt_parser = OptionParser(option_list=option_list);
	opt = parse_args(opt_parser);

if (is.null(opt$abundanceT1)) {
  stop("WARNING: No abundance specified with '-A' flag.")
} else {  cat ("abundance is ", opt$abundanceT1, "\n")
  abundanceT1 <- opt$abundanceT1  
  }

if (is.null(opt$abundanceITST1)) {
  stop("WARNING: No abundance specified with '-I' flag.")
} else {  cat ("abundance is ", opt$abundanceITST1, "\n")
  abundanceITST1 <- opt$abundanceITST1
  }
  
if (is.null(opt$abundanceT0)) {
  stop("WARNING: No abundanceT0 specified with '-A' flag.")
} else {  cat ("abundanceT0 is ", opt$abundanceT0, "\n")
  abundanceT0 <- opt$abundanceT0  
  }

if (is.null(opt$abundanceITST0)) {
  stop("WARNING: No abundance specified with '-I' flag.")
} else {  cat ("abundance is ", opt$abundanceITST0, "\n")
  abundanceITST0 <- opt$abundanceITST0
  }

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

 
  if (is.null(opt$out)) {
  stop("WARNING: No token specified with '-O' flag.")
} else {  cat ("out is ", opt$out, "\n")
  out <- opt$out  
  }
 
  library(data.table)
	library(ape)
	library(openxlsx)
   library("stringr")
  library("dplyr")
  library(janitor)
  library(vegan)
  library(ggplot2)
permanovz<-function()
{
  #library(pairwiseAdonis)
    metadata<-read.xlsx(metafile, sheet = "Campionate_edited")
	metadata$Soil <- as.factor(metadata$Soil)
	metadata$Ster_soil <- as.factor(metadata$Ster_soil)
	metadata$Ster_root <- as.factor(metadata$Ster_root)
	# Generate all combinations for Soil, Ster_soil, and Ster_root
condition_combinations <- expand.grid(
  first = c("Soil", "Ster_soil", "Ster_root"),
  second = c("Soil", "Ster_soil", "Ster_root"),
  third = c("Soil", "Ster_soil", "Ster_root")
)

# Remove invalid combinations where variables are duplicated in the order
condition_combinations <- condition_combinations[
  condition_combinations$first != condition_combinations$second &
  condition_combinations$second != condition_combinations$third &
  condition_combinations$first != condition_combinations$third, ]
# Save all nested combinations
	# Generate all combinations for Soil, Ster_soil
condition_combinationsT0 <- expand.grid(
  first = c("Soil", "Ster_soil"),
  second = c("Soil", "Ster_soil")
)
# Remove invalid combinationsT0 where variables are duplicated in the order
condition_combinationsT0 <- condition_combinationsT0[
  condition_combinationsT0$first != condition_combinationsT0$second, ]
# Save all nested combinations
dir.create(paste0(out, "meta/combos/"), recursive = TRUE)
# with the grouping factors
dir.create(paste0(out,"meta/chem/"),recursive = TRUE)
permanova_results <- adonis2(metadata[, c(18, 20, 24,25)] ~ Soil + Ster_soil + Ster_root, 
                            data = metadata, 
                            method = "bray")
capture.output(permanova_results, file= paste0 (out,"meta/chem/permanova_meta__ph_DOC_ATP.txt"))
permanova_results_interactions <- adonis2(metadata[, c(18, 20, 24,25)] ~ Soil * Ster_soil * Ster_root, 
                            data = metadata, 
                            method = "bray")
capture.output(permanova_results_interactions, file = paste0 (out,"meta/chem/permanova_meta__ph_DOC_ATP_interactions.txt"))
# with the grouping factors
for (i in 1:nrow(condition_combinations)) {
  nested_formula <- paste0(condition_combinations$first[i], "/",
          condition_combinations$second[i], "/",
          condition_combinations$third[i])
  
nested_formula <- as.formula(paste0("metadata[, c(18, 20, 24, 25)]"," ~ ", condition_combinations$first[i], " / ", condition_combinations$second[i], " / ", condition_combinations$third[i]))
  result <- adonis2(nested_formula, 
                    data = metadata, 
                    method = "bray")
  dir.create(paste0(out, "meta/combos/chem/"))
  capture.output(result, file = paste0(
    out, "meta/combos/chem/nested_permanova_",
    condition_combinations$first[i], "_",
    condition_combinations$second[i], "_",
    condition_combinations$third[i], ".txt"
  ))
}

 #############################soilchemT0#############

 metadataT0<-read.xlsx(metafile, sheet = "T0_formatted")
 metadataT0$Soil <- as.factor(metadataT0$Soil)
 metadataT0$Ster_soil <- as.factor(metadataT0$Ster_soil)
 # with the grouping factors
 dir.create(paste0(out,"meta/chem/T0/"),recursive = TRUE)
 permanova_results <- adonis2(metadataT0[, c(5,7,11,12)] ~ Soil + Ster_soil, 
                            data = metadataT0, 
                            method = "bray")
capture.output(permanova_results, file= paste0 (out,"meta/chem/T0/permanova_meta__ph_DOC_ATP.txt"))
permanova_results_interactions <- adonis2(metadataT0[, c(5,7,11,12)] ~ Soil * Ster_soil, 
                        data = metadataT0, 
                           method = "bray")
 capture.output(permanova_results_interactions, file = paste0 (out,"meta/chem/T0/permanova_meta__ph_DOC_ATP_interactions.txt"))
 # with the grouping factors
for (i in 1:nrow(condition_combinationsT0)) {
  nested_formula <- paste0(condition_combinationsT0$first[i], "/",
          condition_combinationsT0$second[i]) 
  
nested_formula <- as.formula(paste0("metadataT0[, c(5,7,11,12)]"," ~ ", condition_combinationsT0$first[i], " / ", condition_combinationsT0$second[i]))
  result <- adonis2(nested_formula, 
                    data = metadataT0, 
                    method = "bray")
   dir.create(paste0(out, "meta/combos/chem/")) 
  capture.output(result, file = paste0(
    out, "meta/combos/chem/T0_nested_permanova_",
    condition_combinationsT0$first[i], "_",
    condition_combinationsT0$second[i], ".txt"
  ))
}
###########################ionomics###############àà
# with the grouping factors
dir.create(paste0(out,"meta/ionomics/"),recursive = TRUE)
permanova_results <- adonis2(metadata[c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Na", "P", "S", "Zn")] ~ Soil + Ster_soil + Ster_root, 
                            data = metadata, 
                            method = "bray")
capture.output(permanova_results, file=paste0 (out,"meta/permanova_meta_ionomics.txt"))
permanova_results_interactions <- adonis2(metadata[c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Na", "P", "S", "Zn")] ~ Soil * Ster_soil * Ster_root, 
                            data = metadata, 
                            method = "bray")
capture.output(permanova_results_interactions, file = paste0 (out,"meta/ionomics/permanova_meta_ionomics_interactions.txt"))

for (i in 1:nrow(condition_combinations)) {
  nested_formula <- paste0(condition_combinations$first[i], "/",
          condition_combinations$second[i], "/",
          condition_combinations$third[i])
  
nested_formula <- as.formula(paste0("metadata[, c(8:17)]"," ~ ", condition_combinations$first[i], " / ", condition_combinations$second[i], " / ", condition_combinations$third[i]))
  result <- adonis2(nested_formula, 
                    data = metadata, 
                    method = "bray")
  dir.create(paste0(out, "meta/combos/ion/"))
  capture.output(result, file = paste0(
    out, "meta/combos/ion/nested_permanova_",
    condition_combinations$first[i], "_",
    condition_combinations$second[i], "_",
    condition_combinations$third[i], ".txt"
  ))
}
###########################multispectral###############àà
# with the grouping factors
dir.create(paste0(out,"meta/multispectral/"),recursive = TRUE)
print("multispectral")
metadataremoved <- metadata[-13, ]

permanova_results <- adonis2(metadataremoved[c( "NDVI", "CVI", "NDWI", "GDVI")] ~ Soil + Ster_soil + Ster_root, 
                            data = metadataremoved, 
                            method = "bray")
capture.output(permanova_results, file=paste0 (out,"meta/multispectral/permanova_meta_multispectral.txt"))
permanova_results_interactions <- adonis2(metadataremoved[c( "NDVI", "CVI", "NDWI", "GDVI")] ~ Soil * Ster_soil * Ster_root, 
                            data = metadataremoved, 
                            method = "bray")
capture.output(permanova_results_interactions, file = paste0 (out,"meta/multispectral/permanova_meta_multispectral_interactions.txt"))
for (i in 1:nrow(condition_combinations)) {
  nested_formula <- paste0(condition_combinations$first[i], "/",
          condition_combinations$second[i], "/",
          condition_combinations$third[i])
  
nested_formula <- as.formula(paste0("metadataremoved[, c(41,46,47,51)]"," ~ ", condition_combinations$first[i], " / ", condition_combinations$second[i], " / ", condition_combinations$third[i]))
  result <- adonis2(nested_formula, 
                    data = metadataremoved, 
                    method = "bray")
  dir.create(paste0(out, "meta/combos/multi/"))
  capture.output(result, file = paste0(
    out, "meta/combos/multi/nested_permanova_",
    condition_combinations$first[i], "_",
    condition_combinations$second[i], "_",
    condition_combinations$third[i], ".txt"
  ))
}

####################16s#################
dir.create(paste0 (out,"rna/"),recursive = TRUE)
print("16s")
browser()
  	countdata<-read.xlsx(abundanceT1)
	rownames(countdata) <- countdata$Taxonomy
	countdata$Taxonomy <- NULL
	flipped<- as.data.frame(t(countdata))
	diss_matrix <- vegdist(flipped, method = 'bray')
	permanova_results <- adonis2(diss_matrix ~ Soil + Ster_soil + Ster_root, data = metadata)
capture.output(permanova_results, file = paste0 (out,"rna/permanova_rna.txt"))
	permanova_results_interactions <- adonis2(diss_matrix ~ Soil * Ster_soil * Ster_root, data = metadata)
capture.output(permanova_results_interactions, file = paste0 (out,"rna/permanova_rna_interactions.txt"))
	permanova_results_nested <- adonis2(diss_matrix ~ Soil / Ster_soil / Ster_root, data = metadata)
capture.output(permanova_results_nested, file = paste0 (out,"rna/permanova_rna_nested.txt"))
for (i in 1:nrow(condition_combinations)) {
  nested_formula <- paste0(condition_combinations$first[i], "/",
          condition_combinations$second[i], "/",
          condition_combinations$third[i])
nested_formula <- as.formula(paste0("diss_matrix"," ~ ", condition_combinations$first[i], " / ", condition_combinations$second[i], " / ", condition_combinations$third[i]))
  result <- adonis2(nested_formula, 
                    data = metadata, 
                    method = "bray")
  dir.create(paste0(out, "rna/combos/"))
  capture.output(result, file = paste0(
    out, "16s/combos/nested_permanova_",
    condition_combinations$first[i], "_",
    condition_combinations$second[i], "_",
    condition_combinations$third[i], ".txt"
  ))
}
#16s impact of doc e ph#
  	# countdata<-read.xlsx(abundance)
	# rownames(countdata) <- countdata$Taxonomy
	# countdata$Taxonomy <- NULL
	# flipped<- as.data.frame(t(countdata))
####################ITS##########################àà
print("ITS")
dir.create(paste0 (out,"ITS/"),recursive = TRUE)
  	countdataITS<-read.xlsx(abundanceITST1)
	rownames(countdataITS) <- countdataITS$Taxonomy
	countdataITS$Taxonomy <- NULL
	flippedITS<- as.data.frame(t(countdataITS))
	diss_matrix <- vegdist(flippedITS, method = 'bray')
	permanova_results <- adonis2(diss_matrix ~ Soil + Ster_soil + Ster_root, data = metadata)
capture.output(permanova_results, file = paste0 (out,"ITS/permanova_ITS.txt"))
	permanova_results_interactions <- adonis2(diss_matrix ~ Soil * Ster_soil * Ster_root, data = metadata)
capture.output(permanova_results_interactions, file = paste0 (out,"ITS/permanova_ITS_interactions.txt"))
for (i in 1:nrow(condition_combinations)) {
  nested_formula <- paste0(condition_combinations$first[i], "/",
          condition_combinations$second[i], "/",
          condition_combinations$third[i])
nested_formula <- as.formula(paste0("diss_matrix"," ~ ", condition_combinations$first[i], " / ", condition_combinations$second[i], " / ", condition_combinations$third[i]))
  result <- adonis2(nested_formula, 
                    data = metadata, 
                    method = "bray")
  dir.create(paste0(out, "ITS/combos/"))
  capture.output(result, file = paste0(
    out, "ITS/combos/nested_permanova_",
    condition_combinations$first[i], "_",
    condition_combinations$second[i], "_",
    condition_combinations$third[i], ".txt"
  ))
}
####################FAPROTAX##########################àà
  	# countdatafapro<-read.xlsx(faprotax)
	# rownames(countdatafapro) <- countdatafapro$Taxonomy
	# countdatafapro$group <- NULL
	# flippedfaprotax<- as.data.frame(t(countdatafapro))
	# diss_matrix <- vegdist(flippedfaprotax, method = 'bray')
	# permanova_results <- adonis2(diss_matrix ~ Soil + Ster_soil + Ster_root, data = metadata)
# capture.output(permanova_results, file = paste0 (out,"faprotax/permanova_faprotax.txt"))
	# permanova_results <- adonis2(diss_matrix ~ Soil * Ster_soil * Ster_root, data = metadata)
# capture.output(permanova_results, file = paste0 (out,"faprotax/permanova_faprotax_interactions.txt"))
# for (i in 1:length(fulllist))
	# {
	# mycycle <- fulllist[i]
	# print(mycycle)
	# pair <- pairwise.adonis2(diss_matrix,metadata[[mycycle]], perm = 999, p.adjust.m='holm')
	# capture.output(pair, file = paste0 (out,"faprotax/",mycycle, "_pairwise_permanova_faprotax.txt"))
  # }
	# permanova_results <- adonis2(diss_matrix ~ Soil / Ster_soil / Ster_root, data = metadata)
# capture.output(permanova_results, file = paste0 (out,"faprotax/permanova_faprotax_nested.txt"))
####################T016s#################
print("T016s")

  	    metadataT0<-read.xlsx(metafile, sheet = "T0")
		dir.create(paste0 (out,"16s/T0/"),recursive = TRUE)
  	countdata<-read.xlsx(abundanceT0)
	rownames(countdata) <- countdata$Taxonomy
	countdata$Taxonomy <- NULL
	flipped<- as.data.frame(t(countdata))
	metadataT0$Soil <- as.factor(metadataT0$Soil)
	metadataT0$Ster_soil <- as.factor(metadataT0$Ster_soil)
	metadataT0$Soil_Autoclave <- paste0(metadataT0$Soil, "_", metadataT0$Ster_soil)
	flipped<- as.data.frame(t(countdata))
	diss_matrix <- vegdist(flipped, method = 'bray')
	permanova_results <- adonis2(diss_matrix ~ Soil + Ster_soil , data = metadataT0)
	dir.create(paste0 (out,"16s/T0/"))
capture.output(permanova_results, file = paste0 (out,"16s/T0/permanova_16s.txt"))
	permanova_results_interactions <- adonis2(diss_matrix ~ Soil * Ster_soil , data = metadataT0)
capture.output(permanova_results_interactions, file = paste0 (out,"16s/T0/permanova_16s_interactions.txt"))
print("nested")
for (i in 1:nrow(condition_combinationsT0)) {
  nested_formula <- paste0(condition_combinationsT0$first[i], "/",
          condition_combinationsT0$second[i]) 
		nested_formula <- as.formula(paste0("diss_matrix"," ~ ", condition_combinationsT0$first[i], " / ", condition_combinationsT0$second[i]))
  result <- adonis2(nested_formula, 
                    data = metadataT0, 
                    method = "bray")
  dir.create(paste0(out, "16s/T0/combos/"))
  capture.output(result, file = paste0(
    out, "16s/T0/combos/nested_permanova_",
    condition_combinationsT0$first[i], "_",
    condition_combinationsT0$second[i], ".txt"
  ))
}
#16s impact of doc e ph#
  	# countdata<-read.xlsx(abundance)
	# rownames(countdata) <- countdata$Taxonomy
	# countdata$Taxonomy <- NULL
	# flipped<- as.data.frame(t(countdata))
####################T0ITS##########################àà
print("T0ITS")
  	countdataITS<-read.xlsx(abundanceITST0)
	dir.create(paste0 (out,"ITS/T0/"),recursive = TRUE)
	rownames(countdataITS) <- countdataITS$Taxonomy
	countdataITS$Taxonomy <- NULL
	flippedITS<- as.data.frame(t(countdataITS))
	diss_matrix <- vegdist(flippedITS, method = 'bray')
	dir.create(paste0 (out,"ITS/T0/"))
	permanova_results <- adonis2(diss_matrix ~ Soil + Ster_soil , data = metadataT0)
capture.output(permanova_results, file = paste0 (out,"ITS/T0/permanova_ITS.txt"))
	permanova_results_interactions <- adonis2(diss_matrix ~ Soil * Ster_soil , data = metadataT0)
capture.output(permanova_results_interactions, file = paste0 (out,"ITS/T0/permanova_ITS_interactions.txt"))
for (i in 1:nrow(condition_combinationsT0)) {
  nested_formula <- paste0(condition_combinationsT0$first[i], "/",
          condition_combinationsT0$second[i]) 
		nested_formula <- as.formula(paste0("diss_matrix"," ~ ", condition_combinationsT0$first[i], " / ", condition_combinationsT0$second[i]))
  result <- adonis2(nested_formula, 
                    data = metadataT0, 
                    method = "bray")
  dir.create(paste0(out, "ITS/T0/combos/"))
  capture.output(result, file = paste0(
    out, "ITS/T0/combos/nested_permanova_",
    condition_combinationsT0$first[i], "_",
    condition_combinationsT0$second[i], ".txt"
  ))
  }
}
permanovz()
