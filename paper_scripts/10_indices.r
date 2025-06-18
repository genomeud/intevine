#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/genus_aggregated.xlsx",
              help="Input file for T0 abundance [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/genus_aggregated.xlsx",
              help="Input file for T1 abundance [default= %default]", metavar="character"),
  make_option(c("-C", "--abundanceITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/genus_aggregated.xlsx",
              help="Input file for ITST0 abundance [default= %default]", metavar="character"),
  make_option(c("-D", "--abundanceITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/genus_aggregated.xlsx",
              help="Input file for ITST1 abundance [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/projects/marroni/intevine/docs/map_grapevine_final.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-O", "--outputT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/07_plot_SILVA/alpha/", 
              help="Output directory for T0 [default= %default]", metavar="character"),
  make_option(c("-P", "--outputT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/07_plot_SILVA/alpha/", 
              help="Output directory for T1 [default= %default]", metavar="character"),
  make_option(c("-Q", "--outputITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/07_plot_unite/alpha/", 
              help="Output directory for ITST0 [default= %default]", metavar="character"),
  make_option(c("-R", "--outputITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/07_plot_unite/alpha/", 
              help="Output directory for ITST1 [default= %default]", metavar="character"),
  make_option(c("-S", "--outtableT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/", 
              help="Output table directory for T0 [default= %default]", metavar="character"),
  make_option(c("-T", "--outtableT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/", 
              help="Output table directory for T1 [default= %default]", metavar="character"),
  make_option(c("-U", "--outtableITST0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/", 
              help="Output table directory for ITST0 [default= %default]", metavar="character"),
  make_option(c("-V", "--outtableITST1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/", 
              help="Output table directory for ITST1 [default= %default]", metavar="character"),
  make_option(c("-W", "--conditionlistT0"), type="character", default="Soil,Ster_soil", 
              help="Condition list for T0 [default= %default]", metavar="character"),
  make_option(c("-X", "--conditionlistT1"), type="character", default="Soil,Ster_soil,Ster_root", 
              help="Condition list for T1 [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$abundanceT0)) {
  stop("WARNING: No abundanceT0 specified with '-A' flag.")
} else {  
  cat("abundanceT0 is ", opt$abundanceT0, "\n")
  aT0 <- opt$abundanceT0  
}

if (is.null(opt$abundanceT1)) {
  stop("WARNING: No abundanceT1 specified with '-B' flag.")
} else {  
  cat("abundanceT1 is ", opt$abundanceT1, "\n")
  aT1 <- opt$abundanceT1  
}

if (is.null(opt$abundanceITST0)) {
  stop("WARNING: No abundanceITST0 specified with '-C' flag.")
} else {  
  cat("abundanceITST0 is ", opt$abundanceITST0, "\n")
  aITST0 <- opt$abundanceITST0  
}

if (is.null(opt$abundanceITST1)) {
  stop("WARNING: No abundanceITST1 specified with '-D' flag.")
} else {  
  cat("abundanceITST1 is ", opt$abundanceITST1, "\n")
  aITST1 <- opt$abundanceITST1  
}

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  
  cat("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
}

if (is.null(opt$conditionlistT0)) {
  stop("WARNING: No conditionlistT0 specified with '-W' flag.")
} else {  
  cat("conditionlistT0 is ", opt$conditionlistT0, "\n")
  conT0 <- opt$conditionlistT0  
}

if (is.null(opt$conditionlistT1)) {
  stop("WARNING: No conditionlistT1 specified with '-X' flag.")
} else {  
  cat("conditionlistT1 is ", opt$conditionlistT1, "\n")
  conT1 <- opt$conditionlistT1  
}

if (is.null(opt$outputT0)) {
  stop("WARNING: No outputT0 specified with '-O' flag.")
} else {  
  cat("outputT0 dir is ", opt$outputT0, "\n")
  oT0 <- opt$outputT0  
}

if (is.null(opt$outputT1)) {
  stop("WARNING: No outputT1 specified with '-P' flag.")
} else {  
  cat("outputT1 dir is ", opt$outputT1, "\n")
  oT1 <- opt$outputT1  
}

if (is.null(opt$outputITST0)) {
  stop("WARNING: No outputITST0 specified with '-Q' flag.")
} else {  
  cat("outputITST0 dir is ", opt$outputITST0, "\n")
  oITST0 <- opt$outputITST0 
}

if (is.null(opt$outputITST1)) {
  stop("WARNING: No outputITST1 specified with '-R' flag.")
} else {  
  cat("outputITST1 dir is ", opt$outputITST1, "\n")
  oITST1 <- opt$outputITST1
}

if (is.null(opt$outtableT0)) {
  stop("WARNING: No outtableT0 specified with '-S' flag.")
} else {  
  cat("outtableT0 dir is ", opt$outtableT0, "\n")
  otableT0 <- opt$outtableT0  
}

if (is.null(opt$outtableT1)) {
  stop("WARNING: No outtableT1 specified with '-T' flag.")
} else {  
  cat("outtableT1 dir is ", opt$outtableT1, "\n")
  otableT1 <- opt$outtableT1  
}

if (is.null(opt$outtableITST0)) {
  stop("WARNING: No outtableITST0 specified with '-U' flag.")
} else {  
  cat("outtableITST0 dir is ", opt$outtableITST0, "\n")
  otableITST0 <- opt$outtableITST0  
}

if (is.null(opt$outtableITST1)) {
  stop("WARNING: No outtableITST1 specified with '-V' flag.")
} else {  
  cat("outtableITST1 dir is ", opt$outtableITST1, "\n")
  otableITST1 <- opt$outtableITST1 
}

indices<-function() 
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
 library(ggpubr)
	inputlist <- c(aT0,aT1,aITST0,aITST1)
  outputlist <- c(oT0,oT1,oITST0,oITST1)
  tablelist <- c(otableT0, otableT1,otableITST0,otableITST1)
  namelist <- c("genus")
  for(myinput in 1:length(inputlist))
  {
   otable <- tablelist[myinput]
  outdir <- outputlist[myinput]
 
 countdata<-read.xlsx(inputlist[myinput])
  row.names(countdata) <- countdata$Taxonomy
	countdata$Taxonomy <- NULL
  if ( inputlist[myinput] == aT0 |  inputlist[myinput] == aITST0 )
  {
  allcond<-unlist(strsplit(conT0,","))
  metadata<-read.xlsx(metafile, sheet = "T0")
  } else {
  allcond<-unlist(strsplit(conT1,","))
  metadata<-read.xlsx(metafile, sheet = "Campionate")
  }
	bigmeta<-metadata
  bigcount <- countdata
  for(bbb in 1:length(allcond))
	{
	
metadata<-bigmeta
  countdata <- bigcount
  condition<-allcond[bbb]
   print(inputlist[myinput])
   print(outputlist[myinput])
   print(tablelist[myinput])
  	metadata<-metadata[,c("Sample_name",condition)]
	names(metadata)[2]<-"condition"
	metadata<-metadata[metadata$Sample_name%in%names(countdata),]
	rownames(metadata)<-metadata$Sample_name
	countdata<-data.frame(t(countdata))
	shannon<-diversity(countdata, index = "shannon", MARGIN = 1, base = exp(1))
	simpson<-diversity(countdata, index = "simpson", MARGIN = 1, base = exp(1))
	richness <- (apply(round(countdata, digits = 3), 1 , function(row) sum(row  > 0)))
write.table(shannon, paste0(otable,namelist[1],"_shannon.txt"), quote=F)
write.table(simpson, paste0(otable,namelist[1],"_simpson.txt"), quote=F)
write.table(richness, paste0(otable,"_",namelist[1],"_richness.txt"), quote=F)
  #x axis Taxonomy
  if (condition == "Soil") {
  namex= "Soil"
} else if (condition == "Ster_root") {
  namex= "Root sterilization"
} else if (condition == "Ster_soil") {
  namex= "Soil sterilization"
}
shannon_frame<-data.frame(shannon)
shannon_frame <- tibble::rownames_to_column(shannon_frame, "Sample_name")
shannon_meta <- merge (shannon_frame, metadata, by="Sample_name")
	if (allcond[bbb] == "Ster_soil" ){
	# Specify your pairwise comparisons
	my_comparisons <- list(
  c("Yes", "No"))
  colorcode <- c("DarkGreen", "LightGreen")
  } else if (allcond[bbb] == "Ster_root"  ){
  	# Specify your pairwise comparisons
	my_comparisons <- list(
  c("Yes", "No"))
  colorcode <- c("DarkBlue", "LightBlue")
  } else if (allcond[bbb] == "Soil" && (inputlist[myinput] == aT1 | inputlist[myinput] == aITST1 )) {
	my_comparisons <- list(
  c("Peat", "Manure"),
  c("Peat", "Sand"),
  c("Manure", "Sand")
)
colorcode <- c("#f5d902","#f50202","#070808")
} else if (allcond[bbb] == "Soil" && ( inputlist[myinput] == aT0 | inputlist[myinput] == aITST0 )) {
	my_comparisons <- list(
  c("Peat", "Manure")
)
colorcode <- c("#070808", "#f50202")
}
# Custom function to format pairwise Wilcoxon test results
format_pairwise_wilcox <- function(data, response, group) {
  pairwise_results <- pairwise.wilcox.test(data[[response]], data[[group]], p.adjust.method = "BH")
  results_df <- as.data.frame(pairwise_results$p.value)
  results_df$Comparison <- rownames(results_df)
  results_long <- tidyr::pivot_longer(results_df, cols = -Comparison, names_to = "Comparison2", values_to = "p.value")
  results_long <- results_long[!is.na(results_long$p.value), ]
  return(results_long)
}
# Apply the custom function
wilcoxon <- format_pairwise_wilcox(shannon_meta, "shannon", "condition")
# Display formatted results
	write.xlsx(wilcoxon,paste0(otable,condition,"_",namelist[1],"_shannon_wilcoxon.xlsx"), row.names = F)
	#wilcoxon_filtered <- wilcoxon[!grepl("ns", wilcoxon$p.value),]
	#if (allcond[bbb] == "Ster_soil" | allcond[bbb] == "Ster_root"  ){
	#wilcoxon_filtered <- wilcoxon_filtered %>% mutate(y.position = c(3.9))
	#} else if (allcond[bbb] == "Soil") {
	#wilcoxon_filtered <- wilcoxon_filtered %>% mutate(y.position = c(4, 4.1))
	if (allcond[bbb] == "Soil")
	{
	shannon_meta$condition <- factor(shannon_meta$condition,levels = c('Sand','Peat','Manure'),ordered = TRUE)
	}
	shannon_graph <- ggplot(data=shannon_meta, aes(x=condition, y=shannon, color=condition)) + 
    labs(y="Shannon index", x=namex) +
    ggtitle("Shannon Diversity Index") +
    theme(axis.line = element_line(colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.margin = margin(100, 10, 100, 10),
          legend.title = element_blank(),
          legend.position = "none", 
          axis.text = element_text(size=10)) +
    scale_color_manual(values=colorcode) +
    stat_compare_means(comparisons = my_comparisons) +
    geom_boxplot(fill = NA)  # Set fill to NA to remove box fill 
	#if(condition =="LC_simpl_2018")
	#{
	#shannon_graph <- shannon_graph + stat_compare_means(comparisons = my_comparisons, label = "p.signif",hide.ns = TRUE)
	#wilcoxon <- compare_means(shannon ~ condition, comparisons = my_comparisons, method='wilcox.test', data = shannon_meta)
	#write.xlsx(wilcoxon,paste0(otable,"shannon_wilcoxon_soil.xlsx"), row.names = F)
	#}
	pdf(paste0(outdir,condition,"_",namelist[1],"_shannon_box_plot_no_kruskal.pdf"))
	print(shannon_graph)
	dev.off()

	shannon_graph <- ggplot(data=shannon_meta, aes(x=condition, y=shannon, color=condition)) + 
    labs(y="Shannon index", x=namex) +
    ggtitle("Shannon Diversity Index") +
    theme(axis.line = element_line(colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.margin = margin(100, 10, 100, 10),
          legend.title = element_blank(),
          legend.position = "none", 
          axis.text = element_text(size=10)) +
    scale_color_manual(values=colorcode) +
    stat_compare_means(comparisons = my_comparisons) +
	 stat_compare_means(method = "kruskal.test")+
    geom_boxplot(fill = NA)  # Set fill to NA to remove box fill 
	#if(condition =="LC_simpl_2018")
	#{
	#shannon_graph <- shannon_graph + stat_compare_means(comparisons = my_comparisons, label = "p.signif",hide.ns = TRUE)
	#wilcoxon <- compare_means(shannon ~ condition, comparisons = my_comparisons, method='wilcox.test', data = shannon_meta)
	#write.xlsx(wilcoxon,paste0(otable,"shannon_wilcoxon_soil.xlsx"), row.names = F)
	#}
	pdf(paste0(outdir,condition,"_",namelist[1],"_shannon_box_plot.pdf"))
	print(shannon_graph)
	dev.off()
	
simpson_frame<-data.frame(simpson)
simpson_frame <- tibble::rownames_to_column(simpson_frame, "Sample_name")
simpson_meta <- merge (simpson_frame, metadata, by="Sample_name")
	# Apply the custom function
wilcoxon <- format_pairwise_wilcox(simpson_meta, "simpson", "condition")
# Display formatted results
	write.xlsx(wilcoxon,paste0(otable,condition,"_",namelist[1],"_simpson_wilcoxon_soil.xlsx"), row.names = F)
	#wilcoxon_filtered <- wilcoxon[!grepl("ns", wilcoxon$p.signif),]
	#if (allcond[bbb] == "Ster_soil" | allcond[bbb] == "Ster_root"  ){
	#wilcoxon_filtered <- wilcoxon_filtered %>% mutate(y.position = c(0.975))
#} else if (allcond[bbb] == "Soil") {
	#wilcoxon_filtered <- wilcoxon_filtered %>% mutate(y.position = c(0.975))
#}
		if (allcond[bbb] == "Soil")
	{
	simpson_meta$condition <- factor(simpson_meta$condition,levels = c('Sand','Peat','Manure'),ordered = TRUE)
	}
	simpson_graph <- ggplot(data=simpson_meta, aes(x=condition, y=simpson, color=condition)) + 
    labs(y="Simpson index", x=namex) +
    ggtitle("Simpson Diversity Index") +
    theme(axis.line = element_line(colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.margin = margin(100, 10, 100, 10),
          legend.title = element_blank(),
          legend.position = "none", 
          axis.text = element_text(size=10)) +
    scale_color_manual(values=colorcode) +
    stat_compare_means(comparisons = my_comparisons) +
    geom_boxplot(fill = NA)  
	#if(condition =="LC_simpl_2018")
	#{
	#simpson_graph <- simpson_graph + stat_compare_means(comparisons = my_comparisons, label = "p.signif",hide.ns = TRUE)
	#wilcoxon <- compare_means(simpson ~ condition, comparisons = my_comparisons, method='wilcox.test', data = simpson_meta)
	#write.xlsx(wilcoxon,paste0(otable,"simpson_wilcoxon_soil.xlsx"), row.names = F)
	#}
	pdf(paste0(outdir,condition,"_",namelist[1],"_simpson_box_plot.pdf"))
	print(simpson_graph)
	dev.off()
richness_frame<-data.frame(richness)
richness_frame <- tibble::rownames_to_column(richness_frame, "Sample_name")
richness_meta <- merge (richness_frame, metadata, by="Sample_name")
	# Apply the custom function
wilcoxon <- format_pairwise_wilcox(richness_meta, "richness", "condition")
# Display formatted results
	write.xlsx(wilcoxon,paste0(otable,condition,"_",namelist[1],"_richness_wilcoxon_soil.xlsx"), row.names = F)
	#wilcoxon_filtered <- wilcoxon[!grepl("ns", wilcoxon$p.signif),]
	#if (allcond[bbb] == "Ster_soil" | allcond[bbb] == "Ster_root"  ){
	#wilcoxon_filtered <- wilcoxon_filtered %>% mutate(y.position = c(59))
	#} else if (allcond[bbb] == "Soil") {
	#wilcoxon_filtered <- wilcoxon_filtered %>% mutate(y.position = c(61, 62))
#}
		if (allcond[bbb] == "Soil")
	{
	richness_meta$condition <- factor(richness_meta$condition,levels = c('Sand','Peat','Manure'),ordered = TRUE)
	}
	richness_graph <- ggplot(data=richness_meta, aes(x=condition, y=richness, color=condition)) + 
    labs(y="Richness", x=namex) +
    ggtitle("Richness") +
    theme(axis.line = element_line(colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.margin = margin(100, 10, 100, 10),
          legend.title = element_blank(),
          legend.position = "none", 
          axis.text = element_text(size=10)) +
    scale_color_manual(values=colorcode) +
    stat_compare_means(comparisons = my_comparisons) +
    geom_boxplot(fill = NA)
	#if(condition =="LC_simpl_2018")
	#{
	#	stat.test <- richness_meta %>%
	#	group_by("condition") %>%
	#	wilcox_test(richness ~ condition) %>%
	#	adjust_pvalue() %>%
	#	add_significance() %>%
	#	filter(p.adj.signif != "ns")
	#	richness_graph <- richness_graph + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = 7700, step.increase = 0.1)
	#wilcoxon <- compare_means(richness ~ condition, comparisons = my_comparisons, method='wilcox.test', data = richness_meta)
		pdf(paste0(outdir,condition,"_",namelist[1],"_richness_box_plot_no_kruskal.pdf"))
	print(richness_graph)
	dev.off()
	
		richness_graph <- ggplot(data=richness_meta, aes(x=condition, y=richness, color=condition)) + 
    labs(y="Richness", x=namex) +
    ggtitle("Richness") +
    theme(axis.line = element_line(colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.margin = margin(100, 10, 100, 10),
          legend.title = element_blank(),
          legend.position = "none", 
          axis.text = element_text(size=10)) +
    scale_color_manual(values=colorcode) +
    stat_compare_means(comparisons = my_comparisons) +
	stat_compare_means(method = "kruskal.test")+
    geom_boxplot(fill = NA)
	#if(condition =="LC_simpl_2018")
	#{
	#	stat.test <- richness_meta %>%
	#	group_by("condition") %>%
	#	wilcox_test(richness ~ condition) %>%
	#	adjust_pvalue() %>%
	#	add_significance() %>%
	#	filter(p.adj.signif != "ns")
	#	richness_graph <- richness_graph + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = 7700, step.increase = 0.1)
	#wilcoxon <- compare_means(richness ~ condition, comparisons = my_comparisons, method='wilcox.test', data = richness_meta)
		pdf(paste0(outdir,condition,"_",namelist[1],"_richness_box_plot.pdf"))
	print(richness_graph)
	dev.off()
	###################################REGRESSION#############################################
	if ( inputlist[myinput] == aT1 | inputlist[myinput] == aITST1 )
  {
  print("regression")
	regressionmeta <- bigmeta
	regressionmeta$shannon <- shannon
	regressionmeta$simpson <- simpson
	regressionmeta$richness <- richness
	indexcycle <- c("shannon","simpson","richness")
	for (myindex in 1:length(indexcycle))
		{
		for (mymeta in 8:25)
			{
			names(regressionmeta) <- gsub("\\..*","",names(regressionmeta))
			response_var <- indexcycle[myindex]
			predictor_var <- names(regressionmeta)[mymeta]
			print(predictor_var)
			# Construct the formula dynamically
			formula <- as.formula(paste(response_var, "~", predictor_var))
			# Fit the model using the dynamically created formula
			model <- lm(formula, data = regressionmeta)
			# Get the summary of the model
			model_summary <- summary(model)
			# Extract the p-value for the F-statistic
			f_p_value <- pf(model_summary$fstatistic[1], 
							model_summary$fstatistic[2], 
							model_summary$fstatistic[3], 
							lower.tail = FALSE)
			# Extract the p-value for the predictor variable 'Ca'
			ca_p_value <- coef(model_summary)[predictor_var, "Pr(>|t|)"]
			significance_level <- 0.05
			if (f_p_value < significance_level && ca_p_value < significance_level) {
			dir.create(paste0(outdir,"regression/",response_var,"/"), recursive = TRUE)
			pdf(paste0(outdir,"regression/",response_var,"/",predictor_var,"_",namelist[1],"_linear_regression.pdf"),width = 10, height = 7)
			regressionmeta$Autoclave_Ster_root <- paste(regressionmeta$Ster_soil,regressionmeta$Ster_root, sep = "_")
			# Define the pch values based on the levels of 'Autoclave_Ster_root'
			pch.levels <- c("No_No" = 1, "No_Yes" = 19, "Yes_Yes" = 17, "Yes_No" = 2)
			soil_colors <- setNames(c("Black", "Red", "darkgoldenrod2"), levels(as.factor(regressionmeta$Soil)))
			# Plot with DOC on the x-axis and Simpson index on the y-axis
			plots <- plot(regressionmeta[,predictor_var], regressionmeta[,response_var],  
				 xlab = predictor_var, 
				 ylab = response_var, 
				 pch = pch.levels[regressionmeta$Autoclave_Ster_root], # Set pch according to the 'Autoclave_Ster_root' level
				 col = soil_colors[regressionmeta$Soil], # Color by 'Soil' type
				 cex = 3,
				 main = paste0("f value: ", f_p_value, " | predictor variable p-value: ", ca_p_value)
			)
						# Add a legend
						if (!is.factor(regressionmeta$Soil)) {
						  regressionmeta$Soil <- as.factor(regressionmeta$Soil)
						}
# Add the legend for Soil types outside the plot area
legend(
  "topright",  # Position of the legend
  inset = c(-0.3, 0),  # Adjust inset to move the legend outside
  legend = levels(regressionmeta$Soil),  # Labels for the legend
  col = soil_colors,  # Colors for the legend
  pch = 16,  # Default point character
  title = "Soil Type",
  xpd = TRUE  # Allow drawing outside the plot region
)

# Add another legend for the point characters outside the plot area
legend(
  "bottomright",  # Position of the legend
  inset = c(-0.3, 0),  # Adjust inset to move the legend outside
  legend = names(pch.levels),  # Labels for the legend
  pch = pch.levels,  # Point characters for the legend
  col = "black",  # Color for the point characters
  title = "Autoclave Ster Root",
  xpd = TRUE  # Allow drawing outside the plot region
)
			# Add the regression line
			abline(model, col = "blue", lwd = 2)
			print(plots)
			dev.off()
			}
		}
				for (mymeta in 27:61)
			{
			names(regressionmeta) <- gsub("\\..*","",names(regressionmeta))
			regressionmeta[!grepl('16', regressionmeta$Sample_name),]
			response_var <- indexcycle[myindex]
			predictor_var <- names(regressionmeta)[mymeta]
			print(predictor_var)
			# Construct the formula dynamically
			formula <- as.formula(paste(response_var, "~", predictor_var))
			# Fit the model using the dynamically created formula
			model <- lm(formula, data = regressionmeta)
			# Get the summary of the model
			model_summary <- summary(model)
			# Extract the p-value for the F-statistic
			f_p_value <- pf(model_summary$fstatistic[1], 
							model_summary$fstatistic[2], 
							model_summary$fstatistic[3], 
							lower.tail = FALSE)
			# Extract the p-value for the predictor variable 'Ca'
			ca_p_value <- coef(model_summary)[predictor_var, "Pr(>|t|)"]
			significance_level <- 0.05
			if (f_p_value < significance_level && ca_p_value < significance_level) {
			dir.create(paste0(outdir,"regression/",response_var,"/"), recursive = TRUE)
			pdf(paste0(outdir,"regression/",response_var,"/",predictor_var,"_",namelist[1],"_linear_regression.pdf"),width = 10, height = 7)
			regressionmeta$Autoclave_Ster_root <- paste(regressionmeta$Ster_soil,regressionmeta$Ster_root, sep = "_")
			# Define the pch values based on the levels of 'Autoclave_Ster_root'
			pch.levels <- c("No_No" = 1, "No_Yes" = 19, "Yes_Yes" = 17, "Yes_No" = 2)
			soil_colors <- setNames(c("Black", "Red", "darkgoldenrod2"), levels(as.factor(regressionmeta$Soil)))
			# Plot with DOC on the x-axis and Simpson index on the y-axis
			plots <- plot(regressionmeta[,predictor_var], regressionmeta[,response_var],  
				 xlab = predictor_var, 
				 ylab = response_var, 
				 pch = pch.levels[regressionmeta$Autoclave_Ster_root], # Set pch according to the 'Autoclave_Ster_root' level
				 col = soil_colors[regressionmeta$Soil], # Color by 'Soil' type
				 cex = 3,
				  main = paste0("f value: ", f_p_value, " | predictor variable p-value: ", ca_p_value)
			)
						# Add a legend
						if (!is.factor(regressionmeta$Soil)) {
						  regressionmeta$Soil <- as.factor(regressionmeta$Soil)
						}
# Add the legend for Soil types outside the plot area
legend(
  "topright",  # Position of the legend
  inset = c(-0.3, 0),  # Adjust inset to move the legend outside
  legend = levels(regressionmeta$Soil),  # Labels for the legend
  col = soil_colors,  # Colors for the legend
  pch = 16,  # Default point character
  title = "Soil Type",
  xpd = TRUE  # Allow drawing outside the plot region
)

# Add another legend for the point characters outside the plot area
legend(
  "bottomright",  # Position of the legend
  inset = c(-0.3, 0),  # Adjust inset to move the legend outside
  legend = names(pch.levels),  # Labels for the legend
  pch = pch.levels,  # Point characters for the legend
  col = "black",  # Color for the point characters
  title = "Autoclave Ster Root",
  xpd = TRUE  # Allow drawing outside the plot region
)
			# Add the regression line
			abline(model, col = "blue", lwd = 2)
			print(plots)
			dev.off()
			}
		}
	}
	} 
	}
	}
	}
indices()