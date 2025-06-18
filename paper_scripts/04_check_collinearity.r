#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
   make_option(c("-M", "--metafile"), type="character", default="/home/massimo.guazzini/map_grapevine_final.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-O", "--output"), type="character", default="/home/massimo.guazzini/intevine/analysis/meta/new_pipeline/collinearity/", 
              help="output dir [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")
  
if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-V' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

  if (is.null(opt$output)) {
  stop("WARNING: No output specified with '-I' flag.")
} else {  cat ("output dir is ", opt$output, "\n")
  outdir <- opt$output  
  }
  
correlation<-function() 
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
	library(reshape2)
	library(vegan)
	library(Hmisc)
	dir.create(outdir, recursive = TRUE)
get_upper_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}

  get_lower_tri <- function(CorMat){
  CorMat[lower.tri(CorMat)]<- NA
  return(CorMat)
}

reorder <- function(CorMat){
  dd <- as.dist((1-CorMat)/2)
  hc <- hclust(dd)
  CorMar <- CorMat[hc$order, hc$order]
}

metadata <- read.xlsx(metafile, sheet = "Campionate_edited")
metacycle <- c("Ionomics","Soil_Chemistry","Multispectral_Imaging")
for (aaa in metacycle) {
  if (aaa == "Ionomics") {
    matrixs <- metadata[, 8:17]
  } else if (aaa == "Soil_Chemistry") {
    matrixs <- metadata[, 18:25]
  } else if (aaa == "Multispectral_Imaging") {
    matrixs <- metadata[-13, c(46,47,51,54)]
  }
  
CorMat <- rcorr(as.matrix(matrixs), type = "pearson")
CorMatcorr <- CorMat$r
CorMatpvalue <- CorMat$P
# Create a new workbook
wb <- createWorkbook()
# Add sheets to the workbook
addWorksheet(wb, "Correlation")
addWorksheet(wb, "Pvalue")
# Write data to the sheets
writeData(wb, "Correlation", CorMatcorr, rowNames = TRUE)
writeData(wb, "Pvalue", CorMatpvalue, rowNames = TRUE)
# Save the workbook to a file
output_file <- paste0(outdir, aaa ,"_correlation_matrix.xlsx")
saveWorkbook(wb, file = output_file, overwrite = TRUE)
CorMatcorr <- reorder(CorMatcorr)
CorMatpvalue <- reorder(CorMatpvalue)
# Initialize the result matrix with the original correlation values
corrmatcorpvalue <- CorMatcorr
corrmatcorpvalue <- round(CorMatcorr, 2)
# Iterate over columns and rows of the p-value matrix
for (i in colnames(CorMatpvalue)) {
  for (j in rownames(CorMatpvalue)) {
    # Determine the appropriate asterisk based on the p-value
    asterisk <- ifelse(is.na(CorMatpvalue[j, i]), "", 
                       ifelse(CorMatpvalue[j, i] < 0.001, "***",
                              ifelse(CorMatpvalue[j, i] < 0.01, "**",
                                     ifelse(CorMatpvalue[j, i] < 0.05, "*", ""))))
    # Append the asterisk to the correlation value
    corrmatcorpvalue[j, i] <- paste0(corrmatcorpvalue[j, i], asterisk)
  }
}
upper_tri <- get_upper_tri(CorMatcorr)
lower_tri <- get_lower_tri(corrmatcorpvalue)
meltNum <- melt(lower_tri, na.rm = T)
meltColor <- melt(upper_tri, na.rm = T)
lower_tri_pvalue <- get_lower_tri(CorMatpvalue)
meltNumpvalue <- melt(lower_tri_pvalue, na.rm = T)
text_size <- if (aaa == "Multispectral_Imaging") 1.5 else 4
xaxis <- if (aaa == "Multispectral_Imaging" |aaa ==  "Soil_Chemistry") 45 else 0 
hjustvalue <- if (aaa == "Multispectral_Imaging" |aaa ==  "Soil_Chemistry") -0.1 else 0 
vjustvalue <- if (aaa == "Multispectral_Imaging" |aaa ==  "Soil_Chemistry") 5 else 0 
ploti <- ggplot() +
  labs(x = NULL, y = NULL) +
  geom_tile(data = meltColor, 
            mapping = aes(Var2, Var1, 
                          fill = value)) +
  geom_text(data = meltNum,
            mapping = aes(Var2, Var1,
                          label = value),
            size = text_size) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "firebrick4",
                      limit = c(-1,1), name = "Pearson\nCorrelation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = xaxis, hjust = hjustvalue, vjust = vjustvalue)) +
  coord_fixed()
  pdf(paste0(outdir, aaa, "_correlation.pdf"))
 print(ploti)
 dev.off()
}
}
correlation()
