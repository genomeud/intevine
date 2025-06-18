#Modified by Fabio Marroni on 2022/03/19

# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundanceT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/genus_aggregated.xlsx",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--abundanceT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/genus_aggregated.xlsx",
              help="Input directory [default= %default]", metavar="character"),
   make_option(c("-M", "--metafile"), type="character", default="/projects/marroni/intevine/docs/map_grapevine_final.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-O", "--outputT0"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/07_plot_SILVA/", 
              help="output dir [default= %default]", metavar="character"),
    make_option(c("-P", "--outputT1"), type="character", default="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/07_plot_SILVA/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-G", "--genescorr"), type="numeric", default=25, 
              help="Comma separate list of samples to be removed, e.g. because of low yield [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

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
  
if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-V' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

if (is.null(opt$genescorr)) {
  stop("WARNING: No genescorr specified with '-V' flag.")
} else {  cat ("genescorr is ", opt$genescorr, "\n")
  genescorr <- opt$genescorr  
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


stacked<-function() 
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
	
	dir.create(outdir,recursive = TRUE)
  inputlist <- c(aT0,aT1)
  outputlist <- c(oT0,oT1)
  for(myinput in 1:length(inputlist))
  {
  	  metadata<-read.xlsx(metafile, sheet = "Campionate")
	  countdata<-read.xlsx(inputlist[myinput])
    #stacked bar plot
    rownames(countdata) <- countdata$Taxonomy
    countdata$Taxonomy <- NULL
    countdata$Average = rowMeans(countdata)
    countdata <- countdata[order(countdata$Average, decreasing = TRUE) ,]
    countdata$Average=NULL
    countdata <- head(countdata,10)
    countdata$Taxonomy = rownames(countdata)
    # Use pivot_longer to reshape the data
full3 <- countdata %>%
  pivot_longer(
    names_to = "Sample",     # Correct argument is names_to
    values_to = "Percentage",
    cols = -Taxonomy        # Exclude the Taxonomy column from being pivoted
  )
  #Full taxonomy graph
     graph= ggplot(full3,aes(fill=Taxonomy, y=Percentage, factor(Sample, level = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "14", "13", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "15", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36","Name")
			))) +
             geom_bar(position="stack", stat="identity") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_brewer(palette="Spectral") +
            #coord_flip() +
            ggtitle("") +
            #scale_x_discrete(limits=colnames(countdata)) +
            xlab("Sample") + ylab("Relative abundance (%)") +
            theme(axis.text.x=element_text(angle=90,margin = margin(0.5, unit = "cm"),vjust =1)) 
    graph <- graph + guides(fill=guide_legend(title="Genus"))
	dir.create(outputlist[myinput], recursive = TRUE)
    ggsave(paste0(outputlist[myinput],"stacked_bar_plot_full_taxonomy.png"),plot=graph, width = 30, units = "cm", device="png")
	  #Genus graph
	  # Modify the Taxonomy column to keep only the part after the last semicolon
	   genus <- full3
	  genus$Taxonomy <- sub(".*;", "", full3$Taxonomy)
    graph= ggplot(genus,aes(fill=Taxonomy, y=Percentage, factor(Sample, level = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "14", "13", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "15", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36","Name")
			))) +
             geom_bar(position="stack", stat="identity") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_brewer(palette="Spectral") +
            #coord_flip() +
            ggtitle("") +
            #scale_x_discrete(limits=colnames(countdata)) +
            xlab("Sample") + ylab("Relative abundance (%)") +
            theme(axis.text.x=element_text(angle=90,margin = margin(0.5, unit = "cm"),vjust =1)) 
    graph <- graph + guides(fill=guide_legend(title="Genus"))
	dir.create(outputlist[myinput], recursive = TRUE)
    ggsave(paste0(outputlist[myinput],"stacked_bar_plot_genus.png"),plot=graph, width = 30, units = "cm", device="png")
		  #Phylum graph
		   phylum <- full3
	phylum$Taxonomy <- sapply(strsplit(full3$Taxonomy, ";"), function(x) x[3])
	graph= ggplot(phylum,aes(fill=Taxonomy, y=Percentage, factor(Sample, level = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "14", "13", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "15", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36","Name")
			))) +
             geom_bar(position="stack", stat="identity") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_brewer(palette="Spectral") +
            #coord_flip() +
            ggtitle("") +
            #scale_x_discrete(limits=colnames(countdata)) +
            xlab("Sample") + ylab("Relative abundance (%)") +
            theme(axis.text.x=element_text(angle=90,margin = margin(0.5, unit = "cm"),vjust =1)) 
    graph <- graph + guides(fill=guide_legend(title="Genus"))
	dir.create(outputlist[myinput], recursive = TRUE)
    ggsave(paste0(outputlist[myinput],"stacked_bar_plot_phylum.png"),plot=graph, width = 30, units = "cm", device="png")
    }
    }
    stacked()

