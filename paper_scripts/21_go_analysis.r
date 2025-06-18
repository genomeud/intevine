# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="/home/massimo.guazzini/transcriptomic/DESeq2/",
		help="Input file (ZHp, HapFLK or other similar format) [default= %default]", metavar="character"), 
  make_option(c("-K", "--keggfile"), type="character", default="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/named_transcript_GO_KEGG_path.txt",
		help="Formatted KEGG+GO file with transcript names [default= %default]", metavar="character") 
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-G' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$keggfile)) {
  stop("WARNING: No keggfile specified with '-K' flag.")
} else {  cat ("keggfile ", opt$keggfile, "\n")
  keggfile <- opt$keggfile  
  }

kegg.enrich<-function(infile,keggfile,outfile,minpos=3,p.value=0.05)
{
suppressPackageStartupMessages(library(dplyr))
library(data.table)
library(openxlsx) 
library(GO.db)
library(UpSetR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(stringr)
 final_dataframe <- data.frame() 
 total_dataframe <- data.frame(class = character(0))
conditionz <- list.files(infile, pattern = "^S")
conditionz <- conditionz[!grepl("auto", conditionz)]
for (ddd in conditionz) {
setwd(paste0(infile,ddd,"/padj"))
comparisonz <- list.files(".")
comparisonz <- comparisonz[grep("vs", comparisonz)]
comparisonz <- comparisonz[!grepl("KEGG", comparisonz)]
comparisonz <- comparisonz[!grepl("GO", comparisonz)]
for (ggg in comparisonz) {
indata <- fread(ggg, data.table = F)
cyclex <- c("positive","negative","all")
bigin <- indata
for (i in cyclex)
{
indata <- bigin
if(nrow(indata)> 2) {
bdata<-fread(keggfile,data.table=F)
alldata <- merge(indata,bdata, by.x = "class", by.y = "transcript", all = FALSE, sort = FALSE)
alldatago <- merge(indata,bdata, by.x = "class", by.y = "transcript", all = FALSE, sort = FALSE)
if (i == "negative")
	{
	alldatago <- alldatago[grepl("-", alldatago$log2FoldChange), ]
	indata <- indata[grepl("-", indata$log2FoldChange), ]
	}else if (i =="positive") {
	alldatago <- alldatago[!grepl("-", alldatago$log2FoldChange), ]	
	indata <- indata[!grepl("-", indata$log2FoldChange), ]
	}
name <- gsub(".txt","",ggg)
write.xlsx(alldata,paste0("KEGG_",name,".xlsx"),row.names=F)
sigdata <- alldata[alldata$padj < 0.05,] 
sigdatago <- alldatago[alldatago$padj < 0.05,] 
#################################INIZIO ENRICHMENT: BDATA è IL DATABASE CON TUTTO, SIGDATA I "MODULI" NEL TUO CASO
 TERM2GENE <- bdata %>%
  tidyr::separate_rows(GO, sep = ";") %>%
  dplyr::select(GO, uniprot) %>%
  dplyr::distinct() %>%
  dplyr::rename(Term = GO, Gene = uniprot) 
# Map GO IDs to descriptions
TERM2NAME <- data.frame(
  Term = keys(GO.db, keytype = "GOID"),
  Description = AnnotationDbi::select(GO.db, keys = keys(GO.db, keytype = "GOID"), 
                                      columns = "TERM", keytype = "GOID")[, "TERM"]
)
# Merge TERM2GENE with TERM2NAME
TERM2GENE <- TERM2GENE %>%
  left_join(TERM2NAME, by = "Term")
gene_list <- sigdatago$uniprot
  go_enrichment <- enricher(
  gene = gene_list,
  TERM2GENE = TERM2GENE,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
# Convert go_enrichment to a data frame
go_results <- as.data.frame(go_enrichment)
go_results$Description <- NULL
formatted_go <- go_results %>%
  left_join(TERM2GENE %>% dplyr::select(Term, Description) %>% distinct(), by = c("ID" = "Term")) %>%
  mutate(
    # Split GeneRatio and BgRatio into numeric numerator/denominator
    GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatioDen = as.numeric(sub(".*/", "", GeneRatio)),
    GeneRatioVal = GeneRatioNum / GeneRatioDen,
    BgRatioNum = as.numeric(sub("/.*", "", BgRatio)),
    BgRatioDen = as.numeric(sub(".*/", "", BgRatio)),
    BgRatioVal = BgRatioNum / BgRatioDen
  ) %>%
  # Select and rename columns for clarity
  transmute(
    class = Description,                  # GO term description
    classratio = paste(GeneRatioNum, "/", GeneRatioDen, sep = ""),
    classbg = paste(BgRatioNum, "/", BgRatioDen, sep = ""),
    classp = pvalue,                      # Raw p-value
    classfdr = p.adjust                   # Adjusted p-value (FDR)
  )
  formatted_go <- na.omit(formatted_go)
  write.table(formatted_go,paste0(i,"_GO_", ggg),quote=F, row.names=F)
 # Select top GO terms based on FDR for bar plot and dot plot
top_go <- formatted_go %>% 
  arrange(classfdr)%>% 
  slice_head(n = 20)   # Select top 20 terms
# Ensure columns are properly formatted for calculations
top_go <- top_go %>%
  mutate(
    GeneRatioVal = as.numeric(sub("/.*", "", classratio)) / as.numeric(sub(".*/", "", classratio)),
    BgRatioVal = as.numeric(sub("/.*", "", classbg)) / as.numeric(sub(".*/", "", classbg))
  )
### 1. Bar Plot
bar_plot <- ggplot(top_go, aes(x = reorder(class, -classfdr), y = GeneRatioVal, fill = classfdr)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    x = "GO Terms",
    y = "Gene Ratio",
    fill = "FDR",
    title = "Enriched GO Terms (Bar Plot)"
  ) +
  theme_minimal()
# Print the bar plot
dir.create("plot")
nameplot <- gsub(".txt","",ggg)
pdf(paste0("plot/",i,"_GO_barplot_",nameplot,".pdf"))
print(bar_plot)
dev.off()
### 2. Dot Plot
# Extract what is before the first "_"
first <- sub("_.*", "", ggg)
# Extract what is after the second "_"
second <- sub(".*?_(.*?_)", "", ggg)
second <- sub("\\.txt$", "", second)  #
first <- sub("Yes", "Heat root", first)
first <- sub("No","Not Heat root", first)
second <- sub("Yes", "Heat root", second)
second <- sub("No","Not Heat root", second)
if (i == "negative")
{ titolo <- paste0("TOP GO terms overexpressed in ", second, " compared to ", first)
}else if (i =="positive") {
titolo <- paste0("TOP GO terms overexpressed in ", first, " compared to ", second)
} else if (i =="all" & ddd == "Soil") {
titolo <- "Soil"
} else if (i =="all" & ddd == "Ster_root") { titolo <- "Heat Root"}
top_go$class <- str_wrap(top_go$class, width = 30)  # Wrap text at 30 characters
library(ggplot2)
library(grid)  # For precise text unit scaling

dot_plot <- ggplot(top_go, aes(x = GeneRatioVal, y = reorder(class, GeneRatioVal), size = BgRatioVal, color = classfdr)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(2, 10)) +
  labs(
    x = "Gene Ratio",
    y = "GO Terms",
    color = "FDR",
    size = "Background Ratio",
    title = titolo
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9, face = "bold", color = "black")
  ) 

nameplot <- gsub(".txt", "", ggg)

# Use ggsave to ensure consistent font scaling
ggsave(filename = paste0("plot/", i, "_GO_dotplot_stab_", nameplot, ".pdf"),
       plot = dot_plot, dpi = 300, units = "in")
### 3. Enrichment Map (Network Visualization)
# Example edge list for GO term relationships based on shared genes
# Assuming shared genes information is available in TERM2GENE
go_edges <- TERM2GENE %>%
  filter(Description %in% top_go$class) %>%
  group_by(Description) %>%
  summarize(Genes = list(Gene)) %>%
  tidyr::unnest(cols = Genes) %>%
  distinct() %>%
  inner_join(., ., by = "Genes") %>%
  filter(Description.x != Description.y) %>%
  group_by(Description.x, Description.y) %>%
  summarise(weight = n(), .groups = "drop") %>%
  rename(from = Description.x, to = Description.y)
# Create a graph object
go_graph <- graph_from_data_frame(go_edges, directed = FALSE)
# Enrichment Map visualization
enrichment_map <- ggraph(go_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
  geom_node_point(size = 5, color = "blue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  labs(title = "GO Enrichment Map") +
  theme_void()
nameplot <- gsub(".txt","",ggg)
pdf(paste0("plot/",i,"_enrichment_map_",nameplot,".pdf"))
print(enrichment_map)
dev.off()
#############################################################################################
#Remove lines without KEGG pathway
bdata<-bdata[bdata$pKEGG_class!=""|bdata$pKEGG_map!="",]
#Create the universe of KEGG classes and maps
allkeggclass<-unlist(strsplit(bdata$pKEGG_class,"\\|"))
allkeggclass<-allkeggclass[allkeggclass!=""]
allkeggmap<-unlist(strsplit(bdata$pKEGG_map,"\\|"))
allkeggmap<-allkeggmap[allkeggmap!=""]
#indata$mappval<-indata$mapbg<-indata$mapratio<-indata$keggmap<-indata$classpval<-indata$classbg<-indata$classratio<-indata$keggclass<-""
#allclassp<-allfamp<-NULL
mykeggclass<-unlist(strsplit(sigdata$pKEGG_class,"\\|"))
	mykeggclass<-mykeggclass[mykeggclass!=""]
	mykeggmap<-unlist(strsplit(sigdata$pKEGG_map,"\\|"))
	mykeggmap<-mykeggmap[mykeggmap!=""]
	myuniqclass<-unique(mykeggclass)
#Loop over significant results
print(ggg)
pino <- data.frame(class = myuniqclass, classratio = NA, classbg = NA, classp = NA)
for(bbb in 1:length(myuniqclass))
{
			cat("first bbb=",bbb,"\n")
			set1<-sum(mykeggclass==myuniqclass[bbb],na.rm=T)
			set0<-length(mykeggclass)
			bg1<-sum(allkeggclass==myuniqclass[bbb],na.rm=T)
			bg0<-length(allkeggclass)
			pino$classratio[bbb]<-paste(set1,set0,sep="/")
			pino$classbg[bbb]<-paste(bg1,bg0,sep="/")
			pino$classp[bbb]<-fisher.test(matrix(c(set1,set0,bg1,bg0),nrow=2))$p.value
			#pino$ES <- enrichmentScore(pino$classratio[bbb], pino$classbg[bbb])
			#Only write entries with at least minpos terms in the same KEGG category and with a p.value for enrichment lower than p.value 
			
	}
pino$classfdr<-p.adjust(pino$classp,"fdr")
pino <- pino[order(pino$classfdr), ]
 write.table(pino,paste0(i,"_KEGG_", ggg),quote=F, row.names=F)
 }
}

if (ddd == "Soil")
{
# Add a new column 'Comparison' to 'alldata' with the value of 'name'
alldata$Comparison <- name
# Create an empty dataframe with the same column names as 'alldata'
# Append 'alldata' to the empty dataframe
final_dataframe <- rbind(final_dataframe, alldata)

}

}
if (ddd == "Soil")
{
filtered_2_dataframe <- final_dataframe %>%
  group_by(class) %>%            # Group by the 'class' column
  filter(n() == 2) %>%           # Keep rows where the count of 'class' is exactly 3
  ungroup() 
filtered_3_dataframe <- final_dataframe %>%
  group_by(class) %>%            # Group by the 'class' column
  filter(n() == 3) %>%           # Keep rows where the count of 'class' is exactly 3
  ungroup() 
  # Remove duplicate rows based on the 'class' column
  
filtered_3_dataframe <- distinct(filtered_3_dataframe, class, .keep_all = TRUE)
write.xlsx(final_dataframe, paste0(ddd,"_every_significative.xlsx"), row.names = F)
write.xlsx(filtered_2_dataframe, paste0(ddd,"_significative_in_2_comparisons.xlsx"), row.names = F)
write.xlsx(filtered_3_dataframe, paste0(ddd,"_significative_in_all_comparisons.xlsx"), row.names = F)
complete_dataframe <- as.data.frame(final_dataframe$class)
names(complete_dataframe) <- c("class")
complete_dataframe <- as.data.frame(complete_dataframe[!duplicated(complete_dataframe$class), ])
names(complete_dataframe) <- c("class")
forgototal <- as.data.frame(bdata[,1:2])
forgosoil <- merge(complete_dataframe, forgototal, by.x = "class", by.y = "transcript")
listsoil <- forgosoil$uniprot
 soil_enrichment <- enricher(
  gene = listsoil,
  TERM2GENE = TERM2GENE,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
# Convert soil_enrichment to a data frame
soil_results <- as.data.frame(soil_enrichment)
soil_results$Description <- NULL
formatted_soil <- soil_results %>%
  left_join(TERM2GENE %>% dplyr::select(Term, Description) %>% distinct(), by = c("ID" = "Term")) %>%
  mutate(
    # Split GeneRatio and BgRatio into numeric numerator/denominator
    GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatioDen = as.numeric(sub(".*/", "", GeneRatio)),
    GeneRatioVal = GeneRatioNum / GeneRatioDen,
    BgRatioNum = as.numeric(sub("/.*", "", BgRatio)),
    BgRatioDen = as.numeric(sub(".*/", "", BgRatio)),
    BgRatioVal = BgRatioNum / BgRatioDen
  ) %>%
  # Select and rename columns for clarity
  transmute(
    class = Description,                  # GO term description
    classratio = paste(GeneRatioNum, "/", GeneRatioDen, sep = ""),
    classbg = paste(BgRatioNum, "/", BgRatioDen, sep = ""),
    classp = pvalue,                      # Raw p-value
    classfdr = p.adjust                   # Adjusted p-value (FDR)
  )
  write.table(formatted_soil,paste0(i,"_GO_soil.txt"),quote=F, row.names=F)
# Select top GO terms based on FDR for bar plot and dot plot
top_go <- formatted_soil %>% 
  arrange(classfdr) %>% 
  slice_head(n = 20)  # Select top 20 terms
   # Select top 20 terms
# Ensure columns are properly formatted for calculations
top_go <- top_go %>%
  mutate(
    GeneRatioVal = as.numeric(sub("/.*", "", classratio)) / as.numeric(sub(".*/", "", classratio)),
    BgRatioVal = as.numeric(sub("/.*", "", classbg)) / as.numeric(sub(".*/", "", classbg))
  )
### 1. Bar Plot
bar_plot <- ggplot(top_go, aes(x = reorder(class, -classfdr), y = GeneRatioVal, fill = classfdr)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    x = "GO Terms",
    y = "Gene Ratio",
    fill = "FDR",
    title = "Enriched GO Terms (Bar Plot)"
  ) +
  theme_minimal()
# Print the bar plot
dir.create("plot")
nameplot <- "soil"
pdf("plot/GO_barplot_soil.pdf")
print(bar_plot)
dev.off()
### 2. Dot Plot
top_go$class <- str_wrap(top_go$class, width = 30)  # Wrap text at 30 characters
dot_plot <- ggplot(top_go, aes(x = GeneRatioVal, y = reorder(class, GeneRatioVal), size = BgRatioVal, color = classfdr)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(0, 0.09)) +  # Fixed x-axis limits
  scale_size(range = c(2, 10)) +
  labs(
    x = "Gene Ratio",
    y = "GO Terms",
    color = "FDR",
    size = wrap_legend_title("Background Ratio"),
    title = titolo
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),  # Fix text size to real-world points
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  ) +
  coord_fixed(ratio = 0.0097)+  # Keep x-axis width consistent
  guides(
 color = guide_colorbar(order = 1),  # Keep FDR as a gradient bar
    size = guide_legend(order = 2)  # Ensure "Background Ratio" (size legend) is second
  )

nameplot <- gsub(".txt", "", ggg)

# Use ggsave to ensure consistent font scaling
ggsave(filename = "plot/GO_dotplot_soil.pdf",
       plot = dot_plot, width = 8, height = 6, dpi = 300, units = "in")
### 3. Enrichment Map (Network Visualization)
# Example edge list for GO term relationships based on shared genes
# Assuming shared genes information is available in TERM2GENE
go_edges <- TERM2GENE %>%
  filter(Description %in% top_go$class) %>%
  group_by(Description) %>%
  summarize(Genes = list(Gene)) %>%
  tidyr::unnest(cols = Genes) %>%
  distinct() %>%
  inner_join(., ., by = "Genes") %>%
  filter(Description.x != Description.y) %>%
  group_by(Description.x, Description.y) %>%
  summarise(weight = n(), .groups = "drop") %>%
  rename(from = Description.x, to = Description.y)
# Create a graph object
go_graph <- graph_from_data_frame(go_edges, directed = FALSE)
# Enrichment Map visualization
enrichment_map <- ggraph(go_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
  geom_node_point(size = 5, color = "blue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  labs(title = "GO Enrichment Map") +
  theme_void()
nameplot <- gsub(".txt","",ggg)
pdf("plot/GO_enrichment_map_soil.pdf")
print(enrichment_map)
dev.off()

complete_dataframe$Soil <- 1
}
if (ddd != "Soil" && ddd != "Ster_soil")
{
alldata_reduced <- as.data.frame(alldata$class)
names(alldata_reduced) <- c("class")
alldata_reduced[[ddd]] <- 1
complete_dataframe <- merge(complete_dataframe,alldata_reduced, by = "class", all = TRUE)
}
}
complete_dataframe[is.na(complete_dataframe)] <- 0
complete_dataframe$Ster_soil <- 0
write.xlsx(complete_dataframe, "/home/massimo.guazzini/transcriptomic/DESeq2/comparisons_complete_conditions.xlsx", rownames = F)
rownames(complete_dataframe) <- complete_dataframe$class
complete_dataframe$class <- NULL
pdf("/home/massimo.guazzini/transcriptomic/DESeq2/complete_upset.pdf")
# Generate the UpSet plot
ups <- upset(
  complete_dataframe,
  sets = colnames(complete_dataframe),
  sets.bar.color = "#56B4E9",
  main.bar.color = "#009E73",
  matrix.color = "#E69F00",
  order.by = "freq"
)
print(ups)
dev.off()
}
kegg.enrich(infile=infile,keggfile=keggfile,outfile=outfile)
