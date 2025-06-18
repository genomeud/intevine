#!/bin/bash
#This script produce all the graphs and analysis presented in the paper: "Variations in root-soil system influence Grapevine Holobiont by Shaping Plant Physiology and Root Microbiome".  The crucial point is to put the 
# metabarcoding reads (avaiable at NCBI Sequence Read Archive under the accession PRJNA1217140) in the READSDIR folder (separating T0 from T1 and 16s from ITS) and the transcriptomic reads  (avaiable at in GEO under the accession GSE287854)
#in the READDIR folder. You also need the metadata file (avaiable in this github) as well as the databases for SILVA, UNITE,  the Vitis vinifera 12Xv0 genome assembly of the strain PN40024 (GCA_000003745.2), FAPROTAX and FUNGUILD, which are publicly avaiable.

# Define the base path for the R scripts
BASE_PATH="/projects/marroni/functions/intevine/new_pipeline/"

source /iga/scripts/dev_modules/miniconda3.test/bin/activate
conda deactivate
conda activate r_4.1.3

#!/bin/bash


#####################DADA2 CLASSIFICATION############
# -------------
# Base R script paths
# -------------
BASE_PATH="/projects/marroni/functions/intevine/new_pipeline"
DADA16ST0="${BASE_PATH}/00_dada2_silva_T0.r"
DADA16T1="${BASE_PATH}/00_dada2_silva"
DADA2ITST0="${BASE_PATH}/00_dada2_T0_unite.r"
DADA2ITS="${BASE_PATH}/00_dada2_unite.r"

# -------------
# 16S T1
# -------------
READSDIR16T1="/projects/marroni/intevine/sequences/dna/16s/T1/Novaseq/"
BASEDIR16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1_testp/"
DBDIR16="/projects/marroni/databases/dada2/"
FILEPATTERN16=""

# -------------
# 16S T0
# -------------
READSDIR16T0="/projects/marroni/intevine/sequences/dna/16s/T0/Novaseq/"
BASEDIR16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/"

# -------------
# ITS T1
# -------------
READSDIRITST1="/projects/marroni/intevine/sequences/dna/ITS/T1/"
BASEDIRITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/"
DBDIRITS="/projects/marroni/intevine/database/dada2/"
FILEPATTERNITS=""

# -------------
# ITS T0
# -------------
READSDIRITST0="/projects/marroni/intevine/sequences/dna/ITS/T0/"
BASEDIRITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/"

# --------------------
# Execute each R script
# --------------------

echo "Running 16S T1..."
Rscript "$DADA16T1" \
  --readsdir "$READSDIR16T1" \
  --basedir "$BASEDIR16T1" \
  --DBdir "$DBDIR16" \
  --previously_run FALSE \
  --max.avs 10000 \
  --use_silva TRUE \
  --use_rdp FALSE \
  --use_greengenes FALSE \
  --use_GTDB FALSE \
  --filepattern "$FILEPATTERN16"

echo "Running 16S T0..."
Rscript "$DADA16ST0" \
  --readsdir "$READSDIR16T0" \
  --basedir "$BASEDIR16T0" \
  --DBdir "$DBDIR16" \
  --previously_run FALSE \
  --max.avs 10000 \
  --use_silva TRUE \
  --use_rdp FALSE \
  --use_greengenes FALSE \
  --use_GTDB FALSE \
  --filepattern "$FILEPATTERN16"

echo "Running ITS T1..."
Rscript "$DADA2ITS" \
  --readsdir "$READSDIRITST1" \
  --basedir "$BASEDIRITST1" \
  --DBdir "$DBDIRITS" \
  --previously_run FALSE \
  --max.avs 10000 \
  --use_unite TRUE \
  --use_rdp FALSE \
  --use_greengenes FALSE \
  --use_GTDB FALSE \
  --filepattern "$FILEPATTERNITS"

echo "Running ITS T0..."
Rscript "$DADA2ITST0" \
  --readsdir "$READSDIRITST0" \
  --basedir "$BASEDIRITST0" \
  --DBdir "$DBDIRITS" \
  --previously_run FALSE \
  --max.avs 10000 \
  --use_unite TRUE \
  --use_rdp FALSE \
  --use_greengenes FALSE \
  --use_GTDB FALSE \
  --filepattern "$FILEPATTERNITS"

echo "classification completed"


#################### RAREFACTION CURVES ###############

#!/bin/bash

# ---------------------------
# Path to R script
# ---------------------------
RAREFACTION_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/01_rarefaction.r"

# ---------------------------
# Abundance table inputs
# ---------------------------
ABUNDANCE16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/_Species_merged_count.txt"
ABUNDANCE16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/_Species_merged_count.txt"
ABUNDANCEITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/_Species_merged_count.txt"
ABUNDANCEITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/_Species_merged_count.txt"

# ---------------------------
# Output directories
# ---------------------------
OUTPLOT16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/07_plot_SILVA/"
OUTPLOT16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/07_plot_SILVA/"
OUTPLOTITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/07_plot_unite/"
OUTPLOTITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/07_plot_unite/"

# ---------------------------
# Execute the script
# ---------------------------
echo "Running rarefaction analysis..."

Rscript "$RAREFACTION_SCRIPT" \
  --abundance "$ABUNDANCE16T1" \
  --abundanceT0 "$ABUNDANCE16T0" \
  --abundanceITS "$ABUNDANCEITST1" \
  --abundanceITST0 "$ABUNDANCEITST0" \
  --output "$OUTPUT16T1" \
  --outputT0 "$OUTPUT16T0" \
  --outputITS "$OUTPUTITST1" \
  --outputITST0 "$OUTPUTITST0"

echo "Rarefaction analysis completed."



######### NORMALIZATION AND AGGREGATION TO GENUS LEVEL, the results of these script will be used as the input of the other scripts

# ---------------------------
# R script path
# ---------------------------
PERCENTAGE_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/02_percentage_pipeline.r"

# ---------------------------
# Output directories
# ---------------------------
OUTPUT16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/"
OUTPUT16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/"
OUTPUTITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/"
OUTPUTITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/"

# ---------------------------
# Count fragments
# ---------------------------
COUNT16T1="/projects/marroni/intevine/docs/count_reads_16s_output_T1.txt"
COUNT16T0="/projects/marroni/intevine/docs/count_reads_16s_T0_output.txt"
COUNTITST1="/projects/marroni/intevine/docs/count_reads_ITS_output_T1.txt"
COUNTITST0="/projects/marroni/intevine/docs/count_reads_ITS_T0_output.txt"

# ---------------------------
# Run the R script
# ---------------------------
echo "Running percentage pipeline..."

Rscript "$PERCENTAGE_SCRIPT" \
  --abundanceT1 "$ABUNDANCE16T1" \
  --abundanceT0 "$ABUNDANCE16T0" \
  --abundanceITST1 "$ABUNDANCEITST1" \
  --abundanceITST0 "$ABUNDANCEITST0" \
  --output "$OUTPUT16T1" \
  --outputT0 "$OUTPUT16T0" \
  --outputITS "$OUTPUTITST1" \
  --outputITST0 "$OUTPUTITST0" \
  --countfragments "$COUNT16T1" \
  --countfragmentsT0 "$COUNT16T0" \
  --countfragmentsITS "$COUNTITST1" \
  --countfragmentsITST0 "$COUNTITST0"

echo "Percentage pipeline completed."

# ---------------------------
# Path to R script
# ---------------------------
GENUS_AGG_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/03_genus_aggregation.r"

# ---------------------------
# Inputs from previous percentage script
# ---------------------------
PERCENTAGE16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/0.001_average_filtered_percentage_counts_edited_16sT0.xlsx"
PERCENTAGE16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/0.001_average_filtered_percentage_counts_edited_16sT1.xlsx"
PERCENTAGEITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_silva/0.001_average_filtered_percentage_counts_edited_ITST0.xlsx"
PERCENTAGEITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_silva/0.001_average_filtered_percentage_counts_edited_ITST1.xlsx"

# ---------------------------
# Run the script
# ---------------------------
echo "Running genus aggregation..."

Rscript "$GENUS_AGG_SCRIPT" \
  --abundanceT0 "$PERCENTAGE16T0" \
  --abundanceT1 "$PERCENTAGE16T1" \
  --abundanceITST0 "$PERCENTAGEITST0" \
  --abundanceITST1 "$PERCENTAGEITST1" \
  --outputT0 "$OUTPUT16T0" \
  --outputT1 "$OUTPUT16T1" \
  --outputITST0 "$OUTPUTITST0" \
  --outputITST1 "$OUTPUTITST1"

echo "Genus aggregation completed."


################ METADATA ANALYSIS #########

#!/bin/bash

# ---------------------------
# Path to R script
# ---------------------------
COLLINEARITY_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/04_check_collinearity.r"

# ---------------------------
# Input metadata file
# ---------------------------
METAFILE="/projects/marroni/intevine/docs/map_grapevine_final.xlsx"
# ---------------------------
# Output directory
# ---------------------------
OUTPUT_COLLINEARITY="/home/massimo.guazzini/intevine/analysis/meta/new_pipeline/collinearity/"

# ---------------------------
# Run the script
# ---------------------------
echo "Running collinearity check..."

Rscript "$COLLINEARITY_SCRIPT" \
  --metafile "$METAFILE" \
  --output "$OUTPUT_COLLINEARITY"

echo "Collinearity check completed."

#!/bin/bash

# ---------------------------
# Path to R script
# ---------------------------
BOXPLOT_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/05_meta_boxplot_analysis.r"

# ---------------------------
# Output directory
# ---------------------------
OUTPUT_BOXPLOTS="/projects/marroni/intevine/analysis/meta/"

# ---------------------------
# Conditions for analysis
# ---------------------------
CONDITIONS="Soil,Autoclave,Heat root"

# ---------------------------
# Run the script
# ---------------------------
echo "Running metadata boxplot analysis..."

Rscript "$BOXPLOT_SCRIPT" \
  --metafile "$METAFILE" \
  --out "$OUTPUT_BOXPLOTS" \
  --conditions "$CONDITIONS"

echo "Boxplot analysis completed."

#!/bin/bash

# ---------------------------
# Path to R script
# ---------------------------
NMDS_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/06_nmds.r"


# ---------------------------
# Run the script
# ---------------------------
echo "Running NMDS analysis..."

Rscript "$NMDS_SCRIPT" \
  --abundanceT0 "$PERCENTAGE16T0" \
  --abundanceT1 "$PERCENTAGE16T1" \
  --abundanceITST0 "$PERCENTAGEITST0" \
  --abundanceITST1 "$PERCENTAGEITST1" \
  --metafile "$METAFILE" \
  --outputT0 "$OUTPLOT16T0" \
  --outputT1 "$OUTPLOT16T1" \
  --outputITST0 "$OUTPLOTITST0" \
  --outputITST1 "$OUTPLOTITST1" \
  --outtableT0 "$OUTPUT16T0" \
  --outtableT1 "$OUTPUT16T1" \
  --outtableITST0 "$OUTPUTITST0" \
  --outtableITST1 "$OUTPUTITST1" 

echo "NMDS analysis completed." REMOVE GENESCORR

#!/bin/bash

# ---------------------------
# Path to R script
# ---------------------------
STACKED_BAR_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/07_stacked_bar_plot.r"

# ---------------------------
# Input genus-aggregated abundance files
# ---------------------------
GENUS16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/genus_aggregated.xlsx"
GENUS16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/genus_aggregated.xlsx"
GENUSITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_silva/genus_aggregated.xlsx"
GENUSITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_silva/genus_aggregated.xlsx"




# ---------------------------
# Gene correlation threshold
# ---------------------------
GENES_CORR=25

# ---------------------------
# Run the script
# ---------------------------
echo "Running stacked bar plot generation..."

Rscript "$STACKED_BAR_SCRIPT" \
  --abundanceT0 "$GENUS16T0" \
  --abundanceT1 "$GENUS16T1" \
  --abundanceT0 "$GENUSITST0" \
  --abundanceT1 "$GENUSITST1" \
  --metafile "$METAFILE" \
  --outputT0 "$OUTPLOT16T0" \
  --outputT1 "$OUTPLOT16T1" \
  --outputITST0 "$OUTPLOTITST0" \
  --outputITST1 "$OUTPLOTITST1" \
  --genescorr "$GENES_CORR"

echo "Stacked bar plot generation completed."


# This script prepares the raw ASV count tables used for
# alpha diversity and DESeq2 differential abundance analysis.

# ---------------------------
# Path to R script
# ---------------------------
RAW_ASV_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/08_raw_asv_after_filtering.r"

# ---------------------------
# Run the script
# ---------------------------
echo "Running raw ASV export for alpha diversity and DESeq2..."

Rscript "$RAW_ASV_SCRIPT" \
  --abundanceT0 "$PERCENTAGE16T0" \
  --abundanceT1 "$PERCENTAGE16T1" \
  --abundanceITST0 "$PERCENTAGEITST0" \
  --abundanceITST1 "$PERCENTAGEITST1" \
  --raw_allT0 "$ABUNDANCE16T0" \
  --raw_allT1 "$ABUNDANCE16T1" \
  --raw_allITST0 "$ABUNDANCEITST0" \
  --raw_allITST1 "$ABUNDANCEITST1" \
  --outtableT0 "$OUTPUT16T0" \
  --outtableT1 "$OUTPUT16T1" \
  --outtableITST0 "$OUTPUTITST0" \
  --outtableITST1 "$OUTPUTITST1"

echo "Raw ASV export completed."

# This script extracts and compares microbiome samples using PERMANOVA
# based on filtered abundance tables and selected metadata conditions.
# Useful for identifying significant group differences in beta diversity.

# ---------------------------
# Path to R script
# ---------------------------
PERMANOVA_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/09_permanova.r"


# ---------------------------
# Output folder
# ---------------------------
OUTPUT_PERMANOVA="/projects/marroni/intevine/analysis/metagenomic/permanova/"

# ---------------------------
# Run the script
# ---------------------------
echo "Running PERMANOVA analysis..."

Rscript "$PERMANOVA_SCRIPT" \
  --abundanceT0 "$PERCENTAGE16T0" \
  --abundanceT1 "$PERCENTAGE16T1" \
  --abundanceITST0 "$PERCENTAGEITST0" \
  --abundanceITST1 "$PERCENTAGEITST1" \
  --metafile "$METAFILE" \
  --token "$TOKEN" \
  --out "$OUTPUT_PERMANOVA"

echo "PERMANOVA analysis completed."


# This script calculates alpha diversity indices from raw ASV counts
# using filtered and raw abundance tables. Useful for comparing richness
# and evenness across conditions using metadata.

# ---------------------------
# Path to R script
# ---------------------------
INDICES_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/10_indices.r"

# ---------------------------
# Input: raw counts after filtering (used for alpha diversity)
# ---------------------------
RAWTABLE16T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/0.001_raw_counts_edited.xlsx"
RAWTABLE16T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/0.001_raw_counts_edited.xlsx"
RAWTABLEITST0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/06_tables_unite/0.001_raw_counts_edited.xlsx"
RAWTABLEITST1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/0.001_raw_counts_edited.xlsx"

# ---------------------------
# Conditions for grouping
# ---------------------------
CONDITIONS_T0="Soil,Ster_soil"
CONDITIONS_T1="Soil,Ster_soil,Ster_root"

# ---------------------------
# Run the script
# ---------------------------
echo "Running alpha diversity index calculation..."

Rscript "$INDICES_SCRIPT" \
  --abundanceT0 "$RAWTABLE16T0" \
  --abundanceT1 "$RAWTABLE16T1" \
  --abundanceITST0 "$RAWTABLEITST0" \
  --abundanceITST1 "$RAWTABLEITST1" \
  --metafile "$METAFILE" \
  --outputT0 "$OUTPLOT16T0" \
  --outputT1 "$OUTPLOT16T1" \
  --outputITST0 "$OUTPLOTITST0" \
  --outputITST1 "$OUTPLOTITST1" \
  --outtableT0 "$OUTPUT16T0" \
  --outtableT1 "$OUTPUT16T1" \
  --outtableITST0 "$OUTPUTITST0" \
  --outtableITST1 "$OUTPUTITST1" \
  --conditionlistT0 "$CONDITIONS_T0" \
  --conditionlistT1 "$CONDITIONS_T1"

echo "Alpha diversity analysis completed."

###ADD DESEQ2###

# ---------------------------
# Path to R script
# ---------------------------
AGG_DESEQ2_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/11_deseq2_aggregated_genus.r"


# ---------------------------
# Output directories for plots
# ---------------------------
PLOT_OUT_16S_T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/deseq2/Genus_aggregated/plot/"
PLOT_OUT_16S_T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/deseq2/Genus_aggregated/plot/"
PLOT_OUT_ITS_T0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/Genus_aggregated/plot/"
PLOT_OUT_ITS_T1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/Genus_aggregated/plot/"

# ---------------------------
# Output directories for tables
# ---------------------------
TABLE_OUT_16S_T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/deseq2/Genus_aggregated/"
TABLE_OUT_16S_T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/deseq2/Genus_aggregated/"
TABLE_OUT_ITS_T0="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T0/Genus_aggregated/"
TABLE_OUT_ITS_T1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/Genus_aggregated/"

# ---------------------------
# Conditions
# ---------------------------
CONDITION_T0="Soil,Ster_soil"
CONDITION_T1="Soil,Ster_soil,Ster_root"

# ---------------------------
# Lists of relevant microbes
# ---------------------------
RELEVANT_LIST="/home/massimo.guazzini/Supplementary_file_8_list_of_agricultural_relevant_microorganisms.xlsx"
BAC_BEN="Bacteria_beneficial"
BAC_DAM="Damaging_bacteria"
FUN_BEN="Fungal_beneficial"
FUN_DAM="Damaging_fungi"

# ---------------------------
# Run the script
# ---------------------------
echo "Running aggregated DESeq2 genus-level analysis..."

Rscript "$AGG_DESEQ2_SCRIPT" \
  --abundanceT0 "$RAWTABLE16T0" \
  --abundanceT1 "$RAWTABLE16T1" \
  --abundanceITST0 "$RAWTABLEITST0" \
  --abundanceITST1 "$RAWTABLEITST1" \
  --metafile "$METAFILE" \
  --outputT0 "$PLOT_OUT_16S_T0" \
  --outputT1 "$PLOT_OUT_16S_T1" \
  --outputITST0 "$PLOT_OUT_ITS_T0" \
  --outputITST1 "$PLOT_OUT_ITS_T1" \
  --outtableT0 "$TABLE_OUT_16S_T0" \
  --outtableT1 "$TABLE_OUT_16S_T1" \
  --outtableITST0 "$TABLE_OUT_ITS_T0" \
  --outtableITST1 "$TABLE_OUT_ITS_T1" \
  --conditionT0 "$CONDITION_T0" \
  --condition "$CONDITION_T1" \
  --listrelevant "$RELEVANT_LIST" \
  --bacben "$BAC_BEN" \
  --bacdam "$BAC_DAM" \
  --funben "$FUN_BEN" \
  --fundam "$FUN_DAM"

echo "Aggregated genus-level DESeq2 analysis completed."

# ---------------------------
# Path to R script
# ---------------------------
FAPROTAX_PREP_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/12_faprotax_species_file_preparation.r"

# ---------------------------
# Input abundance file (T1 species level)
# ---------------------------
FAPROTAX_T1_SPECIES="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/0.001_filtered_percentage_counts.xlsx"

# ---------------------------
# Output directories
# ---------------------------
FAPROTAX_OUT_T0="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T0/06_tables_silva/"
FAPROTAX_OUT_T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/"

# ---------------------------
# Gene correlation parameter
# ---------------------------
GENESCORR=25

# ---------------------------
# Run the script
# ---------------------------
echo "Running FAPROTAX species file preparation..."

Rscript "$FAPROTAX_PREP_SCRIPT" \
  --abundanceT1species "$FAPROTAX_T1_SPECIES" \
  --outputT0 "$FAPROTAX_OUT_T0" \
  --outputT1 "$FAPROTAX_OUT_T1" \
  --genescorr "$GENESCORR"

echo "FAPROTAX file preparation completed."


# ---------------------------
# Path to R script
# ---------------------------
FAPROTAX_BOXPLOT_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/14_faprotax_boxplot.r"

# ---------------------------
# Input abundance file (T1, FAPROTAX percentages)
# ---------------------------
FAPROTAX_PERCENT_T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/06_tables_silva/percentage_faprotax.xlsx"

# ---------------------------
# Metadata file
# ---------------------------
FAPROTAX_META="/projects/marroni/intevine/docs/map_grapevine_final.xlsx"

# ---------------------------
# Output directory
# ---------------------------
FAPROTAX_OUT_T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/09_full_ASV/07_plot_SILVA/001/"

# ---------------------------
# Condition list (T1)
# ---------------------------
FAPROTAX_CONDITIONS_T1="Soil,Ster_soil,Ster_root"

# ---------------------------
# Run the script
# ---------------------------
echo "Running updated FAPROTAX boxplot analysis..."

Rscript "$FAPROTAX_BOXPLOT_SCRIPT" \
  --abundanceT1 "$FAPROTAX_PERCENT_T1" \
  --metafile "$FAPROTAX_META" \
  --outputT1 "$FAPROTAX_OUT_T1" \
  --conditionlistT1 "$FAPROTAX_CONDITIONS_T1"

echo "Updated FAPROTAX boxplot analysis completed."


funguild 14

# ---------------------------
# Path to R script
# ---------------------------
FUNGUILD_BOXPLOT_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/15_boxplot_funguild.r"

# ---------------------------
# Input abundance file (T1, FUNGuild percentages)
# ---------------------------
FUNGUILD_ABUND_T1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/06_tables_unite/0.001_filtered_percentage_counts_funguild.guilds.reviewed.xlsx"

# ---------------------------
# Metadata file
# ---------------------------
FUNGUILD_META="/projects/marroni/intevine/docs/map_grapevine_final.xlsx"

# ---------------------------
# Output directory
# ---------------------------
FUNGUILD_OUT_T1="/projects/marroni/intevine/analysis/metagenomic/ITS/ASV/T1/09_full_ASV/07_plot_unite/001/"

# ---------------------------
# Condition list (T1)
# ---------------------------
FUNGUILD_CONDITIONS_T1="Soil,Ster_soil,Ster_root"

# ---------------------------
# Run the script
# ---------------------------
echo "Running FUNGuild functional group boxplot analysis..."

Rscript "$FUNGUILD_BOXPLOT_SCRIPT" \
  --abundanceT1 "$FUNGUILD_ABUND_T1" \
  --metafile "$FUNGUILD_META" \
  --outputT1 "$FUNGUILD_OUT_T1" \
  --conditionlistT1 "$FUNGUILD_CONDITIONS_T1"

echo "FUNGuild boxplot analysis completed."


############################# TRANSCRIPTOMIC ###############


#Variables you need to set
BASEDIR=/path/to/your/results #Exiting folder that you want to be the top folder of your results
LOGDIR=/path/to/your/logs
###########################


BASEDIR=/projects/marroni/intevine
REFDIR=/genomes/vitis_vinifera
NEWREF=/projects/marroni/intevine/references/vitis_STAR



#Generate genome including vitis vinifera and candidatus phytoplasma solani
module load it/aligners/star/2.7.2b

STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ${NEWREF} \
--genomeSAindexNbases 12 \
--genomeFastaFiles ${REFDIR}/assembly/reference/12xCHR/vitis_vinifera_12xCHR.fasta \
--sjdbGTFfile ${REFDIR}/genes/assembly_build_12xCHR/V2.1/V2.1.gtf




###############################
#
# Align experiments
#
################################



#In case we want to discover novel transcripts we can follow protocol 8 of the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/#S50title
READDIR=${BASEDIR}/sequences/rna
ALDIR=${BASEDIR}/STAR
NTHREAD=8
for READ1 in ${READDIR}/*_R1_*fastq.gz
do
echo $READ1
READ2=${READ1/_R1_/_R2_}
FNAME=$(basename $READ1)
OUTFILE=${READ1/$READDIR/$ALDIR}
OUTFILE=${OUTFILE/_L003_R1_001.fastq.gz/}
echo $OUTFILE
echo "module load it/aligners/star/2.7.2b; STAR --runThreadN $NTHREAD --genomeDir $NEWREF --outFileNamePrefix $OUTFILE --limitBAMsortRAM 15000000000 \
--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --readFilesIn ${READ1} ${READ2}" \
| qsub -N ST_inte -l vmem=16G,walltime=12:00:00,nodes=1:ppn=$NTHREAD
done

#!/bin/bash

# ---------------------------
# Path to R script
# ---------------------------
COMBINE_RNA_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/17_combine_rna_tables.r"

# ---------------------------
# Input directory with ReadsPerGene files
# ---------------------------
INPUT_DIR="/projects/marroni/intevine/STAR/"

# ---------------------------
# Output file path
# ---------------------------
OUTPUT_FILE="/projects/marroni/intevine/STAR/rna_table.txt"

# ---------------------------
# Run the script
# ---------------------------
echo "Running RNA table combination from STAR outputs..."

Rscript "$17_COMBINE_RNA_SCRIPT" \
  --input_dir "$INPUT_DIR" \
  --output_file "$OUTPUT_FILE"

echo "RNA table combination completed."

# ---------------------------
# Path to R script
# ---------------------------
VIRUS_IDENTIFY_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/18_merge_virusdetect_TPM_FPKM.r"


# ---------------------------
# Input folder
# ---------------------------
VIRUS_INPUT_DIRECTORY="/projects/marroni/intevine/alignment_vv_virus/"

# ---------------------------
# Output folder
# ---------------------------
VIRUS_OUTPUT_DIRECTORY="/projects/marroni/intevine/alignment_vv_virus/tables/RNAseq_VirDet"
# ---------------------------
# Run the script
# ---------------------------
echo "Running PERMANOVA analysis..."

Rscript "$VIRUS_IDENTIFY_SCRIPT" \
  --indir "$VIRUS_INPUT_DIRECTORY" \
  --outfile "$VIRUS_OUTPUT_DIRECTORY" 

# ---------------------------
# Path to R script
# ---------------------------
VIRUS_BOX_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/19_virus_boxplot.r"

# ---------------------------
# Input: Virus abundance table (T1 only)
# ---------------------------
VIRUS_ABUND_T1="/projects/marroni/intevine/alignment_vv_virus/tables/RNAseq_VirDet.xlsx"

# ---------------------------
# Metadata file
# ---------------------------
METAFILE="/projects/marroni/intevine/docs/map_grapevine_final.xlsx"

# ---------------------------
# Output: directory for virus boxplot
# ---------------------------
VIRUS_OUTPUT_T1="/projects/marroni/intevine/analysis/metagenomic/16s/ASV/SILVA/T1/07_plot_SILVA/"

# ---------------------------
# Conditions for grouping (T1)
# ---------------------------
VIRUS_CONDITIONS_T1="Soil,Ster_soil,Ster_root"

# ---------------------------
# Run the script
# ---------------------------
echo "Running virus abundance boxplot analysis..."

Rscript "$VIRUS_BOX_SCRIPT" \
  --abundanceT1 "$VIRUS_ABUND_T1" \
  --metafile "$METAFILE" \
  --outputT1 "$VIRUS_OUTPUT_T1" \
  --conditionlistT1 "$VIRUS_CONDITIONS_T1"

echo "Virus boxplot analysis completed."

# ---------------------------
# Path to R script
# ---------------------------
DESEQ2_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/20_deseq2_analysis.r"

# ---------------------------
# Input: abundance table
# ---------------------------
DESEQ2_ABUND="/projects/marroni/intevine/analysis/transcriptomic/rna_table.txt"

# ---------------------------
# Metadata file
# ---------------------------
DESEQ2_META="/projects/marroni/intevine/docs/metadata.xlsx"

# ---------------------------
# Samples to remove (comma-separated)
# ---------------------------
DESEQ2_REMOVE="L471,L593"

# ---------------------------
# Read count summary file
# ---------------------------
DESEQ2_READS="/mnt/vol1/projects/LUCAS/JRC_result/metadata/countreads.txt"

# ---------------------------
# Output directory for results and plots
# ---------------------------
DESEQ2_OUT="/projects/marroni/intevine/analysis/transcriptomic/DESeq2/"
DESEQ2_GRAPHDIR="/projects/marroni/intevine/analysis/transcriptomic/DESeq2/"

# ---------------------------
# Conditions for grouping
# ---------------------------
DESEQ2_CONDITION="Ster_soil,Root_soil,Bioristor"

# ---------------------------
# Run the script
# ---------------------------
echo "Running DESeq2 transcriptomic analysis..."

Rscript "$DESEQ2_SCRIPT" \
  --abundance "$DESEQ2_ABUND" \
  --metafile "$DESEQ2_META" \
  --removeme "$DESEQ2_REMOVE" \
  --readfile "$DESEQ2_READS" \
  --out "$DESEQ2_OUT" \
  --condition "$DESEQ2_CONDITION" \
  --graphdir "$DESEQ2_GRAPHDIR"

echo "DESeq2 analysis completed."


# ---------------------------
# Path to R script
# ---------------------------
GO_ANALYSIS_SCRIPT="/projects/marroni/functions/intevine/new_pipeline/21_go_analysis.r"

# ---------------------------
# Input: directory or file with DESeq2 results
# ---------------------------
GO_INPUT="/home/massimo.guazzini/transcriptomic/DESeq2/"

# ---------------------------
# KEGG+GO annotation file
# ---------------------------
GO_KEGG_FILE="/home/massimo.guazzini/wmm/Giovanni/new_pipeline/named_transcript_GO_KEGG_path.txt"

# ---------------------------
# Run the script
# ---------------------------
echo "Running GO/KEGG enrichment analysis..."

Rscript "$GO_ANALYSIS_SCRIPT" \
  --infile "$GO_INPUT" \
  --keggfile "$GO_KEGG_FILE"

echo "GO/KEGG analysis completed."
