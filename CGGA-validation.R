library(readxl)
library(dplyr)
source("functions.R")

# read GeoMx DSP data
geomx_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/NormTargQ3_1_AL.xlsx", sheet="TargetCountMatrix")
geomx_data_base <- data.frame(geomx_data_base)
annotation_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/TMA Annotation File (EDIT).xlsx")
annotation_data_base <- data.frame(annotation_data_base)

# filter for edge (peripheral) and core regions
edge_regions <- c(2:4, 11:13, 20:21, 28:30, 37:39)
geomx_data_edge <- cbind(geomx_data_base[1], geomx_data_base[edge_regions])
geomx_data_core <- geomx_data_base[-edge_regions]

# read CGGA data
clinical_data_base <- read.csv(file="/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/CGGA693-bulkRNA-clinical.txt", sep="\t")
count_data_base <- read.csv(file="/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/CGGA693-bulkRNA-counts.txt", sep="\t")

# filter for primary GBM cases
clinical_data_pGBM <- filter(clinical_data_base, clinical_data_base$PRS_type == "Primary" & clinical_data_base$Histology == "GBM")

# compile CGGA IDs for primary GBM cases
primaryGBM_IDs <- clinical_data_pGBM$CGGA_ID

gene_IDs <- count_data_base$Gene_Name
count_data_pGBM <- subset(count_data_base, select = c(primaryGBM_IDs))
count_data_pGBM <- cbind(TargetName = gene_IDs, count_data_pGBM)

# a1: geomx_data and CGGA_data (primary GBM)
a1_CGGA_pGBM_edge <- list_concordant_genes(geomx_data_edge, count_data_pGBM)
a1_CGGA_pGBM_core <- list_concordant_genes(geomx_data_core, count_data_pGBM)

