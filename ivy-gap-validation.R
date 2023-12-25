library(readxl)
library(dplyr)
library(utils)
source("functions.R")

# read excel data
geomx_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/NormTargQ3_1_AL.xlsx", sheet="TargetCountMatrix")
geomx_data_base <- data.frame(geomx_data_base)
annotation_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/TMA Annotation File (EDIT).xlsx")
annotation_data_base <- data.frame(annotation_data_base)
hu_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/hu_brain_count_results/Export4_NormalizationQ3.xlsx", sheet="TargetCountMatrix")
hu_data_base <- data.frame(hu_data_base)

# filter for edge (peripheral) and core regions
edge_regions <- c(2:4, 11:13, 20:21, 28:30, 37:39)
geomx_data_edge <- cbind(geomx_data_base[1], geomx_data_base[edge_regions])
geomx_data_core <- geomx_data_base[-edge_regions]

# read IvyGAP data
ivy_data_base <- read.csv("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/ivy-gap-database/fpkm_table.csv")
ivy_data_gene_names <- read.csv("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/ivy-gap-database/rows-genes.csv")
ivy_data_sample_info <- read.csv("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/ivy-gap-database/columns-samples.csv")

# assign gene names
for (i in (1:length(ivy_data_base$gene_id.rna_well_id))) {
  if (ivy_data_base$gene_id.rna_well_id[i] == ivy_data_gene_names$gene_id[i]) {
    ivy_data_base$gene_id.rna_well_id[i] <- ivy_data_gene_names$gene_symbol[i]
  }
}

# remove "X" from sample names
for (i in (1:ncol(ivy_data_base))) {
  if (substring(colnames(ivy_data_base)[i], 1, 1) == "X") {
    colnames(ivy_data_base)[i] <- sub("X", "", colnames(ivy_data_base)[i])
  }
}

# remove unneeded text from sample classification names
for (i in (1:nrow(ivy_data_sample_info))) {
  ivy_data_sample_info$structure_abbreviation[i] <- sub("-.*", "", ivy_data_sample_info$structure_abbreviation[i])
}

# create data subsets for Core Tumor, Infiltrating Tumor, Leading Edge, and Infiltrating Tumor + Leading Edge
ivy_data_sample_info_CT <- filter(ivy_data_sample_info, ivy_data_sample_info$structure_abbreviation == "CT")
ivy_data_sample_info_IT <- filter(ivy_data_sample_info, ivy_data_sample_info$structure_abbreviation == "IT")
ivy_data_sample_info_LE <- filter(ivy_data_sample_info, ivy_data_sample_info$structure_abbreviation == "LE")
ivy_data_sample_info_ITLE <- filter(ivy_data_sample_info, ivy_data_sample_info$structure_abbreviation == "IT" | ivy_data_sample_info$structure_abbreviation == "LE")

# a1: geomx core and IvyGap CT
ivy_data_CT <- cbind(ivy_data_base[1], ivy_data_base[, names(ivy_data_base) %in% ivy_data_sample_info_CT$rna_well_id])
colnames(ivy_data_CT)[1] <- "TargetName"
a1_ivy_CT <- list_concordant_genes(geomx_data_core, ivy_data_CT)



