library(readxl)
library(dplyr)
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

# a1: geomx_data and hu_data (unfiltered)
a1_hu_unfiltered_edge <- list_concordant_genes(geomx_data_edge, hu_data_base)
a1_hu_unfiltered_core <- list_concordant_genes(geomx_data_core, hu_data_base)

# read segment properties for hu_data
hu_data_segment_properties <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/hu_brain_count_results/Export4_NormalizationQ3.xlsx", sheet="SegmentProperties")
hu_data_segment_properties <- data.frame(hu_data_segment_properties)
ID <- seq(1, 252, 1)
hu_data_segment_properties <- cbind(ID, hu_data_segment_properties)

# filter for segments that are cortical (and not hippocampal)
hu_data_subset_hippo <- hu_data_base
colnames(hu_data_subset_hippo)[-1] <- seq(1, 252, 1)
hu_data_segment_properties <- filter(hu_data_segment_properties, hu_data_segment_properties$Region == "Cortex")
hu_data_subset_hippo <- hu_data_subset_hippo[,(names(hu_data_subset_hippo) %in% hu_data_segment_properties$ID)]
hu_data_subset_hippo <- cbind(hu_data_base[1], hu_data_subset_hippo)

# a2: geomx_data and hu_data (without hippocampal regions)
a2_hu_subset_hippo_edge <- list_concordant_genes(geomx_data_edge, hu_data_subset_hippo)
a2_hu_subset_hippo_core <- list_concordant_genes(geomx_data_core, hu_data_subset_hippo)


# filter for segments that are cortical and NOT neuronal or neuropilic
hu_data_segment_properties <- filter(hu_data_segment_properties, hu_data_segment_properties$Iba1 == "True" | hu_data_segment_properties$GFAP == "True")
hu_data_subset_neuronal <- hu_data_base
colnames(hu_data_subset_neuronal)[-1] <- seq(1, 252, 1)
hu_data_subset_neuronal <- hu_data_subset_neuronal[,(names(hu_data_subset_neuronal) %in% hu_data_segment_properties$ID)]
hu_data_subset_neuronal <- cbind(hu_data_base[1], hu_data_subset_neuronal)

# a3: geomx_data and hu_data (only GFAP+/Iba1+ regions)
a3_hu_subset_neuronal_edge <- list_concordant_genes(geomx_data_edge, hu_data_subset_neuronal)
a3_hu_subset_neuronal_core <- list_concordant_genes(geomx_data_core, hu_data_subset_neuronal)


