# Main GeoMx Data Set Analysis

library(tidyverse)
library(readxl)
library(utils)
library(BiocManager)
library(ComplexHeatmap)
library(SpatialDecon)
library(EnhancedVolcano)

# read excel data
geomx_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/NormTargQ3_1_AL.xlsx", sheet="TargetCountMatrix")
geomx_data_base <- data.frame(geomx_data_base)
annotation_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/TMA Annotation File (EDIT).xlsx")
annotation_data_base <- data.frame(annotation_data_base)
