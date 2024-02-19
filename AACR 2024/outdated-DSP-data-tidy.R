# 01: DSP data import and tidy


# Load packages
packages <- c("tidyverse", "readxl", "utils", "stats", "BiocManager", "NanoStringNCTools", "GeomxTools", "GeoMxWorkflows", "Seurat", "SpatialExperiment", "standR", "edgeR", "ComplexHeatmap", "SpatialDecon", "EnhancedVolcano", "GSVA")
lapply(packages, library, character.only = TRUE)

# Read GeoMx Q3-normalized count & annotation excel data
count_data_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Q3_Normalization.xlsx", sheet = "TargetCountMatrix") |> 
  as.data.frame()
seg_properties_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Q3_Normalization.xlsx", sheet = "SegmentProperties") |> 
  as.data.frame()
target_properties_base <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Q3_Normalization.xlsx", sheet = "TargetProperties") |> 
  as.data.frame()

# Read excel data into SpatialExperiment object (standR package)
spe_standR <- readGeoMx(count_data_base, seg_properties_base, target_properties_base)


# GBM subtyping algorithm (created by me)

# Z-normalize count data to prepare for subtyping
spe_standR@assays@data@listData$zcounts <- scale(spe_standR@assays@data@listData[["counts"]])

# Read subtype profile matrix from Wang et al.
subtype_profile_matrix <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Reorganized") |> 
  as.data.frame()

# Cut and paste gene target names to row names
rownames(subtype_profile_matrix) <- subtype_profile_matrix$TargetName
subtype_profile_matrix <- subtype_profile_matrix |>
  select(!TargetName)

# Filter rows by concordant genes in Z-normalized count data & subtype profile matrices
spe_standR@assays@data@listData$zcounts <- spe_standR@assays@data@listData$zcounts |>
  as.data.frame() |>
  filter(rownames(spe_standR@assays@data@listData$zcounts) %in% rownames(subtype_profile_matrix))
subtype_profile_matrix <- subtype_profile_matrix |>
  filter(rownames(subtype_profile_matrix) %in% rownames(spe_standR@assays@data@listData$zcounts))

# Match order of gene names across matrices
spe_standR@assays@data@listData$zcounts <- spe_standR@assays@data@listData$zcounts[match(rownames(subtype_profile_matrix), rownames(spe_standR@assays@data@listData$zcounts)),]

# Transpose subtype profile matrix
subtype_profile_matrix <- subtype_profile_matrix |>
  t()

# Create new matrix with scores for Mesenchymal, Proneural, and Classical subtypes for each ROI
MES <- apply(X = spe_standR@assays@data@listData[["zcounts"]], MARGIN = 2, FUN = `%*%`, subtype_profile_matrix["Class1",])
PN <- apply(X = spe_standR@assays@data@listData[["zcounts"]], MARGIN = 2, FUN = `%*%`, subtype_profile_matrix["Class2",])
CL <- apply(X = spe_standR@assays@data@listData[["zcounts"]], MARGIN = 2, FUN = `%*%`, subtype_profile_matrix["Class3",])
spe_standR@assays@data@listData$subtypescores <- rbind(MES, PN, CL)

# Get numeric indices for maximum scores within each ROI
max_score <- apply(X = spe_standR@assays@data@listData$subtypescores, MARGIN = 2, FUN = max)
max_index <- match(max_score, spe_standR@assays@data@listData$subtypescores)

# Initialize empty character vector to store assigned subtypes
Subtype <- character()

# Fill "Subtype" vector
for (i in 1:length(max_index)) {
  temp <- max_index[i] %% 3
  if (temp == 1) {
    Subtype <- append(Subtype, rownames(spe_standR@assays@data@listData$subtypescores)[1])
  }
  else if (temp == 2) {
    Subtype <- append(Subtype, rownames(spe_standR@assays@data@listData$subtypescores)[2])
  }
  else if (temp == 0) {
    Subtype <- append(Subtype, rownames(spe_standR@assays@data@listData$subtypescores)[3])
  }
}

# Append "Subtype" vector to appropriate list within spe object
spe_standR@colData@listData$Subtype <- Subtype


# GBM subtyping algorithm (GSVA)

# Read subtype gene sets (3) from Wang et al.
gene_set_MES <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Subtype Signatures", range = "A5:A54", col_names = FALSE) |>
  pull()
gene_set_PN <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Subtype Signatures", range = "I5:I54", col_names = FALSE) |>
  pull()
gene_set_CL <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Subtype Signatures", range = "Q5:Q54", col_names = FALSE) |>
  pull()

# Create gene list from gene sets
gs <- list(gene_set_MES, gene_set_PN, gene_set_CL)
names(gs) <- c("MES", "PN", "CL")

# subtype_profile_matrix <- subtype_profile_matrix |>
#   filter(rownames(subtype_profile_matrix) %in% rownames(spe_standR@assays@data@listData$zcounts))

# GSVA
par <- gsvaParam(as.matrix(spe_standR@assays@data@listData[["counts"]]), gs)
es <- gsva(par)


# Create more operable data frame from spe object
dsp_df <- spe_standR@colData@listData |>
  as.data.frame() |>
  select(c("PatientID", "SpotID", "Type", "VascularProximity")) |>
  slice(rep(1:n(), each = spe_standR@rowRanges@elementMetadata@nrows))
counts_temp <- spe_standR@assays@data@listData[["counts"]] |>
  arrange(rownames(spe_standR@assays@data@listData[["counts"]]))
Gene <- rownames(counts_temp) |>
  rep(ncol(counts_temp))
counts_temp <- counts_temp |>
  stack() |>
  cbind(Gene)
colnames(counts_temp) <- c("Q3_Normalized", "Source", "Gene")

dsp_df <- cbind(dsp_df, counts_temp)
dsp_df <- dsp_df |>
  relocate(Source,.after = SpotID) |>
  relocate(Gene, .before = Q3_Normalized)


# Figure: counts per ROI (cpROI)
cpROI <- dsp_df |>
  ggplot(mapping = aes(x = Source, y = Q3_Normalized, fill = Type)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10')
cpROI

