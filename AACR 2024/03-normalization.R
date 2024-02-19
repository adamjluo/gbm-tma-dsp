# 03*: Normalization and standard analysis (UMAP, tSNE, PCA, and DE)

# Load relevant packages
packages <- c("DESeq2")
lapply(packages, library, character.only = TRUE)

# Normalize data using NanoString recommended methods

# Q3 normalization
target_dsp_data <- normalize(target_dsp_data,
                             norm_method = "q3",
                             toElt = "Q3_Norm")
# Negative Normalization
target_dsp_data <- normalize(target_dsp_data,
                             norm_method = "neg",
                             averageType = "mean",
                             toElt = "Neg_Norm")

# Housekeeping Normalization
# target_dsp_data <- normalize(target_dsp_data,
#                              norm_method = "hk",
#                              toElt = "HK_Norm")

# Normalize data using other methods (DESeq2, Limma, RUV4, quantile)
# cts <- round(target_dsp_data@assayData[["exprs"]])
# coldata <- target_dsp_data@phenoData@data
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ Type)
# dds <- DESeq(dds)
# res <- results(dds)
# resOrdered <- res[order(res$pvalue), ]





