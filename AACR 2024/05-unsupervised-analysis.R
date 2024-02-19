# Instructions

# Load relevant packages
packages <- c("umap", "Rtsne", "pheatmap", "RColorBrewer", "PoiClaClu")
lapply(packages, library, character.only = TRUE)

Subtype <- target_dsp_data@phenoData@data[["Verhaak et al. Subtype"]]

# PART ONE: UMAP

custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
umap_out <-
  umap(t(log2(assayDataElement(target_dsp_data, elt = "Q3_Norm"))),  
       config = custom_umap)
pData(target_dsp_data)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

# Figure 05a: UMAP 1
fig05a <- ggplot(pData(target_dsp_data),
                 aes(x = UMAP1, y = UMAP2,
                     color = PatientID, shape = Type)) +
  geom_point(size = 3) +
  theme_bw()
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/05a.png", 
    height=2400, 
    width=3200, 
    res=300)
fig05a
dev.off()

# Figure 05b: UMAP 2
fig05b <- ggplot(pData(target_dsp_data),
                 aes(x = UMAP1, y = UMAP2,
                     color = PatientID, shape = Subtype)) +
  geom_point(size = 3) +
  theme_bw()
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/05b.png", 
    height=2400, 
    width=3200, 
    res=300)
fig05b
dev.off()


# PART TWO: TSNE

set.seed(42)
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_dsp_data, elt = "Q3_Norm"))),
        perplexity = ncol(target_dsp_data)*.15)
pData(target_dsp_data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]

# Figure 05c: tSNE 1
fig05c <- ggplot(pData(target_dsp_data),
                 aes(x = tSNE1, y = tSNE2,
                     color = PatientID, shape = Type)) +
  geom_point(size = 3) +
  theme_bw()
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/05c.png", 
    height=2400, 
    width=3200, 
    res=300)
fig05c
dev.off()

# Figure 05d: tSNE 2
fig05d <- ggplot(pData(target_dsp_data),
                 aes(x = tSNE1, y = tSNE2,
                     color = PatientID, shape = Subtype)) +
  geom_point(size = 3) +
  theme_bw()
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/05d.png", 
    height=2400, 
    width=3200, 
    res=300)
fig05d
dev.off()

# PART THREE: UNSUPERVISED CLUSTERING

# Create a log2 transform of the data for analysis
assayDataElement(object = target_dsp_data, elt = "log_2") <-
  assayDataApply(target_dsp_data, 2, FUN = log, base = 2, elt = "Q3_Norm")

# Create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_dsp_data,
                         elt = "log_2", MARGIN = 1, calc_CV)
# Show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>      ZPBP   ANGPTL4    TTC21A      PLP1    CARNS1 
#> 0.5310768 0.4384522 0.4135341 0.4132069 0.4073844

# Identify genes in the top 20 percent of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
fig05e <- pheatmap(assayDataElement(target_dsp_data[GOI, ], elt = "log_2"),
                   scale = "row", 
                   show_rownames = FALSE, show_colnames = FALSE,
                   border_color = NA,
                   clustering_method = "average",
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   breaks = seq(-3, 3, 0.05),
                   color = colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))(120),
                   annotation_col = 
                     pData(target_dsp_data)[, c("PatientID", "Type", "Verhaak et al. Subtype")])
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/05e.png", 
    height=2400, 
    width=3200, 
    res=300)
fig05e
dev.off()

# Distance Matrix
sampleDists <- stats::dist(t(target_dsp_data@assayData[["exprs"]]))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_row = pData(target_dsp_data)[, c("PatientID", "Type", "Verhaak et al. Subtype")],
         annotation_col = pData(target_dsp_data)[, c("PatientID", "Type", "Verhaak et al. Subtype")],
         col = colors)
poisd <- PoissonDistance(t(target_dsp_data@assayData[["exprs"]]))

# Multidimensional scaling
verhaak_genes <- c(verhaak_gs[["VERHAAK_GLIOBLASTOMA_CLASSICAL"]]@geneIds,
                   verhaak_gs[["VERHAAK_GLIOBLASTOMA_MESENCHYMAL"]]@geneIds,
                   verhaak_gs[["VERHAAK_GLIOBLASTOMA_NEURAL"]]@geneIds,
                   verhaak_gs[["VERHAAK_GLIOBLASTOMA_PRONEURAL"]]@geneIds)
mds_goi <- verhaak_genes[verhaak_genes %in% GOI]
mds_dsp_data <- target_dsp_data[mds_goi, ]
mds_dist <- stats::dist(t(mds_dsp_data@assayData[["log_2"]]))
mds <- as.data.frame(mds_dsp_data@phenoData@data) |>
  cbind(cmdscale(mds_dist))
ggplot(mds, aes(x = `1`, y = `2`, color = Subtype)) +
  geom_point(size = 3) + coord_fixed()


mds_dist <- poisd$dd
mds <- as.data.frame(mds_dsp_data@phenoData@data) |>
  cbind(cmdscale(mds_dist))
ggplot(mds, aes(x = `1`, y = `2`, color = Subtype)) +
  geom_point(size = 3) + coord_fixed()


