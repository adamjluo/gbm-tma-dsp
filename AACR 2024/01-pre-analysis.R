# 01*: Standard pre-analysis workflow

# Load relevant packages

packages <- c("tidyverse", "readxl", "utils", "stats", 
              "scales", "BiocManager", "NanoStringNCTools", 
              "GeomxTools", "GeoMxWorkflows")
lapply(packages, library, character.only = TRUE)

# PART ONE: PRE-PROCESSING AND QUALITY CONTROL (QC)

# Load data into GeoMxSet object

datadir <- file.path("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBM_TMA_RawFiles")
DCCFiles <- dir(file.path(datadir, "dcc"),
                pattern = ".dcc$",
                full.names = TRUE,
                recursive = TRUE)

PKCFile <- dir(file.path(datadir, "pkc"),
                pattern = ".pkc$",
                full.names = TRUE)

AnnotationFile <- dir(file.path(datadir, "annotation"),
                      pattern = ".xlsx$",
                      full.names= TRUE)

dsp_data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                       pkcFiles = PKCFile,
                       phenoDataFile = AnnotationFile,
                       phenoDataSheet = "Template",
                       phenoDataDccColName = "Sample_ID",
                       protocolDataColNames = c("AOI", "ROI"),
                       experimentDataColNames = c("Panel")
                       )

# Store PKC file information in appropriate variables
pkcs <- annotation(dsp_data)
modules <- gsub(".pkc", "", pkcs)

# Shift zero count values to ones to enable downstream processing
dsp_data <- shiftCountsOne(dsp_data, useDALogic = TRUE)

# Generate segment QC flags using recommended parameters
SegQCParams <- list(minSegmentReads=1000, percentTrimmed=80, percentStitched=80, 
                 percentAligned=80, percentSaturation=50, minNegativeCount=10, 
                 maxNTCCount=1000, minNuclei=200, minArea=16000, minProbeCount=10, 
                 minProbeRatio=0.1, outlierTestAlpha=0.01, percentFailGrubbs=20, 
                 loqCutoff=2.0, highCountCutoff=10000)
dsp_data <- setSegmentQCFlags(dsp_data, qcCutoffs = SegQCParams)

# Consolidate and summarize segment QC results
QCResults <- protocolData(dsp_data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QCSummary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QCSummary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# Calculate and store negative probe geometric average data
negativeGeoMeans <- 
  esBy(negativeControlSubset(dsp_data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       })
protocolData(dsp_data)[["NegGeoMean"]] <- negativeGeoMeans

# Run probe QC using recommended parameters
dsp_data <- setBioProbeQCFlags(dsp_data,
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = FALSE)
ProbeQCResults <- fData(dsp_data)[["QCFlags"]]
ProbeQCPassed <- 
  subset(dsp_data, 
         fData(dsp_data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(dsp_data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dsp_data <- ProbeQCPassed

# Create gene-level count matrix for downstream processing
target_dsp_data <- aggregateCounts(dsp_data)
target_dsp_data <- target_dsp_data[, pData(target_dsp_data)$'Slide Name' != "No Template Control"]

# Calculate the limit of quantitation (LOQ)
cutoff <- 2
minLOQ <- 2
LOQ <- data.frame(row.names = colnames(target_dsp_data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_dsp_data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_dsp_data)[, vars[1]] *
             pData(target_dsp_data)[, vars[2]] ^ cutoff)
  }
}
pData(target_dsp_data)$LOQ <- LOQ


# PART 2: FILTERING

# Generate Boolean matrix of gene-level counts greater than LOQ per segment
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_dsp_data)$Module == module
  Mat_i <- t(esApply(target_dsp_data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# Ensure ordering of LOQ_Mat (stored outside target_dsp_data)
LOQ_Mat <- LOQ_Mat[fData(target_dsp_data)$TargetName, ]

# Save segment-level gene detection rate data to phenoData
pData(target_dsp_data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_dsp_data)$GeneDetectionRate <-
  pData(target_dsp_data)$GenesDetected / nrow(target_dsp_data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_dsp_data)$DetectionThreshold <- 
  cut(pData(target_dsp_data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.25, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", "15-25%", ">25%"))

# Fig. 01a: Gene detection rate by segment type (pre-filter)
fig01a <- ggplot(pData(target_dsp_data),
                 aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Type)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/01a.png", 
    height=2400, 
    width=3200, 
    res=300)
fig01a
dev.off()

# Filter segments 
target_dsp_data <-
  target_dsp_data[, pData(target_dsp_data)$GeneDetectionRate >= 0.1]

# Fig. 01b: Gene detection rate by segment type (post-filter)
fig01b <- ggplot(pData(target_dsp_data),
                 aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Type)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/01b.png", 
    height=2400, 
    width=3200, 
    res=300)
fig01b
dev.off()

# Generate feature-level segment detection rate data
LOQ_Mat <- LOQ_Mat[, colnames(target_dsp_data)]
fData(target_dsp_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_dsp_data)$SegmentDetectionRate <-
  fData(target_dsp_data)$DetectedSegments / nrow(pData(target_dsp_data))

# Fig. 01c: Segment detection rate by gene (pre-filter)
plot_df <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_df$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_dsp_data)$SegmentDetectionRate >= x)}))
plot_df$Rate <- plot_df$Number / nrow(fData(target_dsp_data))
rownames(plot_df) <- plot_df$Freq
fig01c <- ggplot(plot_df, 
                 aes(x = as.factor(Freq), 
                     y = Rate, 
                     fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "% of WTA Panel > LOQ, # Genes Detected")
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/01c.png", 
    height=2400, 
    width=3200, 
    res=300)
fig01c
dev.off()

# Perform target filter operation (making sure to include NegProbe-WTX)
negativeProbefData <- subset(fData(target_dsp_data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_dsp_data <- 
  target_dsp_data[fData(target_dsp_data)$SegmentDetectionRate >= 0.1 |
                    fData(target_dsp_data)$TargetName %in% neg_probes, ]

# Figure 01d: Q3 value & negative probe geometric average by segment type (core vs. edge)
# A: Histogram
# B: Scatter Plot
library(reshape2)
library(cowplot)
ann_of_interest <- "Type"
Stat_data <-
  data.frame(row.names = colnames(exprs(target_dsp_data)),
             Segment = colnames(exprs(target_dsp_data)),
             Annotation = pData(target_dsp_data)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_dsp_data), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_dsp_data)[neg_probes, ])
Stat_data_m <- melt(Stat_data,
                    measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic",
                    value.name = "Value")
plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) +
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")
plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43, 0.57))
fig01d <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/01d.png", 
    height=2400, 
    width=3200, 
    res=300)
fig01d
dev.off()




