# INSTRUCTIONS

# Load relevant packages
packages <- c("ggrepel", "EnhancedVolcano", "viridisLite")
lapply(packages, library, character.only = TRUE)

# Differential Expression Analysis: Edge vs. Core
# LMM: log2(gene) ~ Type + (1 + Type | SpotID)
# LMM w/ random slope (comparing ROIs within the same tissue)

pData(target_dsp_data)$DEType <- 
  factor(pData(target_dsp_data)$Type)
pData(target_dsp_data)$DESpotID <- 
  factor(pData(target_dsp_data)$SpotID)

LMM_type <- function(object) {
  mixedOutmc <-
    mixedModelDE(object,
                 elt = "log_2",
                 modelFormula = ~ DEType + (1 + DEType | DESpotID),
                 groupVar = "DEType",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE,
                 pAdjust = "BH", 
                 pairwise = TRUE)
  dgea_type <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(dgea_type)
  dgea_type <- as.data.frame(dgea_type)
  dgea_type$Contrast <- tests
  dgea_type$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  dgea_type$FDR <- p.adjust(dgea_type$`Pr(>|t|)`, method = "BH")
  dgea_type <- dgea_type[, c("Gene", "Contrast", "Estimate", 
                             "Pr(>|t|)", "FDR")]
  row.names(dgea_type) <- NULL
  colnames(dgea_type) <- c("Gene", "Contrast", "FC", 
                           "P-value", "FDR")
  return(dgea_type)
}

all_patients <- target_dsp_data
dgea_result <- LMM_type(all_patients)

# Figure 06a:
fig06a <- EnhancedVolcano(dgea_result,
                          lab = dgea_result$Gene,
                          x = "FC",
                          y = "P-value",
                          title = "Linear Mixed-Effects Model, Edge vs. Core",
                          subtitle = "",
                          caption = "",
                          drawConnectors = TRUE,
                          pCutoff = 0.05,
                          FCcutoff = 1,
                          pointSize = 3,
                          labSize = 6, 
                          boxedLabels = TRUE, 
                          col = turbo(4),
                          max.overlaps = 40)
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/06a.png", 
    height=2400, 
    width=3200, 
    res=300)
fig06a
dev.off()


patient_one <- target_dsp_data[, pData(target_dsp_data)$PatientID == "JHU-113188965"]
dgea_result <- LMM_type(patient_one)

# Figure 06b:
fig06b <- EnhancedVolcano(dgea_result,
                          lab = dgea_result$Gene,
                          x = "FC",
                          y = "P-value",
                          title = "Linear Mixed-Effects Model, Edge vs. Core: Patient JHU-113188965",
                          subtitle = "",
                          caption = "",
                          drawConnectors = TRUE,
                          pCutoff = 0.05,
                          FCcutoff = 2,
                          pointSize = 3,
                          labSize = 6, 
                          boxedLabels = TRUE, 
                          col = turbo(4),
                          max.overlaps = 40)
png(file = "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/AACR 2024/Figures/06b.png", 
    height=2400, 
    width=3200, 
    res=300)
fig06b
dev.off()


