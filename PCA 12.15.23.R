library("scatterplot3d")
library(tidyverse)
library(readxl)
library(viridisLite)

# Read data
PCA_PoV <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/NormQ3.PCA.xlsx", range = "A1:B6") |>
  as.data.frame()
PCA_data <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/NormQ3.PCA.xlsx", range = "A9:J53") |>
  as.data.frame()
annot <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/TMA Annotation File (EDIT).xlsx") |>
  as.data.frame()

# Layer PatientID and Type Factors onto PCA_data
annot <- annot[-c(19, 46, 47), ]
PCA_data <- cbind(PCA_data, as.factor(annot$PatientID), as.factor(annot$Type))

# Layer Collated PatientID/Type Factor onto PCA_data
PCA_data$collate <- paste(PCA_data$`as.factor(annot$PatientID)`, PCA_data$`as.factor(annot$Type)`, sep = "+") |>
  as.factor()

# Create Plots
png("PCA1.png", height=1766, width=2000, res=300)
shapes = c(15, 17)
shapes <- shapes[PCA_data$`as.factor(annot$Type)`]
colors = viridis(2)
colors <- colors[PCA_data$`as.factor(annot$Type)`]
plot1 <- scatterplot3d(PCA_data[, 6:8], pch = shapes, color = colors)
dev.off()

png("PCA2.png", height=1766, width=2000, res=300)
shapes = c(15:20)
shapes <- shapes[PCA_data$`as.factor(annot$PatientID)`]
colors <- plasma(5)
colors <- colors[PCA_data$`as.factor(annot$PatientID)`]
plot2 <- scatterplot3d(PCA_data[, 6:8], pch = shapes, color = colors)
dev.off()

png("PCA3.png", height=1766, width=2000, res=300)
shapes = c(0:9)
shapes <- shapes[PCA_data$collate]
colors <- turbo(10)
colors <- colors[PCA_data$collate]
plot3 <- scatterplot3d(PCA_data[, 6:8], pch = shapes, color = colors, angle = 55)
dev.off()
