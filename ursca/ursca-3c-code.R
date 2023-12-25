# Volcano Plot for Fig. 3c

library(readxl)
library(viridisLite)
library(tidyverse)
library(EnhancedVolcano)

fig_3c_data <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/ursca/ursca-3c.xlsx", range = "C7:H8537") |>
  as.data.frame()

fig_3c <- EnhancedVolcano(fig_3c_data,
                          title = "Tumor Periphery (Left) versus Core (Right)",
                          subtitle = "",
                          caption = "",
                          lab = fig_3c_data$`Target name`,
                          selectLab = c("CD276", "DAAM2", "NTRK2", "FABP7", "LAMC1"),
                          x = "Log2",
                          y = "Pvalue",
                          ylim = c(0, -log10(10e-6)),
                          pCutoff = 0.05,
                          pointSize = 3.0,
                          labSize = 6.0,
                          col=plasma(4),
                          colAlpha = 1,
                          drawConnectors = TRUE,
                          boxedLabels = TRUE,
                          legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                         'p-value & Log (base 2) FC'),
                          legendPosition = 'right',
                          legendLabSize = 16,
                          legendIconSize = 5.0
                          )

png("fig3c.png", height=2400, width=3200, res=300)
fig_3c
dev.off()