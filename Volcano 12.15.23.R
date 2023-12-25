library(readxl)
library(viridisLite)
library(tidyverse)
library(EnhancedVolcano)

DGEA_type <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/Volcano.Data.Q3Norm.Filtered.xlsx") |>
  as.data.frame()
DGEA_vessels <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/Volcano.Data.SubsetCore.xlsx") |>
  as.data.frame()

V1 <- EnhancedVolcano(DGEA_type,
                title = "Edge versus Core",
                subtitle = "",
                caption = "",
                lab = DGEA_type$`Target name`,
                x = "Log2",
                y = "Pvalue",
                ylim = c(0, -log10(10e-6)),
                pCutoff = 0.05,
                pointSize = 3.0,
                labSize = 6,
                col=turbo(4),
                shape = c(1, 4, 23, 25),
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'top',
                legendLabSize = 16,
                legendIconSize = 5.0,
                max.overlaps = 15)

png("V1.png", height=2400, width=3200, res=300)
V1
dev.off()

V2 <- EnhancedVolcano(DGEA_vessels,
                            title = "Vessel-Low versus Vessel-High",
                            subtitle = "",
                            caption = "",
                            lab = DGEA_vessels$`Target name`,
                      x = "Log2",
                      y = "Pvalue",
                      ylim = c(0, -log10(10e-6)),
                      pCutoff = 0.05,
                      pointSize = 3.0,
                      labSize = 6,
                      col=viridis(2),
                      shape = c(1),
                      drawConnectors = TRUE,
                      boxedLabels = TRUE,
                      legendPosition = 'top',
                      legendLabSize = 16,
                      legendIconSize = 5.0,
                      max.overlaps = 15)

png("V2.png", height=2400, width=3200, res=300)
V2
dev.off()
