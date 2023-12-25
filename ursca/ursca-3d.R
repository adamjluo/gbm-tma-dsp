# Bar Plot for Fig. 3d

library(readxl)
library(extrafont)
loadfonts()
library(viridisLite)
library(tidyverse)

fig_3d_data <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/ursca/ursca-3d.xlsx", range = "D7:AU14") |>
  as.data.frame()

fig_3d_data <- fig_3d_data[c(1, 3, 7),]
fig_3d_data <- t(fig_3d_data)
colnames(fig_3d_data) <- c("ROI", "Type", "Count")
fig_3d_data <- data.frame(fig_3d_data)

fig_3d <- ggplot(fig_3d_data, aes(x = as.numeric(ROI), y = as.numeric(Count), fill = Type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  xlab("ROI") + ylab("CD276 (Count)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 16, family = "serif"))

png("fig3d.png", height=2400, width=3200, res=300)
fig_3d
dev.off()

