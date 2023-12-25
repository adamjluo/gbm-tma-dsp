library(tidyverse)

plot.data <- data.frame(Condition = c("B7-H3 Wildtype", "Scramble", "KO-Unsorted", "KO-Sorted"), Intensity = c(1, 1.013102872, 0.156079351, 0.000570524))

ggplot(plot.data, aes(x = Condition, y = Intensity)) +
  geom_bar(stat = "identity")
