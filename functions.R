# functions for data analysis

list_concordant_genes <- function(geomx_data, other_data){
  # compute median counts
  geomx_data_subset_counts <- geomx_data[-1]
  geomx_data_median_counts <- data.frame(geomx_data[1], MedianCounts=apply(geomx_data_subset_counts, 1, median, na.rm=TRUE))
  other_data_subset_counts <- other_data[-1]
  other_data_median_counts <- data.frame(other_data[1], MedianCounts=apply(other_data_subset_counts, 1, median, na.rm=TRUE))
  
  # reorder targets (descending)
  geomx_data_median_counts <- arrange(geomx_data_median_counts, desc(geomx_data_median_counts$MedianCounts))
  other_data_median_counts <- arrange(other_data_median_counts, desc(other_data_median_counts$MedianCounts))
  
  # match targets (remove unique targets)
  geomx_data_median_counts <- geomx_data_median_counts %>%
    filter(geomx_data_median_counts$TargetName %in% other_data_median_counts$TargetName)
  other_data_median_counts <- other_data_median_counts %>%
    filter(other_data_median_counts$TargetName %in% geomx_data_median_counts$TargetName)
  
  # restrict targets to top 100
  geomx_data_median_counts <- geomx_data_median_counts[1:1000,]
  other_data_median_counts <- other_data_median_counts[1:1000,]
  
  # calculate percent concordance
  concordant_genes <- as.data.frame(geomx_data_median_counts$TargetName) %>%
    filter(geomx_data_median_counts$TargetName %in% other_data_median_counts$TargetName)
  names(concordant_genes)[names(concordant_genes) == "geomx_data_median_counts$TargetName"] <- "TargetName"
  
  return(concordant_genes)
}