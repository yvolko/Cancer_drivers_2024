# METASTATIC CANCERS - INTERVAL MERGING

library(dplyr)
library(ggplot2)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df <- read.csv('../2_annotation_with_gene_names/Metastatic_cancers_with_gene_names_from_gencode.csv')

interval_df <- df %>%
  arrange(start) %>%
  group_by(sample.id.tumor, chrom) %>%
  mutate(dif_interval = ifelse(lead(sample.id.tumor) == sample.id.tumor & 
                                 lead(start) != start &
                                 lead(chrom) == chrom, lead(start) - end, NA))

# DETERMINING THE INTERVALS BETWEEN REGIONS TO COMBINE

ggplot(interval_df, aes(x = log10(dif_interval), y = after_stat(density))) +
  geom_histogram(position = "dodge", bins = 25, color = "black", fill = "steelblue") +
  geom_density(color = "black", size = 1, alpha = 0.5) +
  labs(title = "Interval length distribution", x = "log10 Interval length", y = "Density") +
  theme_classic() 

# INTERVAL DETERMINED AS 20000 nb (same as used for TCGA)

for (i in unique(interval_df$sample.id.tumor)){
  for (j in unique(interval_df[interval_df$sample.id.tumor == i,]$chrom)){
    test_df <- filter(interval_df, sample.id.tumor == i & chrom == j) 
    if (nrow(test_df) > 1){
      for (k in 1:(nrow(test_df)-1)){
        if (test_df[k,]$start != test_df[k+1,]$start & 
            (test_df[k+1,]$start - test_df[k,]$end) <= 20000){
          interval_df$start[interval_df$start == test_df[k+1,]$start &
                                 interval_df$chrom == test_df[k+1,]$chrom &
                                 interval_df$sample.id.tumor == test_df[k+1,]$sample.id.tumor] =  test_df[k,]$start
          interval_df$end[interval_df$end == test_df[k,]$end & 
                               interval_df$sample.id.tumor == test_df[k,]$sample.id.tumor &
                               interval_df$chrom == test_df[k,]$chrom] =  test_df[k+1,]$end
        }
      }
    }
  }
}

new_interval_df <- interval_df %>%
  group_by(sample.id.tumor, chrom, end) %>%
  mutate(start = min(start))

new_interval_df <- new_interval_df %>%
  group_by(sample.id.tumor, chrom, end, start) %>%
  mutate(tcn.em = round(mean(tcn.em), digits = 0))

new_interval_df$X <- NULL

write.csv(new_interval_df, 'Metastatic_cancers_with_gene_names_merged.csv', quote = FALSE, row.names = FALSE)
