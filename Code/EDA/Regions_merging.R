library(dplyr)
library(ggplot2)
library(tidyverse)


my_df <- read.csv('../../Data/Processed_data/Table_with_gene_names_from_gencode.csv')

interval_df <- my_df %>%
  arrange(startpos) %>%
  group_by(sample, chr) %>%
  mutate(dif_interval = ifelse(lead(sample) == sample & 
                                 lead(startpos) != startpos &
                                 lead(chr) == chr, lead(startpos) - endpos, NA))

# DETERMINING THE INTERVALS BETWEEN REGIONS TO COMBINE

ggplot(interval_df, aes(x = log10(dif_interval), y = after_stat(density))) +
  geom_histogram(position = "dodge", bins = 25, color = "black", fill = "steelblue") +
  geom_density(color = "black", size = 1, alpha = 0.5) +
  labs(title = "Interval length distribution", x = "log10 Interval length", y = "Density") +
  theme_classic() 

# INTERVAL DETERMINED AS 20000 nb

for (i in unique(interval_df$sample)){
  for (j in unique(interval_df[interval_df$sample == i,]$chr)){
    test_df <- filter(interval_df, sample == i & chr == j) 
    if (nrow(test_df) > 1){
      for (k in 1:(nrow(test_df)-1)){
        if (test_df[k,]$startpos != test_df[k+1,]$startpos & 
            (test_df[k+1,]$startpos - test_df[k,]$endpos) <= 20000){
          interval_df$startpos[interval_df$startpos == test_df[k+1,]$startpos &
                                 interval_df$chr == test_df[k+1,]$chr &
                                 interval_df$sample == test_df[k+1,]$sample] =  test_df[k,]$startpos
          interval_df$endpos[interval_df$endpos == test_df[k,]$endpos & 
                               interval_df$sample == test_df[k,]$sample &
                               interval_df$chr == test_df[k,]$chr] =  test_df[k+1,]$endpos
        }
      }
    }
  }
}

new_interval_df <- interval_df %>%
  group_by(sample, chr, endpos) %>%
  mutate(startpos = min(startpos))

new_interval_df <- new_interval_df %>%
  group_by(sample, chr, endpos, startpos) %>%
  mutate(nMajor = round(mean(nMajor), digits = 0))

new_interval_df <- new_interval_df[,-11]
write.csv(new_interval_df, '../../Data/Processed_data/Merge_table_with_gene_names_from_gencode.csv', quote = FALSE, row.names = FALSE)
  