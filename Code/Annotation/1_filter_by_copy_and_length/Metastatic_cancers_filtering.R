# METASTATIC CANCERS

library(ggplot2)
library(dplyr)
library(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Get data for MET500 and METAPRISM together in one table
met_df <- read.csv("../../../Data/Raw_data/data/MET500/CNV_table_metaprism_met500.txt", sep="\t")
MET_500 <- met_df[met_df$dataset == 'MET500',]
METAPRISM <- read.csv('../../../Data/Raw_data/data/METAPRISM/Data_Table_5.somatic_cna_segments_METAPRISM.tsv', sep="\t")
METAPRISM$dataset <- 'METAPRISM'
MET_500 <- MET_500 %>%
  separate(sample, into = c("sample.id.tumor", "sample.id.normal"), sep = "_Vs_")
MET_500$svtype <- NA
MET_500$svlen <- NA
MET_500$copy.number.simple <- NA
MET_500$copy.number.category <- NA
met_df <- rbind(MET_500, METAPRISM)

# Thresholding by length and copy number

# Visualize regions count by length and copy number

met_df$length <- met_df$end - met_df$start


ggplot(met_df, aes(x=length)) + 
  geom_histogram() +
  ggtitle("Distributions of length of regions")

ggplot(met_df, aes(x=log2(length))) + 
  geom_histogram() +
  ggtitle("Distributions of length of regions")


ggplot(met_df, aes(x=log10(tcn.em)))+
  geom_histogram() +
  ggtitle("Distributions of copy number of regions")

# Visualize sample count by length and copy number

summary_by_sample <- met_df %>%
  group_by(sample.id.tumor) %>%
  summarize(n_high_copy_regions = sum(tcn.em >= 15)) %>%
  ungroup()

ggplot(summary_by_sample, aes(x=n_high_copy_regions)) + 
  geom_histogram() +
  ggtitle("Distributions of samples by copy number") 

# This is without length filtration (but with filtration it is not different)
# How many samples have no regions of high copy number (< 15) --> 767 (71%)
length(summary_by_sample$sample.id.tumor[summary_by_sample$n_high_copy_regions == 0])

# How many samples have regions of high copy number (>= 15) --> 306 (29%)
length(summary_by_sample$sample.id.tumor[summary_by_sample$n_high_copy_regions != 0])



# Filter data (this leaves only 866 rows out of 84954)

met_df_filtered <- filter(met_df, tcn.em >= 15)
met_df_filtered_length <- filter(met_df_filtered, length >= 280000 & length <= 1900000)

write.csv(met_df_filtered_length, "Metastatic_cancers_filtered_by_copy_and_length.csv", row.names = F)
