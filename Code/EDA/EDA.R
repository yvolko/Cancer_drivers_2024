#install.packages('readxl')
#install.packages('ggplot2')
#install.packages('dplyr')
#install.packages('tidyverse')


#Libraries
library (readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)


# IMPORT AND MERGE DATASETS
path = '../../Data/Raw_data/41586_2022_4738_MOESM3_ESM.xlsx'
barcode_annotation <- read_xlsx(path, sheet = 3)
df <- read_xlsx(path, sheet = 2)
cancer_type_annotation <- read_xlsx(path, sheet = 1)
df <- inner_join(df, barcode_annotation, by=c("sample"="name")) 
df <- inner_join(df, cancer_type_annotation, by = c('cancer_type' = 'Abbreviation'))

# SHOW DISTRIBUTIONS OF REGIONS WITH DIFFERENT COPY NUMBERS

nMajor_before_filt <- ggplot(data=df, aes(x=nMajor))+ 
  geom_histogram(fill = '#5F9EA0', color = 'white') +
  theme_classic() +
  scale_y_continuous(trans = 'log10') +
  labs(x = "Number of the region's copies", y = "Number of cases (log10)") +
  geom_vline(xintercept = 15, color = "red", linetype = "dashed", size = 1)

# SHOW DISTRIBUTIONS OF REGIONS WITH DIFFERENT LENGTHS

df$reg_length <- df$endpos - df$startpos
length_before_filt <- ggplot(data=df, aes(x=log2(reg_length)))+ 
  geom_histogram(fill = '#5F9EA0', color = 'white') +
  theme_classic() +
  labs(x = "Length (log2)", y = "Number of cases") +
  geom_vline(xintercept = log2(280000), color = "red", linetype = "dashed", size = 1)


# BASED ON LITERARY DATA FILTERING REGIONS WITH THE NUMBER OF COPIES LESS THAN 15

df_filt_by_nMajor <- filter(df, nMajor >= 15)

# LENGTH DISTRIBUTION DEPENDING ON THE nMajor's BIN

length_nMajor_summary <- aggregate(reg_length ~ nMajor, data = df_filt_by_nMajor, 
                                       FUN = function(x) c(mean = mean(x), sd = sd(x), median = median(x)))
length_nMajor_summary <- length_nMajor_summary %>% 
  mutate(bin = cut2(nMajor, g = 5))
ggplot(length_nMajor_summary, aes(x = log(reg_length[,'mean'], base = 2), fill = bin)) +
  geom_histogram(position = "dodge", binwidth = 0.8, color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Length distribution depending on the bin", x = "Mean length (log 2)", y = "Frequency") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# BASED ON PREVIOUS PLOT FILTERING REGIONS WITH LENGTH LESS THAN 280000 nb AND MORE THAN 1900000 nb

df_filt_by_nMajor_and_length <- filter(df_filt_by_nMajor, reg_length >= 280000 & reg_length <= 1900000)[,-c(11:23)]

write.csv(df_filt_by_nMajor_and_length, file = "../../Data/Raw_data/Table_filtering_by_nMajor_and_length.csv", row.names = FALSE)

# STATISTICS BY CANCER TYPES

sample_number_with_high_nMajor <- aggregate(sample ~ cancer_type, data = df_filt_by_nMajor_and_length, FUN = function(x) length(unique(x)))
sample_number_all <- aggregate(sample ~ cancer_type, data = df, FUN = function(x) length(unique(x)))
stats_cancer_type <- df_filt_by_nMajor_and_length %>%
  group_by(cancer_type) %>%
  summarise(mean = mean(nMajor), 
            sd = sd(nMajor),
            median = median(nMajor)) %>%
  inner_join(cancer_type_annotation, by = c('cancer_type' = 'Abbreviation')) %>%
  inner_join(sample_number_with_high_nMajor, by = 'cancer_type', suffix = c('', '_high_nMajor')) %>%
  inner_join(sample_number_all, by = 'cancer_type', suffix = c('', '_all')) %>%
  mutate(high_copy_frequency = round(sample / sample_all, digits = 2))
  
ggplot(stats_cancer_type, aes(x = reorder(cancer_type, -high_copy_frequency), y = high_copy_frequency)) +
  geom_bar(stat = "identity", position=position_dodge(), fill = '#5F9EA0') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  labs(x = "Cancer type", y = "Frequency")

# STATISTICS BY CHROMOSOME

sample_number_with_high_nMajor <- aggregate(sample ~ chr, data = df_filt_by_nMajor_and_length, FUN = function(x) length(unique(x)))
sample_number_all <- aggregate(sample ~ chr, data = df, FUN = function(x) length(unique(x)))
stats_chr <- df_filt_by_nMajor_and_length %>%
  group_by(chr) %>%
  summarise(mean = mean(nMajor), 
            sd = sd(nMajor),
            median = median(nMajor)) %>%
  inner_join(sample_number_with_high_nMajor, by = 'chr', suffix = c('', '_high_nMajor')) %>%
  inner_join(sample_number_all, by = 'chr', suffix = c('', '_all')) %>%
  mutate(high_copy_frequency = round(sample / sample_all, digits = 3))

ggplot(stats_chr, aes(x = reorder(chr, -high_copy_frequency), y = high_copy_frequency)) +
  geom_bar(stat = "identity", position=position_dodge(), fill = '#5F9EA0') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()+
  labs(x = "Chromosome", y = "Frequency")





