library(dplyr)
library(Metrics)


new_summary_data <- read.csv("../../Data/Processed_data/Table_with_rank_filtering_by_oncogene.csv")

summary_data_df <- new_summary_data %>%
  group_by(chr) %>%
  mutate(first_onco = ifelse(rank >= quantile(rank, probs = 0.9) & (Role.in.Cancer.COSMIC == "oncogene" | Is.Oncogene.oncoKB == 'Yes'), 1, 0))

summary_data_df_new <- summary_data_df %>%
  group_by(chr) %>%
  mutate(first_onco = ifelse(1 %in% first_onco, 1, 0),
         there_is_onco = ifelse(("oncogene" %in% Role.in.Cancer.COSMIC | 'Yes' %in% Is.Oncogene.oncoKB), 1, 0))

df_for_metrics <- summary_data_df_new[!duplicated(summary_data_df_new$chr), ]                            
df_for_metrics$first_onco <- result_vector


accuracy <- accuracy(df_for_metrics$there_is_onco, df_for_metrics$first_onco)
recall <- recall(df_for_metrics$there_is_onco, df_for_metrics$first_onco) 
