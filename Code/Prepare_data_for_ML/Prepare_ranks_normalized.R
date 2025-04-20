library(tidyr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# PREPARE DATA FOR MACHINE LEARNING
# IMPORT DATA SET
all_data <- read.csv('Table_with_rank_04_03.csv')

# CHECK HOW MANY ONCOGENES PER SAMPLE IS PRESENT 
all_data_summary <- all_data %>%
  group_by(sample) %>%
  summarise(
    n_oncogenes = sum(Role.in.Cancer.COSMIC %in% c("oncogene", "oncogene, fusion") | Is.Oncogene.oncoKB == "Yes")
  )

# ADD COLUMN IS ONCOGENE
all_data$is_oncogene <- ifelse(all_data$Is.Oncogene.oncoKB == "Yes" | 
                                 all_data$Role.in.Cancer.COSMIC %in% 
                                 c("oncogene, fusion", "oncogene"), 1, 0)


# NOW WE NORMALIZE RANKS WITHIN SAMPLE
all_data <- all_data[, colnames(all_data) %in% c('sample', 'chr', 'startpos', 
                                                 'endpos', 'copy_number', 
                                                 'gene_name', 'Role.in.Cancer.COSMIC',
                                                 'Is.Oncogene.oncoKB', 'pbm',
                                                 'tf', 'kin', 'crisp1', 'crisp2',
                                                 'GO', 'PPI', 'expr_rank',
                                                 'expr_fold_rank', 'rank', 'is_oncogene')]

colnames(all_data)[9:18] <- c('pubmed_mean_rank', 'tf_rank', 'kinase_rank', 'crisp_mean_rank',
                              'crisp_min_rank', 'GO_rank', 'PPI_rank', 'expr_rank',
                              'expr_fold_rank', 'sum_of_ranks')


norm_ranked_table <- all_data %>%
  group_by(sample) %>%
  mutate(
    pubmed_mean_rank_norm = round((pubmed_mean_rank / n()), 2),
    tf_rank_norm = round((tf_rank / n()), 2),
    kinase_rank_norm = round((kinase_rank / n()), 2),
    crisp_mean_rank_norm = round((crisp_mean_rank / n()), 2),
    crisp_min_rank_norm = round((crisp_min_rank / n()), 2),
    GO_rank_norm = round((GO_rank / n()), 2),
    PPI_rank_norm = round((PPI_rank / n()), 2),
    expr_rank_norm = round((expr_rank / n()), 2),
    expr_fold_rank_norm = round((expr_fold_rank / n()), 2),
    sum_of_normalized_ranks = (pubmed_mean_rank_norm + tf_rank_norm + kinase_rank_norm
                               + crisp_mean_rank_norm + crisp_min_rank_norm + GO_rank_norm +
                                 PPI_rank_norm + expr_rank_norm + expr_fold_rank_norm)
  ) %>%
  ungroup()



# PREPARE TRAINING DF: SELECT ONLY REGIONS HAVING ONCOGENES
train_df <- norm_ranked_table %>%
  group_by(sample) %>%
  filter(any(is_oncogene == 1)) %>%
  ungroup()

# REMOVE COLUMNS THAT WE WONT NEED IN ML
train_df <- train_df[, c(c(1:8), c(20:28), 29, 19)]

# PREPARE DF FOR PREDICTIONS: SELECT ONLY REGIONS WITHOUT KNOWN ONCOGENES
prediction_df <- norm_ranked_table %>%
  group_by(sample) %>%
  filter(!any(is_oncogene != 0)) %>%
  ungroup()

# REMOVE COLUMNS THAT WE WONT NEED IN ML
prediction_df <- prediction_df[, c(c(1:8), c(20:28), 29)]

# ADD SAMPLE NUMBER TO MAKE IT POSSIBLE TO MAKE PREDICTION SAMPLE BY SAMPLE
prediction_df <- prediction_df %>%
  group_by(sample) %>%
  mutate(sample_num = cur_group_id()) %>%
  ungroup()


# SAVE DATA 
write.csv(train_df, "Train_df_TCGA_HARTWIG_METASTATIC.csv", row.names = F)
write.csv(prediction_df, "Prediction_df_TCGA_HARTWIG_METASTATIC.csv", row.names = F)

# N samples in train df 2100, in prediction df  853
length(unique(train_df$sample))
length(unique(prediction_df$sample))
