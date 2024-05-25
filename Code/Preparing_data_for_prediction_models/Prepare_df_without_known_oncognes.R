library(tidyr)
library(dplyr)
library(ggplot2)


# PREPARE REGIONS WITHOUT KNOWN ONCOGENES FOR PREDICTION

all_data <- read.csv('../../Data/Processed_data/Table_with_rank_filtering_by_oncogene.csv')

all_data$is_oncogene <- ifelse(all_data$Is.Oncogene.oncoKB == "Yes" | 
                                 all_data$Role.in.Cancer.COSMIC %in% 
                                 c("oncogene, fusion", "oncogene"), 1, 0)

# SELECT ONLY REGIONS NOT HAVING KNOWN ONCOGENES
data <- all_data %>%
  group_by(chr) %>%
  filter(all(is_oncogene == 0)) %>%
  ungroup()

# REMOVE COLUMNS THAT WE WONT NEED IN ML
data <- data[, -c(1, 3, 4, 5, c(7:21))]

# RENAME SOME COLUMNS
colnames(data) <- c("ID", "gene_name", "rank_expr", "rank_pubmed_mean", "rank_tf", "rank_kinase",
                    "rank_crisp_mean", "rank_crisp_min", "rank_expr_fold_change", "rank_is_TSG",
                    "rank_GO_term", "rank_PPI", "sum_of_ranks", "is_oncogene")

# WE ALREADY HAVE RANKS
# NOW WE NORMALIZE RANKS WITHIN REGION
norm_ranked_table <- data %>%
  group_by(ID) %>%
  mutate(
    rank_expr_norm = round((rank_expr / n()), 2),
    rank_pubmed_mean_norm = round((rank_pubmed_mean / n()), 2),
    rank_tf_norm = round((rank_tf / n()), 2),
    rank_kinase_norm = round((rank_kinase / n()), 2),
    rank_crisp_mean_norm = round((rank_crisp_mean / n()), 2),
    rank_crisp_min_norm = round((rank_crisp_min / n()), 2),
    rank_expr_fold_change_norm = round((rank_expr_fold_change / n()), 2),  
    rank_is_TSG_norm = round((rank_is_TSG / n()), 2), 
    rank_GO_term_norm = round((rank_GO_term / n()), 2),
    rank_PPI_norm = round((rank_PPI / n()), 2)
  ) %>%
  ungroup()

# REMOVE UNNECESSARY COLUMNS
norm_ranked_table <- norm_ranked_table[, -c(3:14)]

# WE HAVE FEATURE SUM_OF_RANKS, HOWEVER TO MAKE IT ALSO NORMALIZED WITHIN THE REGIONL LET'S MAKE SUM OF NORMALIZED RANKS
norm_ranked_table$SUM_RANK <- rowSums(norm_ranked_table[, c(3:12)])

# GIVE NUMERICAL ID TO EACH REGION TO ALLOW PREDICTION BY REGION
norm_ranked_table <- norm_ranked_table %>%
  group_by(ID) %>%
  mutate(group_num = cur_group_id()) %>%
  ungroup()

# SAVE DATA 
write.csv(norm_ranked_table, "../../Data/Processed_data/Selection_of_regions_for_ML_without_known_oncogenes.csv")
