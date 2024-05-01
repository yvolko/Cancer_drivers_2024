library(tidyr)
library(dplyr)
library(ggplot2)


# PREPARE REGIONS WITHOUT KNOWN ONCOGENES FOR PREDICTION

all_data <- read.csv("../result_files/df_with_GO_and_PPI.csv", sep="\t")

all_data$is_oncogene <- ifelse(all_data$Is.Oncogene.oncoKB == "Yes" | 
                                        all_data$Role.in.Cancer.COSMIC %in% 
                                 c("oncogene, TSG, fusion", "oncogene, fusion", 
                                   "oncogene, TSG", "oncogene"), 1, 0)



data <- all_data %>% 
  group_by(chr, startpos, endpos) %>%
  dplyr::filter(!any(is_oncogene == 1)) %>%
  ungroup()

# REMOVE ALL TSG FROM DF

data <- data[data$Is.Tumor.Suppressor.Gene.oncoKB != "Yes", ]
data <- data[(data$Role.in.Cancer.COSMIC != "TSG"), ]
data <- data[(data$Role.in.Cancer.COSMIC != "TSG, fusion"), ]

# REMOVE UNNECESSARY COLUMNS

data$X <- NULL
data$nMinor <- NULL
data$CNclass <- NULL
data$CNsignatureMapping <- NULL
data$ENS_id <- NULL
data$sample_type <- NULL
data$Tier.COSMIC <- NULL
data$Hallmark.COSMIC <- NULL
data$Role.in.Cancer.COSMIC <- NULL
data$Primary_disease <- NULL
data$Mutation.Types.COSMIC <- NULL
data$Translocation.Partner.COSMIC <- NULL
data$Is.Oncogene.oncoKB <- NULL
data$Is.Tumor.Suppressor.Gene.oncoKB <- NULL
data$OncoKB.Annotated.oncoKB <- NULL

# MAKE FOLD CHANGE OF EXPRESSION. Current units are log2(TPM + 0.001)
data$expr_fold_change <- round((2^data$expr)/(2^data$norm_expr), digits = 2)
data$norm_expr <- NULL
data$expr <- NULL


# TURN FEATURE INTO RANKS
data$is_kinase=ifelse(data$Family_kinase != "No_data","Yes","No")
data$Group_kinase <- NULL
data$Family_kinase <- NULL


data$pubmed_sum <- rowSums(data[, c("pubmed_cancer", "pubmed_growth", "pubmed_prolifiration", "pubmed_migration", "pubmed_invasion")])
data$pubmed_cancer <- NULL
data$pubmed_growth <- NULL
data$pubmed_prolifiration <- NULL
data$pubmed_migration <- NULL
data$pubmed_invasion <- NULL


data$Is.TF. <- ifelse(data$Is.TF. != "No_data","Yes","No")
data$DBD <- NULL


ranked_table <- data %>%
  group_by(sample, chr, startpos, endpos) %>%
  mutate(
    rank_expr_fold_change = round((rank(expr_fold_change) / n()), 2),
    rank_pubmed_sum = round((rank(pubmed_sum) / n()), 2),
    rank_TF = round((rank(Is.TF.) / n()), 2),
    rank_median.CRISPR = round((rank(desc(median.CRISPR)) / n()), 2),
    rank_mean.CRISPR = round((rank(desc(mean.CRISPR)) / n()), 2),
    rank_min.CRISPR = round((rank(desc(min.CRISPR)) / n()), 2),
    rank_kinase = round((rank(is_kinase) / n()), 2),
    rank_pubmed_mean = round((rank(pubmed_mean) / n()), 2),
    rank_GO_terms = round((rank(GO_terms) / n()), 2),
    rank_PPI = round((rank(PPI_count) / n()), 2)
  ) %>%
  ungroup()


# REMOVE UNNECESSARY COLUMNS

ranked_table$Is.TF. <- NULL
ranked_table$median.CRISPR <- NULL
ranked_table$mean.CRISPR <- NULL
ranked_table$min.CRISPR <- NULL
ranked_table$sd.CRISPR <- NULL
ranked_table$expr_fold_change <- NULL
ranked_table$is_kinase <- NULL
ranked_table$pubmed_sum <- NULL
ranked_table$PPI_count <- NULL
ranked_table$GO_terms <- NULL
ranked_table$pubmed_mean <- NULL

ranked_table$is_oncogene <- NULL

# Removed as it is the same as pubmed mean
ranked_table$rank_pubmed_sum <- NULL

ranked_table$SUM_RANK <- rowSums(ranked_table[, c(8:16)])

# GIVE ID TO EACH REGION TO ALLOW PREDICTION BY REGION
ranked_table <- ranked_table %>%
  group_by(sample, chr, startpos, endpos) %>%
  mutate(group_num = cur_group_id()) %>%
  ungroup()

write.csv(ranked_table, "../Dats/Processed_data/Selection_of_regions_without_known_oncognes.csv")