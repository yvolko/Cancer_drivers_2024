library(tidyr)
library(dplyr)


# PREPARE DATA FOR MACHINE LEARNING
# IMPORT DATA SET

all_data <- read.csv("../../Data/Processed_data/df_with_GO_and_PPI.csv", sep="\t")

all_data$is_oncogene <- ifelse(all_data$Is.Oncogene.oncoKB == "Yes" | 
                                 all_data_GO_PPI$Role.in.Cancer.COSMIC %in% 
                                 c("oncogene, TSG, fusion", "oncogene, fusion", 
                                   "oncogene, TSG", "oncogene"), 1, 0)

# SELECT ONLY REGIONS HAVING ONCOGENES
data <- all_data %>%
  group_by(chr, startpos, endpos) %>%
  filter(any(is_oncogene == 1)) %>%
  ungroup()

# REMOVE COLUMNS THAT WE WONT NEED IN ML
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
data$oncogene_sum <- NULL

# MAKE FOLD CHANGE OF EXPRESSION. Current units are log2(TPM + 0.001)
data$expr_fold_change <- round((2^data$expr)/(2^data$norm_expr), digits = 2)

# SELECT ONLY ONE REGION PER ONCOGENE
# ADDITIONALLY REMOVE REGIONS WITH MORE THAN 1 ONCOGENE
data_unique <- data %>% 
  group_by(chr, startpos, endpos) %>%
  mutate(num_oncogenes = sum(is_oncogene)) %>%
  ungroup()

data <- data_unique %>% 
  group_by(chr, startpos, endpos) %>%
  filter(any(num_oncogenes == 1)) %>%
  ungroup()

data$num_oncogenes <- NULL

# Leave only one region per unique oncogene
oncogenes <- data[data$is_oncogene == 1, ]
oncogene_count <- table(oncogenes$gene_name)
oncogenes <- unique(oncogenes$gene_name)

result_list <- list()

# For all the oncogenes select only first region
for(oncogene in oncogenes) {
  # Filter data for the current gene
  gene_data <- data %>%
    group_by(chr, startpos, endpos) %>%
    filter(gene_name == oncogene) %>%
    ungroup()
  
  # Keep only the first row
  gene_first <- gene_data %>%
    slice(1)
  
  # Store the result in the list
  result_list[[oncogene]] <- gene_first
}

result_df <- do.call(rbind, result_list)
 
# Select respective regions from our data
result_df$coordinates <- paste(result_df$chr, result_df$startpos, result_df$endpos)
coordinates <- result_df$coordinates
data$coordinates <- paste(data$chr, data$startpos, data$endpos)

final_df <- data[data$coordinates %in% coordinates, ]
final_df$coordinates <- NULL


# TURN FEATURE INTO RANKS
# ADDITIONALLY NORMALIZE RANKS WITHIN REGIONS
final_df$is_kinase=ifelse(final_df$Family_kinase != "No_data","Yes","No")
final_df$Group_kinase <- NULL
final_df$Family_kinase <- NULL

final_df$pubmed_sum <- rowSums(final_df[, c("pubmed_cancer", "pubmed_growth", "pubmed_prolifiration", "pubmed_migration", "pubmed_invasion")])
final_df$pubmed_cancer <- NULL
final_df$pubmed_growth <- NULL
final_df$pubmed_prolifiration <- NULL
final_df$pubmed_migration <- NULL
final_df$pubmed_invasion <- NULL

final_df$Is.TF. <- ifelse(final_df$Is.TF. != "No_data","Yes","No")
final_df$DBD <- NULL

final_df$norm_expr <- NULL
final_df$expr <- NULL

ranked_table <- final_df %>%
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

# Removed as it is the same as pubmed mean
ranked_table$rank_pubmed_sum <- NULL

# MAKE ADDITIONAL FEATURE THAT IS SUM OF ALL RANKS
ranked_table$SUM_RANK <- rowSums(ranked_table[, c(9:17)])


# SAVE DATA 
write.csv(ranked_table, "../../Data/Processed_data/Selection_of_regions_for_ML.csv")
