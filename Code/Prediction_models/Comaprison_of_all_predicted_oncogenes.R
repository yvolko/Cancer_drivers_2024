library(tidyr)
library(dplyr)
library(VennDiagram)

# From all the prediction genes labeled as TSG in COSMIC and oncoKB are removed

# RANKS
ranks <- read.csv("../../Data/Result_data/Oncogenes_from_ranks.txt", header = F)
oncogenes_from_ranks <- unique(ranks$V1)

# RANDOM FOREST CLASSIFIER
prediction_prob_RFclf <- read.csv("../../Data/Result_data/Predictions_by_RandomForestClassifier.csv")
prediction_prob_RFclf <- prediction_prob_RFclf[, c(2, 3, 16)] 

# Select only genes that are having probability of >= 0.89 
prediction_prob_RFclf <- prediction_prob_RFclf[prediction_prob_RFclf$prob_of_onco >= 0.89, ]

# Count how many time each gene is present 
summary_prob_pred <- prediction_prob_RFclf %>%
  group_by(gene_name) %>%
  summarise(count = n()) %>%
  ungroup()

oncogenes_from_RF <- unique(prediction_prob_RFclf$gene_name)
oncogenes_from_RF <- oncogenes_from_RF[oncogenes_from_RF != ""]

# Also count how many time gene was predicted as oncogenes in different samples
RFclf_samples <- prediction_prob_RFclf %>%
  separate(ID, into = c("sample", "rest"), sep = "_", extra = "drop") %>%
  select(sample) %>%
  rename(column = sample)

prediction_prob_RFclf$sample <- RFclf_samples

summary_per_sample_2 <- prediction_prob_RFclf %>%
  group_by(gene_name) %>%
  summarise(count = n_distinct(sample)) %>%
  ungroup()

# writeLines(oncogenes_from_RF, '../../Data/Result_data/Oncogenes_from_RF.txt')

# ENSEMBLE
ensemble_predictions <- read.csv('../../Data/Result_data/Prediction_by_ensemble.csv')

ens_summary <- ensemble_predictions %>%
  group_by(gene_name) %>%
  summarise(count = n()) %>%
  ungroup()

# Extract gene names

oncogenes_from_ensemble <- unique(ensemble_predictions$gene_name)
# writeLines(oncogenes_from_ensemble, '../../Data/Result_data/Oncogenes_from_ensemble.txt')

# COMAPARE SETS OF PUTATIVE ONCOGENES FROM DIFFERENT TECHNIQUES
oncogenes_from_ranks
oncogenes_from_RF
oncogenes_from_ensemble

VennDiagramm <- venn.diagram(
  x = list(oncogenes_from_ranks, oncogenes_from_RF, oncogenes_from_ensemble),
  category.names = c("Ranks" , "RFClassifier", "Ensemble"),
  filename = NULL,
  output=TRUE,
  fill = c('#7ac4cd', '#fde725ff', "#fed9a5"),
  cat.cex = 0.7,
  cat.fontsize = 20,   
  cex.line = 2  
)
grid.draw(VennDiagramm)


# Find overlap of all 3
overlap <- intersect(oncogenes_from_ranks, oncogenes_from_RF)

# Save this list of oncognes
# writeLines(overlap, '../../Data/Result_data/Oncogenes_predicted_by_all_3_models.txt')

# GET TABLE FOR RANKS WITH COUNTS
ranks_with_counts <- read.csv('../../Data/Result_data/Ranks_no_TSG.csv')

# SELECT ONLY RELEVANT GENES AND COLUMNS
overlap_of_ranks <- ranks_with_counts[ranks_with_counts$gene_name %in% overlap, c(2,4)]
overlap_of_RF <- summary_prob_pred[summary_prob_pred$gene_name %in% overlap, ]
overlap_of_Ens <- ens_summary[ens_summary$gene_name %in% overlap, ]

colnames(overlap_of_ranks)[2] <- "count_Ranks"
colnames(overlap_of_RF)[2] <- "count_RF"
colnames(overlap_of_Ens)[2] <- "count_Ens"

overlap_counts <- left_join(overlap_of_ranks, overlap_of_RF)
overlap_counts <- left_join(overlap_counts, overlap_of_Ens)

overlap_counts$total_counts <- with(overlap_counts, count_Ranks + count_RF + count_Ens)


# Add data on how often genes were present in the inintial dataset
initial_intervals <- read.csv('../../Data/Processed_data/Selection_of_regions_for_ML_without_known_oncogenes.csv')
initial_intervals <- initial_intervals[initial_intervals$gene_name != "",]

count_of_genes_in_initial_int <- initial_intervals %>%
  group_by(gene_name) %>%
  summarise(frequency_in_df = n()) %>%
  ungroup()
  
overlap_counts_with_gene_freq <- left_join(overlap_counts, count_of_genes_in_initial_int)
overlap_counts_with_gene_freq$count_Ranks_frequency_in_df <- 
  round(overlap_counts_with_gene_freq$count_Ranks / overlap_counts_with_gene_freq$frequency_in_df, 2)

overlap_counts_with_gene_freq$count_RF_frequency_in_df <- 
  round(overlap_counts_with_gene_freq$count_RF / overlap_counts_with_gene_freq$frequency_in_df, 2)

overlap_counts_with_gene_freq$count_Ens_frequency_in_df <- 
  round(overlap_counts_with_gene_freq$count_Ens / overlap_counts_with_gene_freq$frequency_in_df, 2)

write.csv(overlap_counts_with_gene_freq, '../../Data/Result_data/Overlap_counts_with_gene_freq.csv')

