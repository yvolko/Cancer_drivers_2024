library(tidyr)
library(dplyr)
library(VennDiagram)

# From all the prediction genes labeled as TSG in COSMIC and oncoKB are removed

# RANKS
ranks <- read.csv("../../Data/Result_data/Ranks_no_TSG.csv", sep=" ")

oncogenes_from_ranks <- unique(ranks$gene_name)
#writeLines(oncogenes_from_ranks, '../../Data/Result_data/Oncogenes_from_ranks.txt')

# RANDOM FOREST CLASSIFIER WITH PROBABILITIES
prediction_prob_df <- read.csv("../../Data/Result_data/Predictions_by_RandomForestClassifier.csv")
prediction_prob_df$Unnamed..0 <- NULL  

# Select only genes that are having probability of >= 0.9 
prediction_prob_df <- prediction_prob_df[prediction_prob_df$prob_of_onco >= 0.9, ]

# Count how many time each gene is present 
summary_prob_pred <- prediction_prob_df %>%
  group_by(gene_name) %>%
  summarise(count = n()) %>%
  ungroup()

oncogenes_from_RF <- unique(prediction_prob_df$gene_name)
oncogenes_from_RF <- oncogenes_from_RF[oncogenes_from_RF != ""]

# Also count how many time gene was predicted as ocogenes in different samples
summary_per_sample <- prediction_prob_df %>%
  group_by(sample)  %>%
  summarise(count = n(),
            gene_names = paste(unique(gene_name), collapse = ", ")) %>%
  ungroup()
  
summary_per_sample_2 <- prediction_prob_df %>%
  group_by(gene_name) %>%
  summarise(count = n_distinct(sample)) %>%
  ungroup()

#writeLines(oncogenes_from_RF, '../../Data/Result_data/Oncogenes_from_RF.txt')

# ENSEMBLE
ensemble_predictions <- read.csv('../../Data/Result_data/Prediction_by_ensemble.csv')

ens_summary <- ensemble_predictions %>%
  group_by(gene_name) %>%
  summarise(count = n()) %>%
  ungroup()

# Extract gene names

oncogenes_from_ensemble <- unique(ensemble_predictions$gene_name)
#writeLines(oncogenes_from_ensemble, '../../Data/Result_data/Oncogenes_from_ensemble.txt')

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
#writeLines(overlap, 'Oncogenes_predicted_by_all_3_models.txt')

overlap_of_ranks <- ranks[ranks$gene_name %in% overlap, c(1,3)]
overlap_of_RF <- summary_prob_pred[summary_prob_pred$gene_name %in% overlap, ]
overlap_of_Ens <- ens_summary[ens_summary$gene_name %in% overlap, ]

colnames(overlap_of_ranks)[2] <- "count_Ranks"
colnames(overlap_of_RF)[2] <- "count_RF"
colnames(overlap_of_Ens)[2] <- "count_Ens"

overlap_counts <- left_join(overlap_of_ranks, overlap_of_RF)
overlap_counts <- left_join(overlap_counts, overlap_of_Ens)

overlap_counts$total_counts <- with(overlap_counts, count_Ranks + count_RF + count_Ens)

# Sort by counts
overlap_counts <- overlap_counts[order(overlap_counts$total_counts), ]

# Reset row indexes
rownames(overlap_counts) <- NULL

# Add data on how often genes were present in the inintial dataset
initial_intervals <- read.csv('../../Data/Processed_data/Selection_of_regions_without_known_oncognes.csv')
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

# write.csv(overlap_counts_with_gene_freq, '../../Data/Result_data/Overlap_counts_with_gene_freq.csv')

