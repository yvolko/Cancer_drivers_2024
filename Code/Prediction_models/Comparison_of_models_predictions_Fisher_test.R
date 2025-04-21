library(tidyr)
library(dplyr)
library(VennDiagram)

# COMBINE PREDICTIONS FROM ALL MODEL AND TRY TO FIND A CRITERIA TO REDUCE FALSE CALLS
# LOAD ALL THE DATA
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

predictions_RF <- read.csv('../RF/Predictions_by_RandomForestClassifier_TCGA_METASTATIC_HARTWIG.csv')
predictions_NN_counts <- read.csv('../NN/Oncogenes_predicted_by_NN_with_loc_and_count.tsv', sep = '\t')

# Keep only genes that were predicted as oncognes with probability >= 0.95
# Add localizations where oncognes were predicted
predictions_RF <- predictions_RF[predictions_RF$prob_of_onco >= 0.95, ]

localizations_TCGA <- read.csv('../Paper_writing/TCGA_with_cancer_types.csv')
localizations_MET500 <- read.csv('../Paper_writing/MET500_with_cancer_types.csv')
localizations_METAPRISM <- read.csv('../Paper_writing/METAPRISM_with_cancer_types.csv')
localizations_Hartwig <- read.csv('../Paper_writing/Hartwig_with_cancer_types.csv')

localizations_MET500$sample <- paste(localizations_MET500$tumor_sample, localizations_MET500$blood_sample, sep = '_Vs_')
colnames(localizations_MET500)[colnames(localizations_MET500) == "Project_TCGA_More"] <- "cancer_type"
  
localizations_METAPRISM$sample <- paste(localizations_METAPRISM$sample.id.tumor, localizations_METAPRISM$sample.id.normal, sep = '_Vs_')
colnames(localizations_METAPRISM)[colnames(localizations_METAPRISM) == "Cancer.Cohort"] <- "cancer_type"

colnames(localizations_Hartwig)[colnames(localizations_Hartwig) == "cancer_type_code"] <- "cancer_type"

localizations_TCGA <- unique(localizations_TCGA[, c('sample', 'cancer_type')])
localizations_MET500 <- unique(localizations_MET500[, c('sample', 'cancer_type')])
localizations_METAPRISM <- unique(localizations_METAPRISM[, c('sample', 'cancer_type')])
localizations_Hartwig <- unique(localizations_Hartwig[, c('sample', 'cancer_type')])

all_localizations <- rbind(localizations_TCGA, localizations_MET500, localizations_METAPRISM, localizations_Hartwig)
  
predictions_RF <- left_join(predictions_RF, all_localizations, join_by('sample'))

predictions_RF_counts <- predictions_RF %>%
  group_by(gene_name) %>%
  summarise(
    count = n(),
    localizations = paste(unique(cancer_type), collapse = ', ')
  ) %>%
  ungroup()
  
write.csv(predictions_RF_counts, 'Oncogene_predicted_by_RF_0.95_with_loc_and_count.csv', row.names = F)

# SEE OVERLAP BETWEEN PREDICTIONS
VennDiagramm <- venn.diagram(
  x = list(unique(NN_predictions_95$X0), unique(RF_predictions_95$gene_name)),
  category.names = c("Neural\nNetwork", "Random\nForest\nClassifier"),
  filename = NULL,
  output = TRUE,
  fill = c('#3C5488FF', '#DC0000FF'),
  cat.cex = 2.5,              # Adjust category label size
  cat.fontsize = 20,          # Adjust category label font size
  cex = 2,                    # Adjust text size inside the circles
  scaled = TRUE,              # Automatic scaling of circles
  lwd = 0,
  col = NA
)


grid.newpage()
grid.draw(VennDiagramm)


# COMBINE ALL PREDICTIONS
colnames(predictions_RF_counts) <- c("gene_name", "count_RF", "localizations_RF")
colnames(predictions_NN_counts) <- c("gene_name", "count_NN", "localizations_NN")

all_predictions <- full_join(predictions_RF_counts, predictions_NN_counts, by = 'gene_name')

all_predictions$count_RF <- ifelse(is.na(all_predictions$count_RF), 0, all_predictions$count_RF)
all_predictions$count_NN <- ifelse(is.na(all_predictions$count_NN), 0, all_predictions$count_NN)

all_predictions$total_count <- all_predictions$count_RF + all_predictions$count_NN

write.csv(all_predictions, 'Combined_oncogene_predictions.csv', row.names = F)






# CALCULATE CO_OCURENCE WITH KNOWN ONCOGENES
full_df <- read.csv("../Add_ranks/from Yurii/Table_with_rank_without_dupl.csv")
full_df$is_oncogene <-  ifelse(full_df$Role.in.Cancer.COSMIC %in% c("oncogene", "oncogene, fusion") | full_df$Is.Oncogene.oncoKB == "Yes", 1, 0)

samples_with_known_oncognes <- full_df %>%
  group_by(sample) %>%
  summarise(
    n_known_oncogenes = sum(is_oncogene),
    known_oncogene_names = paste(gene_name[is_oncogene == 1], collapse = ", ")
  ) %>%
  ungroup()

summary_by_gene <- full_df %>%
  group_by(gene_name) %>%
  summarise(
    found_in_samples = paste(unique(sample), collapse = ",")
  ) %>%
  ungroup()

# Make table that helps to see if gene is often found in the samples with known oncogenes
new_df <- data.frame(matrix(ncol = 5, nrow = 0))

for (i in 1:nrow(summary_by_gene)) {
  # Iterate row by row
  cur_df <- as.data.frame(summary_by_gene[i,])
  
  # Take one gene and all the samples it is present in
  samples <- cur_df[,c(2)]
  samples <- unlist(strsplit(samples, ","))
  
  cur_gene <- cur_df[,c(1)]
  
  # Calculate how many samples in total, how many have and don't have oncogene(s)
  cur_sample_table <- samples_with_known_oncognes[samples_with_known_oncognes$sample %in% samples,]
  n_samples_with_known_oncognes <- length(cur_sample_table$sample[cur_sample_table$n_known_oncogenes>0])
  n_samples_without_known_oncognes <- length(cur_sample_table$sample[cur_sample_table$n_known_oncogenes==0])
  oncogene_names <- paste(cur_sample_table$known_oncogene_names[cur_sample_table$known_oncogene_names != ''], collapse = '/')
  
  # Make a table that summarize information about each gene
  cur_new_df <- as.data.frame(t(c(cur_gene, 
                                  length(cur_sample_table$sample), 
                                  n_samples_with_known_oncognes, 
                                  n_samples_without_known_oncognes, 
                                  oncogene_names)))
  
  new_df <- rbind(new_df, cur_new_df)
}

colnames(new_df) <- c('gene', 'total_n_samples_in_4_dfs', 'n_samples_with_known_oncognes', 'n_samples_without_known_oncognes', 'oncogene_names_that_gene_co_found_with')


# RUN FISHER TEST
samples_without_any_oncognes <- length(samples_with_known_oncognes$sample[samples_with_known_oncognes$n_known_oncogenes == 0])
samples_with_oncogenes <- length(samples_with_known_oncognes$sample[samples_with_known_oncognes$n_known_oncogenes > 0])
predictions_all$fisher_result_p_values <- 0

for (i in 1:nrow(predictions_all)) {
  cur_df <- predictions_all[i,]
  gene = cur_df$gene_name
  cur_new_df <- new_df[new_df$gene == gene, ]
  
  N_pr_onco = as.numeric(cur_new_df$n_samples_with_known_oncognes)
  N_pr_only = as.numeric(cur_new_df$n_samples_without_known_oncognes)
  N_onco_only = samples_with_oncogenes - N_pr_onco
  N_none = samples_without_any_oncognes - N_pr_only
  contingency_table <- matrix(c(N_pr_onco, N_pr_only, N_onco_only, N_none), 
                              nrow = 2, 
                              byrow = TRUE,
                              dimnames = list(
                                A = c("Prdicted_yes", "Prdicted_no"),
                                B = c("Known_onco_yes", "Known_onco_no")
                              ))
  
  fisher_result <- fisher.test(contingency_table)
  predictions_all[i,3] <- fisher_result$p.value
}

# HOW MANY POTENTIAL ONCOGENES LEFT
length(predictions_all$gene_name[predictions_all$fisher_result_p_values < 0.01])

predictions_all$gene_name[predictions_all$fisher_result_p_values < 0.01]

# DO SAME FISHER TEST FOR KNOWN ONCOGENES
oncogenes <- unique(all_data$gene_name[all_data$is_oncogene == 1])
oncogenes_df <- as.data.frame(oncogenes)
oncogenes_df$fisher_result_p_values <- 0

for (i in 1:nrow(oncogenes_df)) {
  cur_onco_df <- oncogenes_df[i,]
  oncogene = cur_onco_df$oncogenes
  cur_new_onco_df <- new_df[new_df$gene == oncogene, ]
  
  N_pr_onco = as.numeric(cur_new_onco_df$n_samples_with_known_oncognes)
  N_pr_only = as.numeric(cur_new_onco_df$n_samples_without_known_oncognes)
  N_onco_only = samples_with_oncogenes - N_pr_onco
  N_none = samples_without_any_oncognes - N_pr_only
  contingency_table_onco <- matrix(c(N_pr_onco, N_pr_only, N_onco_only, N_none), 
                              nrow = 2, 
                              byrow = TRUE,
                              dimnames = list(
                                A = c("Prdicted_yes", "Prdicted_no"),
                                B = c("Known_onco_yes", "Known_onco_no")
                              ))
  
  fisher_result <- fisher.test(contingency_table)
  predictions_all[i,3] <- fisher_result$p.value
}

