library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)
library(enrichplot)


top_not_oncogene <- read.csv('../../Data/Result_data/Overlap_counts_with_gene_freq.csv')
df_with_annotation <- read.csv('../../Data/Processed_data/Table_expression_COSMIC_oncoKB_merged_intervals.csv')

idx_not <- grch38$symbol %in% unique(top_not_oncogene$gene_name)
ids_not <- grch38[idx_not, ]
non_duplicates_not <- which(duplicated(ids_not$symbol) == FALSE)
ids_not <- ids_not[non_duplicates_not, ] 
sig_genes <- as.character(ids_not$ensgene)

idx_all <- grch38$symbol %in% unique(df_with_annotation$gene_name)
ids_all <- grch38[idx_all, ]
non_duplicates <- which(duplicated(ids_all$symbol) == FALSE)
ids_all <- ids_all[non_duplicates, ]
all_genes <- as.character(ids_all$ensgene)

ego_ALL_not_onco <- enrichGO(gene = sig_genes,
                             universe = all_genes,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             ont = "ALL", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)

cluster_summary_not_onco <- data.frame(ego_ALL_not_onco)
write.csv(cluster_summary_not_onco, file = "../../Data/Result_data/GO_terms_for_potential_oncogenes.csv", row.names = FALSE)
