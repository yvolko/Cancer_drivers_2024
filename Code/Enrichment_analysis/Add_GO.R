library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)
library(enrichplot)
library(dplyr)


df <- read.csv('../../Data/Processed_data/Table_expression_COSMIC_oncoKB_merged_intervals.csv')

idx_all <- grch38$symbol %in% df$gene_name
ids_all <- grch38[idx_all, ]
non_duplicates <- which(duplicated(ids_all$symbol) == FALSE)
ids_all <- ids_all[non_duplicates, ]
all_genes <- as.character(ids_all$ensgene)

unique_oncogene <- unique(dplyr::filter(df, Role.in.Cancer.COSMIC == 'oncogene' | Role.in.Cancer.COSMIC == 'oncogene, fusion' |
                                          Is.Oncogene.oncoKB == 'Yes')$gene_name)
idx_onco <- grch38$symbol %in% unique_oncogene 
ids_onco <- grch38[idx_onco, ]
non_duplicates_onco <- which(duplicated(ids_onco$symbol) == FALSE)
ids_onco <- ids_onco[non_duplicates_onco, ] 
sig_genes_onco <- as.character(ids_onco$ensgene)


ego_ALL_onco <- enrichGO(gene = sig_genes_onco,
                         universe = all_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         ont = "ALL", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)

cluster_summary_onco <- data.frame(ego_ALL_onco)
write.table(cluster_summary_onco, '../../Data/Processed_data/cluster_summary_oncogene.txt', sep = ' ', row.names = FALSE)
