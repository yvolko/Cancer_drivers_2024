library(dplyr)
library(tidyr)
library(rvest)
library(tibble)
library(data.table)
library(stringr)
library(Hmisc)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)
library(enrichplot)
library(rbioapi)
library(GO.db)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# LOAD DATASETS

TCGA <- read.csv('TCGA_PPI.csv')
HARTWIG <- read.csv('HARTWIG_PPI.csv')
METASTATIC <- read.csv('../Metastatic_cancers_TF_kinase_CRISPR_oncoKB_COSMIC_pubmed_PPI.csv')

# GET LIST OF KNOWN ONCOGENES

COSMIC <- read.csv("../../../../Data/Raw_data/Cancer_Gene_census_COSMIC.tsv", sep = "\t")
oncoKB <- read.csv("../../../../Data/Raw_data/cancerGeneList_OncoKB.tsv", sep = "\t")

COSMIC <- COSMIC[grepl("oncogene", COSMIC$Role.in.Cancer), ]
oncoKB <- oncoKB[oncoKB$Is.Tumor.Suppressor.Gene!="Yes" & oncoKB$Is.Oncogene!="No", ]
known_oncogenes <- unique(c(COSMIC$Gene.Symbol, oncoKB$Hugo.Symbol))

# ADD PATHWAY ENRICHMENT
# Get all the genes in the analysed DFs
genes_from_4_df <- unique(c(TCGA$gene_name, HARTWIG$gene_name, METASTATIC$gene_name))
idx_all <- grch37$symbol %in% genes_from_4_df
ids_all <- grch37[idx_all, ]
non_duplicates <- which(duplicated(ids_all$symbol) == FALSE)
ids_all <- ids_all[non_duplicates, ]
all_genes <- as.character(ids_all$ensgene)

# Get all the known oncogenes
idx_onco <- grch37$symbol %in% known_oncogenes 
ids_onco <- grch37[idx_onco, ]
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


new_feature <- cluster_summary_onco[,c(3,9)]
for (i in cluster_summary_onco$Description){
  go_id = GOID( GOTERM[ Term(GOTERM) == i])
  get(go_id, org.Hs.egGO2ALLEGS)
  allegs = get(go_id, org.Hs.egGO2ALLEGS)
  genes = unlist(mget(allegs,org.Hs.egSYMBOL), use.names = FALSE)
  new_feature[new_feature$Description == i, 'geneID'] <- paste(unique(genes), sep = ',', collapse = '/')
}

new_feature <- new_feature %>%
  mutate(gene_list = strsplit(geneID, "/"))

GO_df <- data.frame()
for (i in genes_from_4_df){
  gene = i
  GO_terms = sum(sapply(new_feature$gene_list, function(x) i %in% x))
  GO_df <- rbind(GO_df, c(gene, GO_terms))
}

colnames(GO_df) <- c('gene', 'GO_terms')
GO_df$GO_terms <- as.numeric(GO_df$GO_terms)

TCGA_GO <- left_join(TCGA, GO_df, join_by(gene_name == gene))
HARTWIG_GO <- left_join(HARTWIG, GO_df, join_by(gene_name == gene))
METASTATIC_GO <- left_join(METASTATIC, GO_df, join_by(gene_name == gene))


write.csv(TCGA_GO, 'TCGA_with_GO.csv', row.names = FALSE)
write.csv(HARTWIG_GO, 'HARTWIG_with_GO.csv', row.names = FALSE)
write.csv(METASTATIC_GO, 'METASTATIC_with_GO.csv', row.names = FALSE)
