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

# IMPORT DATA

HARTWIG <- read.csv("HARTWIG_TF_kinase_CRISPR_expr.csv")

TCGA <- read.csv('df_with_GO_and_PPI.csv', sep = '\t', row.names = NULL)
TCGA$X <- NULL
TCGA$GO_terms <- NULL
TCGA$PPI_count <- NULL

# ADD PPI (number means - how many interactions this gene has with the known oncogenes excluding TSG)
COSMIC <- read.csv("../../../../Data/Raw_data/Cancer_Gene_census_COSMIC.tsv", sep = "\t")
oncoKB <- read.csv("../../../../Data/Raw_data/cancerGeneList_OncoKB.tsv", sep = "\t")

COSMIC <- COSMIC[grepl("oncogene", COSMIC$Role.in.Cancer), ]
oncoKB <- oncoKB[oncoKB$Is.Tumor.Suppressor.Gene!="Yes" & oncoKB$Is.Oncogene!="No", ]
known_oncogenes <- unique(c(COSMIC$Gene.Symbol, oncoKB$Hugo.Symbol))

known_oncogenes_string_ids <- rba_string_map_ids(known_oncogenes, species = 9606)


# TAKE PPI FOR THE GENES THAT ARE ALREADY KNOWN FOR OTHER DATASETS
METASTATIC <- read.csv('../Metastatic_cancers_TF_kinase_CRISPR_oncoKB_COSMIC_pubmed_PPI.csv')
METASTATIC_PPI <- unique(METASTATIC[, c(24, 46)])


# MAKE QUERY LIST EXCLUDING ALREADY EXISTING DATA IN TCGA AND METASTATIC CANCERS
query <- unique(c(HARTWIG$gene_name, TCGA$gene_name))
query <- setdiff(query, METASTATIC_PPI$gene_name)
query_string_ids <- rba_string_map_ids(query, species = 9606)

# SERACH FOR INTERACTIONS
interactions_df <- data.frame()
for (id in query_string_ids$stringId){
  ids = c(known_oncogenes_string_ids$stringId, id)
  int_net <- rba_string_interactions_network(ids = ids,
                                             species = 9606,
                                             required_score = 700)

  int_net_final <- subset(int_net, stringId_A == id | stringId_B == id)
  
  print(c(id, nrow(int_net_final)))
  temp_df <- data.frame(gene_id = id, n_interactions = nrow(int_net_final))
  interactions_df <- rbind(interactions_df, temp_df)
}

# Add information about PPI count to TCGA, HARTWIG
query_string_ids <- query_string_ids[,c(2,3)]
interactions_df_final <- left_join(interactions_df, query_string_ids, by=c("gene_id" = "stringId"))
interactions_df_final$gene_id <- NULL
colnames(METASTATIC_TCGA_PPI) <- c("queryItem", "n_interactions")
interactions_df_final <- rbind(interactions_df_final, METASTATIC_TCGA_PPI)

HARTWIG <- left_join(HARTWIG, interactions_df_final, by=c("gene_name" = "queryItem"))
colnames(HARTWIG)[32] <- 'PPI_count'
HARTWIG$PPI_count <- ifelse(is.na(HARTWIG$PPI_count), 0, HARTWIG$PPI_count)

TCGA <- left_join(TCGA, interactions_df_final, by=c("gene_name" = "queryItem"))
colnames(TCGA)[38] <- 'PPI_count'
TCGA$PPI_count <- ifelse(is.na(TCGA$PPI_count), 0, TCGA$PPI_count)


write.csv(HARTWIG, 'HARTWIG_PPI.csv', row.names = F)
write.csv(TCGA, 'TCGA_PPI.csv', row.names = F)
