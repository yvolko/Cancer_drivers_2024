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

data <- read.csv("../3_interval_merging/Metastatic_cancers_with_gene_names_merged.csv")

# ADD INFORMATION ABOUT TRANSCRIPTIONAL FACTORS

TF <- read.csv("../../../github_structure/Data/Raw_data/Human_TFs.txt", sep="\t")
TF <- TF[, c("HGNC.symbol", "Is.TF.", "DBD")]
TF <- TF[TF$Is.TF.=="Yes", ]

data_TF <- left_join(data, TF, join_by("gene_name" == "HGNC.symbol"))
data_TF <- data_TF %>%
  mutate_at(c("Is.TF.", "DBD"), ~replace_na(.x, "No_data"))

# ADD INFORMATION ABOUT KINASES

webpage <- read_html("http://www.kinhub.org/kinases.html")
kinase <- html_table(webpage)[[1]]
kinase <- kinase[, c("xName", "HGNC Name", "Group", "Family")]
colnames(kinase)[2] <- "HGNC_Name"

# Some kinases have no name in "HGNC Name" column, we take it than from "xName" column
# Impute missing values

kinase <- kinase %>%
  mutate(HGNC_Name = ifelse(HGNC_Name == "", xName, HGNC_Name))

# Some of values were not unique due to different columns (with their name) being different in the full data set
# Let make them unique by using other column for naming of these

kinase$Name <- ifelse(duplicated(kinase$HGNC_Name) & !duplicated(kinase$xName), kinase$xName, kinase$HGNC_Name)
kinase <- kinase[, -c(1,2)]
colnames(kinase)[1] <- "Group_kinase"
colnames(kinase)[2] <- "Family_kinase"

# Add kinase to our data

data_TF_kinase <- left_join(data_TF, kinase, join_by("gene_name" == "Name"))
data_TF_kinase <- data_TF_kinase %>%
  mutate_at(c("Group_kinase", "Family_kinase"), ~replace_na(.x, "No_data"))

# ADD INFORMATION ABOUT CRISPR SCREENING (Columns: Gene - Rows: ScreenID)

CRISPR <- read.csv("../../../Task_3_TF_kinase_CRISPR_pubmed/CRISPRGeneEffect.csv")
CRISPR <- CRISPR[,-1]

# Calculate means/median/sd per column
median <- apply(CRISPR, 2, median, na.rm = TRUE)
mean <- colMeans(CRISPR, na.rm = TRUE)
sd <- apply(CRISPR, 2, sd, na.rm = TRUE)
min <- apply(CRISPR, 2, min, na.rm = TRUE)

# Assembl new CRISPR table
CRISPR_new <- data.frame(median, mean, sd, min)

CRISPR_new$median <- round(CRISPR_new$median, 3)
CRISPR_new$mean <- round(CRISPR_new$mean, 3)
CRISPR_new$sd <- round(CRISPR_new$sd, 3)
CRISPR_new$min <- round(CRISPR_new$min, 3)

CRISPR_new <- rownames_to_column(CRISPR_new, var = "Name")
CRISPR_new$Name <- sub("\\..*", "", CRISPR_new$Name)

colnames(CRISPR_new) <- paste(colnames(CRISPR_new), ".CRISPR", sep = "")

# Merge CRISPR screening data into our data

data_TF_kinase_CRISPR <- left_join(data_TF_kinase, CRISPR_new, join_by("gene_name" == "Name.CRISPR"))

# ADD COSMIC and oncoKB
COSMIC <- read.csv("../../../Task_2_statistics/Cancer_Gene_census_COSMIC.tsv", sep = "\t")
oncoKB <- read.csv("../../../Task_2_statistics/cancerGeneList_OncoKB.tsv", sep = "\t")

# KEEP ONLY RELATED COLUMNS

COSMIC <- COSMIC[, c("Gene.Symbol", "Tier", "Hallmark", "Role.in.Cancer", "Mutation.Types", 
                     "Translocation.Partner")]

colnames(COSMIC) <- paste(colnames(COSMIC), ".COSMIC", sep = "")
COSMIC <- COSMIC[grepl("TSG|oncogene", COSMIC$Role.in.Cancer), ]


oncoKB <- oncoKB[, c("Hugo.Symbol", "Is.Oncogene", "Is.Tumor.Suppressor.Gene", "OncoKB.Annotated")]

colnames(oncoKB) <- paste(colnames(oncoKB), ".oncoKB", sep = "")
oncoKB <- oncoKB[(oncoKB$Is.Tumor.Suppressor.Gene.oncoKB=="Yes"|oncoKB$Is.Oncogene=="Yes"), ]

# MERGE TABLES INTO OUR DATA

data_with_COSMIC <- left_join(data_TF_kinase_CRISPR, COSMIC, join_by("gene_name" == "Gene.Symbol.COSMIC"))
data_with_COSMIC_oncoKB <- left_join(data_with_COSMIC, oncoKB, join_by("gene_name" == "Hugo.Symbol.oncoKB"))

data_with_COSMIC_oncoKB$Tier.COSMIC <- as.character(data_with_COSMIC_oncoKB$Tier.COSMIC)
colnames_to_replace_na <- c("Tier.COSMIC", "Hallmark.COSMIC", 
                            "Role.in.Cancer.COSMIC", "Mutation.Types.COSMIC",
                            "Translocation.Partner.COSMIC", "Is.Oncogene.oncoKB",
                            "Is.Tumor.Suppressor.Gene.oncoKB", "OncoKB.Annotated.oncoKB")   
data_with_COSMIC_oncoKB <- data_with_COSMIC_oncoKB %>%
  mutate_at(colnames_to_replace_na, ~replace_na(.x, "No_data"))


# ADD PUBMED MINING DATA

pubmed <- read.csv("../../../Task_3_TF_kinase_CRISPR_pubmed/pubmed_data_260224.txt", sep="\t", header = TRUE)
data_pubmed <- left_join(data_with_COSMIC_oncoKB, pubmed, join_by("gene_name" == "gene"))

# ADD PPI (number means - how many interactions this gene has with the known oncogenes excluding TSG)
# If we would like to get all the interaction of the protein with every other protein in STRING db - use rba_string_interaction_partners
COSMIC <- read.csv("../../../Task_2_statistics/Cancer_Gene_census_COSMIC.tsv", sep = "\t")
oncoKB <- read.csv("../../../Task_2_statistics/cancerGeneList_OncoKB.tsv", sep = "\t")

COSMIC <- COSMIC[grepl("oncogene", COSMIC$Role.in.Cancer), ]
oncoKB <- oncoKB[oncoKB$Is.Tumor.Suppressor.Gene!="Yes" & oncoKB$Is.Oncogene!="No", ]
known_oncogenes <- unique(c(COSMIC$Gene.Symbol, oncoKB$Hugo.Symbol))

df <- data_pubmed
known_oncogenes_string_ids <- rba_string_map_ids(known_oncogenes, species = 9606)

query_string_ids <- rba_string_map_ids(unique(df$gene_name), species = 9606)


# NOTE: protein "METTL21B/EEF1AKMT3" with ID "9606.ENSP00000300209" was causing rba_string_interactions_network function to fail for unknown reason
# therefore it is removed from the query and substitute its number of interaction by 0
query_string_ids <- query_string_ids[-1572,]

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

temp_df <- data.frame(gene_id = "9606.ENSP00000300209", n_interactions = 0)
interactions_df <- rbind(interactions_df, temp_df)


# Add information about PPI count to our df
query_string_ids <- query_string_ids[,c(2,3)]
interactions_df_final <- left_join(interactions_df, query_string_ids, by=c("gene_id" = "stringId"))
interactions_df_final$gene_id <- NULL
df <- left_join(df, interactions_df_final, by=c("gene_name" = "queryItem"))
colnames(df)[41] <- 'PPI_count'

# File TCGA_Hartwig_PPI.R adds PPI information to TCGA and Hartwig. Use csv HARTWIG_PPI.csv and TCGA_PPI.csv for further GO results


# Additional test about PPI

# Test if known oncogenes have a lot of interactions with other known oncogenes
oncogenes_interactions_df <- data.frame()
for (oncogene in known_oncogenes_string_ids$stringId){
  ids = c(known_oncogenes_string_ids$stringId, oncogene)
  int_net <- rba_string_interactions_network(ids = ids,
                                             species = 9606,
                                             required_score = 700)
  
  int_net_final <- subset(int_net, stringId_A == oncogene | stringId_B == oncogene)
  
  print(c(oncogene, nrow(int_net_final)))
  temp_df <- data.frame(gene_id = oncogene, n_interations = nrow(int_net_final))
  oncogenes_interactions_df <- rbind(oncogenes_interactions_df, temp_df)
}

print(c(mean(oncogenes_interactions_df$n_interations), sd(oncogenes_interactions_df$n_interations)))
# -> 16.44489 20.28577 Mean ans SD of PPI with other known oncogenes for all known oncogenes

interactions_df_non_oncogenes <- anti_join(interactions_df, oncogenes_interactions_df)
print(c(mean(interactions_df_non_oncogenes$n_interactions), sd(interactions_df_non_oncogenes$n_interactions)))
# -> 1.473620 4.552767 Mean ans SD of PPI with other known oncogenes for non oncogenes from our metastatic cancers df

boxplot(oncogenes_interactions_df$n_interations, interactions_df_non_oncogenes$n_interactions,
        names = c("Oncogenes", "Non oncogenes"),
        main = "Number of PPI with known oncogenes",
        ylab = "Count of PPI",
        col = c("lightblue", "lightgreen"))


combined_df <- data.frame(
  group = c(rep("oncogenes", length(oncogenes_interactions_df$n_interations)),  
            rep("non-oncogenes", length(interactions_df_non_oncogenes$n_interactions))), 
  value = c(oncogenes_interactions_df$n_interations, interactions_df_non_oncogenes$n_interactions)  
)

ggplot(combined_df, aes(x = log(value+1), fill = group)) +
  geom_density(alpha = 0.5, position = "identity", color = "grey") +  
  scale_fill_manual(values = c("lightblue", "lightgreen")) + 
  labs(title = "Density Plot PPI with known oncognes for oncogenes and non oncogenes",
       x = "log(count of PPI)",
       y = "Density") +
  theme_minimal()



# ADD PATHWAY ENRICHMENT (!!!important to run it for all the df together as all_genes should be taken from all four datasets)

TCGA <- read.csv('')
HARTWIG <- read.csv('')

all_genes_from_4_df <- unique(c(df$gene_name, TCGA$gene_name, HARTWIG$gene_name))

idx_all <- grch37$symbol %in% all_genes_from_4_df
ids_all <- grch37[idx_all, ]
non_duplicates <- which(duplicated(ids_all$symbol) == FALSE)
ids_all <- ids_all[non_duplicates, ]
all_genes <- as.character(ids_all$ensgene)

unique_oncogene <- known_oncogenes

idx_onco <- grch37$symbol %in% unique_oncogene 
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
for (i in unique(df$gene_name)){
  gene = i
  GO_terms = sum(sapply(new_feature$gene_list, function(x) i %in% x))
  GO_df <- rbind(GO_df, c(gene, GO_terms))
}

colnames(GO_df) <- c('gene', 'GO_terms')
df_with_GO <- left_join(df, GO_df, join_by(gene_name == gene))

# ADD PPI COUNT AND GO TERMS TO OLD TABLE
# GO TERMS
old_df <- read.csv("../result_files/Table_expr_merged_TF_kinase_CRISPR_min.csv")


GO_df <- data.frame()
for (i in unique(old_df$gene_name)){
  gene = i
  GO_terms = sum(sapply(new_feature$gene_list, function(x) i %in% x))
  GO_df <- rbind(GO_df, c(gene, GO_terms))
}

colnames(GO_df) <- c('gene', 'GO_terms')
old_df_with_GO <- left_join(old_df, GO_df, join_by(gene_name == gene))
  
# PPI  
  
query_string_ids <- rba_string_map_ids(unique(old_df_with_GO$gene_name), species = 9606)

# Remove "9606.ENSP00000300209" from query ids

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

temp_df <- data.frame(gene_id = "9606.ENSP00000300209", n_interactions = 0)
interactions_df <- rbind(interactions_df, temp_df)

# Add information about PPI count to our df
query_string_ids <- query_string_ids[,c(2,3)]
interactions_df_final <- left_join(interactions_df, query_string_ids, by=c("gene_id" = "stringId"))
interactions_df_final$gene_id <- NULL
old_df_with_GO_PPI <- left_join(old_df_with_GO, interactions_df_final, by=c("gene_name" = "queryItem"))
colnames(old_df_with_GO_PPI)[37] <- 'PPI_count'


write.csv(old_df_with_GO_PPI, "Non_matastatic_annotated.csv", row.names = FALSE)  








# ADD EXPRESSION DATA TO METASTATIC DF

# MET500

# DATA LOADING
df <- read.csv('Metastatic_annotated.csv')

tpm_MET500 <- read.csv("MET-PRISM_MET500_expression/MET500/MET500_TPM.txt", sep = "\t")
additional_metadata_MET500 <- read.csv('MET-PRISM_MET500_expression/MET500/SraRunTable_MET500_additional_metadata.csv')

# DATA PREPARATION
# Choose from tpm table only genes we have in our data
gene_list_MET500 <- which(tpm_MET500$Gene.name %in% table_with_gene_name$gene_name)
tpm_MET500_selection <- tpm_MET500[gene_list_MET500,]

# Make table long and combine with corresponding annotation data
tpm_MET500_long <- tpm_MET500_selection %>%
  pivot_longer(cols = -c(Gene.name), names_to = 'sample', values_to = 'expr')


# ADD METADATA (seems like we don't need annotation_RNA_MET500 and annotation_WES_MET500 anymore)
tpm_MET500_long_metadata <- left_join(tpm_MET500_long, additional_metadata_MET500, join_by('sample' == 'Run')) 
tpm_MET500_long_metadata <- tpm_MET500_long_metadata[, c(1, 2, 3, 11, 22, 25, 32, 33, 38, 39, 40, 41)]

# FROM METADATA COMBINE SRR ids for DNA and RNA via using submitted_subject_id
selection <- additional_metadata_MET500[, c(1, 2, 36, 37, 38)]

summary <- selection %>%
  group_by(submitted_subject_id) %>%
  summarise(
    n_RNA = sum(analyte_type == 'RNA'),
    n_DNA = sum(analyte_type == 'DNA')
  ) %>%
  ungroup()

# Some of the submitted subject ids don't have DNA sample, some submitted subject ids have 2 DNA, etc.

result <- selection %>%
  group_by(submitted_subject_id) %>%
  summarise(
    DNA_Runs = list(Run[analyte_type == "DNA"]),
    RNA_Runs = list(Run[analyte_type == "RNA"])
  ) %>%
  ungroup()


# Modify results so each DNA id gets corresponding RNA ids

expanded_result <- result %>%
  unnest(DNA_Runs) %>% # Expand DNA_Run_ids so each gets its own row
  unnest_longer(RNA_Runs) # Repeat each DNA_Run_id for all RNA_Run_ids

# Check if multiple RNA seq ids available for the same DNA, select only ones with PolyA LibrarySelection method
LibrarySelection <- additional_metadata_MET500[, c(1, 23)]

expanded_result <- left_join(expanded_result, LibrarySelection, join_by('RNA_Runs' == 'Run'))

# Make summary on which RNA samples belong to same DNA sample, also select the ones with PolyA if possible

expanded_result <- expanded_result %>%
  group_by(DNA_Runs) %>%
  summarise(
    RNA_Runs = if (any(LibrarySelection == "PolyA")) {
      paste(RNA_Runs[LibrarySelection == "PolyA"], collapse = ",")
    } else {
      paste(RNA_Runs, collapse = ",")
    }
  ) %>%
  ungroup()

# PREPARE OUR DATA
# Select only data for MET500 from our df and split sample names into cancer and normal tissue sample names
MET500_df <- df[df$dataset=='MET500', ]
MET500_df <- MET500_df %>%
  separate(sample, into = c("sample_cancer", "sample_normal"), sep = "_Vs_")

# For each cancer WES SRR sample add corresponding RNA seq samples SRR
MET500_df_with_RNA_SRRs <- left_join(MET500_df, expanded_result, join_by('sample_cancer' == 'DNA_Runs'))

# Add tpm values to the MET500 df
MET500_df_with_RNA_SRRs <- left_join(MET500_df_with_RNA_SRRs, tpm_MET500_long_metadata, join_by('RNA_Runs' == 'sample', 'gene_name' == 'Gene.name'))

MET500_df_with_RNA_SRRs$RNA_Runs <- NULL
MET500_df_with_RNA_SRRs$biospecimen_repository_sample_id <- NULL
MET500_df_with_RNA_SRRs$Is_Tumor <- NULL
MET500_df_with_RNA_SRRs$LibrarySelection <- NULL
MET500_df_with_RNA_SRRs$Sample.Name <- NULL
MET500_df_with_RNA_SRRs$SRA.Study <- NULL
MET500_df_with_RNA_SRRs$submitted_subject_id <- NULL
MET500_df_with_RNA_SRRs$sex <- NULL

write.csv(MET500_df_with_RNA_SRRs, 'MET-PRISM_MET500_expression/MET500/MET500_with_expr.csv')



# METAPRISM
# Select only data for METAPRISM from our df and split sample names into cancer and normal tissue sample names
df_METAPRISM <- df[df$dataset=='METAPRISM', ]
df_METAPRISM <- df_METAPRISM %>%
  separate(sample, into = c("sample_cancer", "sample_normal"), sep = "_Vs_")
df_METAPRISM$sample_cancer <- gsub("^Sample_", "", df_METAPRISM$sample_cancer)
df_METAPRISM$sample_normal <- gsub("^Sample_", "", df_METAPRISM$sample_normal)

# Load expression data

tpm_METAPRISM <- read.csv("MET-PRISM_MET500_expression/METAPRISM/METAPRISM_WHOLE_TPM.txt", sep = '\t')
annotation_METAPRISM <- read.csv("MET-PRISM_MET500_expression/METAPRISM/METAPRISM_annotation_281124.tsv", sep = '\t')
annotation_METAPRISM_short <- annotation_METAPRISM[, c(6, 7, 12, 13, 14, 17, 18, 19, 20, 21, 30, 32, 33, 34, 35)]
clinical_data_METAPRISM <- read.csv('MET-PRISM_MET500_expression/METAPRISM/metaprism_2023_clinical_data.tsv', sep = '\t')

# DATA PREPARATION
# Choose only genes we have in our data
gene_list_METAPRISM <- which(tpm_METAPRISM$Gene.name %in% unique(df_METAPRISM$gene_name))
tpm_METAPRISM <- tpm_METAPRISM[gene_list_METAPRISM,]

# Make table long
tpm_METAPRISM_long <- tpm_METAPRISM %>%
  pivot_longer(cols = -c(Gene.name), names_to = 'sample', values_to = 'expr')

# Merge tpm and annotation
tpm_METAPRISM_long$sample <- gsub("\\.", "-", tpm_METAPRISM_long$sample)
annotation_METAPRISM <- annotation_METAPRISM[, c(2, 3, 7, 8, 9, 12, 15, 16, 28, 32, 33, 37, 49)]


# Figure out how samples in our df and annotation and tpm are named!!! Next step

tpm_METAPRISM_long_annotated <- left_join(tpm_METAPRISM_long, annotation_METAPRISM, join_by('sample' == 'RNAseq.Sample.ID'))
tpm_METAPRISM_long_annotated_selected <- tpm_METAPRISM_long_annotated[, c(1, 2, 3, 11, 12, 22)]


# Select only data for METAPRISM from our df
METAPRISM_df <- df[df$dataset=='METAPRISM', ]







# Load GTEX prepared table to get normal tissue expression values
GTEX <- read.csv('MET-PRISM_MET500_expression/GTEX/GTEX_normal_tissues_expr_summary.csv')

# Select genes we need
genes_MET500 <- unique(MET500_df_with_RNA_SRRs$gene_name)
genes_METAPRISM <- unique(tpm_METAPRISM$Gene.name)
genes <- unique(c(genes_MET500, genes_METAPRISM))

genes_indexes <- which(GTEX_tpm_long_gene_names$gene %in% genes)
GTEX_tpm_long_gene_names <- GTEX_tpm_long_gene_names[genes_indexes, ]

GTEX_tpm_long_gene_names <- left_join(GTEX_tpm_long_gene_names, GTEX_pheno, join_by('sample_Id' == 'Sample'))

# Calculate mean and std for each site and each gene in GTEX
GTEX_summary_primary_site <- GTEX_tpm_long_gene_names %>%
  group_by(`_primary_site`, gene) %>%
  summarize(
    mean = mean(expr),
    std = sd(expr),
    n = n()
    ) %>%
  ungroup()

GTEX_summary_primary_site <- GTEX_summary_primary_site[GTEX_summary_primary_site$`_primary_site` != '<not provided>', ]


# Add information from GTEX (norm expr) to MET500 combining by the corresponding tissues
merging_table <- read.csv('correspondance_of_tissues.txt', sep='\t')

df_MET500_with_expr_norm <- left_join(df_MET500_with_expr_norm, )
























# DATA MERGING
table_gene_tpm_data <- table_gene_tpm_data %>%
  left_join(annotation_gene_list, by = 'sample') 
normal_table_gene_tpm_data <- colnames(table_gene_tpm_data)[grep("-11$", colnames(table_gene_tpm_data))]
table_gene_tpm_data_norm <- table_gene_tpm_data[,..normal_table_gene_tpm_data]
table_gene_tpm_data_norm$gene_name <- table_gene_tpm_data$gene
table_gene_tpm_data_norm$sample <- table_gene_tpm_data$sample
table_gene_tpm_data_norm_long <- table_gene_tpm_data_norm %>%
  pivot_longer(cols = -c(gene_name, sample), names_to = 'sample_cut', values_to = 'norm_expr')
table_gene_tpm_data_norm_long <- left_join(table_gene_tpm_data_norm_long, sample_anot, join_by('sample_cut' == 'sample'))


colnames(table_gene_tpm_data_norm_long)[colnames(table_gene_tpm_data_norm_long) == "sample"] <- "ENS_id"
table_gene_tpm_data_norm_long$sample <- table_gene_tpm_data_norm_long$sample_cut
table_gene_tpm_data_norm_long$sample_cut <- str_sub(table_gene_tpm_data_norm_long$sample_cut, 1, 12)





result_df_with_norm <- left_join(result_df, table_gene_tpm_data_norm_long, by = c('sample_cut', 'gene_name', 'ENS_id'))
colnames(result_df_with_norm)[colnames(result_df_with_norm) == "sample_type.x"] <- "sample_type"
colnames(result_df_with_norm)[colnames(result_df_with_norm) == "X_primary_disease.x"] <- "Primary_disease"
#result_df_with_norm <- result_df_with_norm[,-c(17,18,19)]
result_df_with_norm <- result_df_with_norm %>%
  group_by(sample_cut, gene_name) %>%
  filter(expr == max(expr) | is.na(expr))

for (i in 1:4423){
  if (is.na(result_df_with_norm[i,]$norm_expr) & !is.na(result_df_with_norm[i,]$expr)){
    if (result_df_with_norm[i,]$Primary_disease == "brain lower grade glioma"){
      result_df_with_norm[i,]$norm_expr <- mean(filter(table_gene_tpm_data_norm_long, 
                                                       table_gene_tpm_data_norm_long$X_primary_disease ==  'glioblastoma multiforme'& 
                                                         table_gene_tpm_data_norm_long$gene_name == result_df_with_norm[i,]$gene_name & 
                                                         table_gene_tpm_data_norm_long$ENS_id == result_df_with_norm[i,]$ENS_id)$norm_expr)
    }
    if (result_df_with_norm[i,]$Primary_disease == "uterine carcinosarcoma"){
      result_df_with_norm[i,]$norm_expr <- mean(filter(table_gene_tpm_data_norm_long, 
                                                       table_gene_tpm_data_norm_long$X_primary_disease ==  'uterine corpus endometrioid carcinoma'& 
                                                         table_gene_tpm_data_norm_long$gene_name == result_df_with_norm[i,]$gene_name & 
                                                         table_gene_tpm_data_norm_long$ENS_id == result_df_with_norm[i,]$ENS_id)$norm_expr)
    }
    else{
      result_df_with_norm[i,]$norm_expr <- mean(filter(table_gene_tpm_data_norm_long, 
                                                       table_gene_tpm_data_norm_long$X_primary_disease == result_df_with_norm[i,]$Primary_disease & 
                                                         table_gene_tpm_data_norm_long$gene_name == result_df_with_norm[i,]$gene_name & 
                                                         table_gene_tpm_data_norm_long$ENS_id == result_df_with_norm[i,]$ENS_id)$norm_expr)
    }
  }
}

# result_df_with_norm <- result_df_with_norm[,-c(17,18,19)]
result_df_with_norm$norm_expr <- round(result_df_with_norm$norm_expr, digits = 2)
result_df_with_norm_na <- result_df_with_norm %>%
  mutate(norm_expr = ifelse(is.nan(norm_expr), NA, norm_expr))

for (i in 1:48860){
  if (is.na(result_df_with_norm_na[i,]$norm_expr) & !is.na(result_df_with_norm_na[i,]$expr)){
    if (result_df_with_norm_na[i,]$Primary_disease == "brain lower grade glioma"){
      result_df_with_norm_na[i,]$norm_expr <- mean(filter(table_gene_tpm_data_norm_long, 
                                                          table_gene_tpm_data_norm_long$X_primary_disease ==  'glioblastoma multiforme'& 
                                                            table_gene_tpm_data_norm_long$gene_name == result_df_with_norm_na[i,]$gene_name & 
                                                            table_gene_tpm_data_norm_long$ENS_id == result_df_with_norm_na[i,]$ENS_id)$norm_expr)
    }
    if (result_df_with_norm_na[i,]$Primary_disease == "uterine carcinosarcoma"){
      result_df_with_norm_na[i,]$norm_expr <- mean(filter(table_gene_tpm_data_norm_long, 
                                                          table_gene_tpm_data_norm_long$X_primary_disease ==  'uterine corpus endometrioid carcinoma'& 
                                                            table_gene_tpm_data_norm_long$gene_name == result_df_with_norm_na[i,]$gene_name & 
                                                            table_gene_tpm_data_norm_long$ENS_id == result_df_with_norm_na[i,]$ENS_id)$norm_expr)
    }
  }
}

result_df_with_norm_na$norm_expr <- round(result_df_with_norm_na$norm_expr, digits = 2)
result_df_with_norm_na <- result_df_with_norm_na[,-11]
write.csv(result_df_with_norm_na, '../../Data/Processed_data/table_with_gene_expr.csv', row.names = FALSE)


# EXPORT FINAL DATAFRAME

data_TF_kinase_CRISPR$dif_interval <- NULL
write.csv(data_with_COSMIC_oncoKB, file = "Metastatic_cancers_with_gene_names_from_gencode_merged_annotated.csv", row.names = FALSE)
