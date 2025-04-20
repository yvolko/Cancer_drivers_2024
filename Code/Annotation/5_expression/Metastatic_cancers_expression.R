library(dplyr)
library(tidyr)
library(rvest)
library(tibble)
library(data.table)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ADD EXPRESSION DATA TO METASTATIC DF

# IMPORT DATA

df <- read.csv("../4_annotation_with_DBs/TCGA_and_Hartwig/METASTATIC_with_GO.csv")

# MET500
# Select only data for MET500 from our df
MET500_df <- df[df$dataset=='MET500', ]

tpm_MET500 <- read.csv("../../raw_data/expression/MET500/MET500_TPM.txt", sep = "\t")
metadata_MET500 <- read.csv('../../raw_data/annotation/MET500/cln_met500_in_design_curated.tsv', sep = "\t")

# DATA PREPARATION
# Choose from tpm table only genes we have in our data
gene_list_MET500 <- which(tpm_MET500$Gene.name %in% MET500_df$gene_name)
tpm_MET500_selection <- tpm_MET500[gene_list_MET500,]

# Make table long and combine with corresponding annotation data
tpm_MET500_long <- tpm_MET500_selection %>%
  pivot_longer(cols = -c(Gene.name), names_to = 'sample', values_to = 'expr')


# ADD METADATA TO TPM TABLE AND THEN TPM TABLE TO OUR DF
tpm_MET500_long_metadata <- left_join(tpm_MET500_long, metadata_MET500, join_by('sample' == 'Sample_Id_RNA_T')) 

MET500_df <- left_join(MET500_df, tpm_MET500_long_metadata[, c(1, 3, 10, 11, 21, 22)], join_by('sample.id.normal' == 'Sample_Id_DNA_N', 'gene_name' == 'Gene.name'))



# METAPRISM
# Select only data for METAPRISM from our df
METAPRISM_df <- df[df$dataset=='METAPRISM', ]

# Load expression data

tpm_METAPRISM <- read.csv("../../raw_data/expression/METAPRISM/METAPRISM_WHOLE_TPM.txt", sep = '\t')
annotation_METAPRISM <- read.csv("../../raw_data/annotation/METAPRISM/metaprism_2023_clinical_data.tsv", sep = '\t')
annotation_METAPRISM_short <- annotation_METAPRISM[, c(3, 8, 9, 12, 15, 16, 19, 28, 32, 33, 37, 49)]

# DATA PREPARATION
# Choose only genes we have in our data
gene_list_METAPRISM <- which(tpm_METAPRISM$Gene.name %in% unique(METAPRISM_df$gene_name))
tpm_METAPRISM <- tpm_METAPRISM[gene_list_METAPRISM,]

# Make table long and rename samples to match format of the annotation
tpm_METAPRISM_long <- tpm_METAPRISM %>%
  pivot_longer(cols = -c(Gene.name), names_to = 'sample', values_to = 'expr')
tpm_METAPRISM_long$sample <- gsub("\\.", "-", tpm_METAPRISM_long$sample)

# Merge tpm and annotation
annot_samples <- unique(annotation_METAPRISM_short$RNAseq.Sample.ID)
tpm_samples <- unique(tpm_METAPRISM_long$sample)
tpm_METAPRISM_long_annotated <- left_join(tpm_METAPRISM_long, annotation_METAPRISM_short, join_by('sample' == 'RNAseq.Sample.ID'))

# Merge annotated tpm table into our df
METAPRISM_df_with_expr <- left_join(METAPRISM_df, tpm_METAPRISM_long_annotated, join_by('sample.id.tumor' == 'WES.Sample.ID', 'gene_name' == 'Gene.name'))
METAPRISM_df_with_expr$sample <- NULL
METAPRISM_df_with_expr$Sample.ID <- NULL
METAPRISM_df_with_expr$Sample.Sequencing <- NULL



# ADD NORMAL EXPRESSION FOR CORRESPONDING TISSUES

# Load GTEX prepared table to get normal tissue expression values
GTEX <- read.csv('GTEX/GTEX_normal_tissues_expr_summary.csv')


# Add information about primary site
MET500_df$Primary_Site_and_Subsite <- paste(MET500_df$Primary_Site, '-', MET500_df$Primary_Subsite)

unique(MET500_df$Primary_Site_and_Subsite)

# Match with the normal tissues expression
tissues <- c("Bladder - Bladder, NOS", "Breast - Breast, NOS", "Esophagus - Esophagus, NOS", "Colon - Colon, NOS", "Prostate gland - Prostate gland", "Ureter - Ureter", "Other and unspecified parts of biliary tract - Extrahepatic bile duct",  "Unknown - Unknown primary site", "Other and unspecified parts of mouth - Mouth, NOS", "Bronchus and lung - Lung, NOS", "Other and unspecified parts of tongue - Tongue, NOS", "Gallbladder - Gallbladder", "Brain - Brain, NOS", "Adrenal gland - Adrenal gland", "Connective, subcutaneous and other soft tissues - Connective tissue abdomen", "Retroperitoneum and peritoneum - Retroperitoneum", "Accessory sinuses - Accessory sinus, NOS", "Other and unspecified parts of biliary tract - Ampulla of Vater", "Testis - Testis, NOS", "Ovary - Ovary", "Connective, subcutaneous and other soft tissues - Connective, subcutaneous and other soft tissues, NOS", "Stomach - Stomach, NOS", "Bones, joints and articular cartilage - Upper limb long bones, joints")
matching_normal_tissue <- c("Bladder", "Breast", "Esophagus", "Colon", "Prostate", "Bladder", NA, NA, "Esophagus", "Lung", "Esophagus", NA, "Brain", "Adrenal Gland", "Adipose Tissue", "Adipose Tissue", NA, "Small Intestine", "Testis", "Ovary", NA, "Stomach", NA)
tissues_matching_MET500 <- data.frame(tissues, matching_normal_tissue)

MET500_df <- left_join(MET500_df, tissues_matching_MET500, 
                                     join_by('Primary_Site_and_Subsite' == 'tissues'))

MET500_df[MET500_df$Primary_Site_and_Subsite == 'Connective, subcutaneous and other soft tissues - Connective tissue leg' 
                        & MET500_df$Histological_Type == "Undifferentiated pleomorphic sarcoma", c(length(colnames(MET500_df)))] <- NA

MET500_df[MET500_df$Primary_Site_and_Subsite == 'Connective, subcutaneous and other soft tissues - Connective tissue leg' 
                        & MET500_df$Histological_Type == "Dedifferentiated liposarcoma", c(length(colnames(MET500_df)))] <- "Adipose Tissue"

MET500_df[MET500_df$Primary_Site_and_Subsite == 'Connective, subcutaneous and other soft tissues - Connective tissue leg' 
                        & MET500_df$Histological_Type == "Soft tissue myoepithelial carcinoma", c(length(colnames(MET500_df)))] <- NA


MET500_df <- left_join(MET500_df, GTEX, by = c('matching_normal_tissue' = 'X_primary_site', 'gene_name' = 'gene_name')) 


# Remove unnecessary columns
MET500_df$X <- NULL
MET500_df$median <- NULL
MET500_df$std <- NULL
MET500_df$n <- NULL

colnames(MET500_df)[length(colnames(MET500_df))] <- 'normal_tissue_mean_expr'

# In GTEX expression data are presented as log2(TPM + 0.001), in MET500 as raw counts => unify
MET500_df$normal_tissue_mean_expr_delog <- 2^(MET500_df$normal_tissue_mean_expr)
MET500_df$normal_tissue_mean_expr_delog_round <- round(MET500_df$normal_tissue_mean_expr_delog, digits = 6)

MET500_df$normal_tissue_mean_expr <- NULL
MET500_df$normal_tissue_mean_expr_delog <- NULL
colnames(MET500_df)[length(colnames(MET500_df))] <- 'normal_tissue_mean_expr'

write.csv(MET500_df, 'MET500/MET500_with_expr_and_normal_expr.csv', row.names = FALSE)

# METAPRISM
cancer_types <- sort(unique(METAPRISM_df_with_expr$Cancer.Type))
matching_normal_tissue <- c("Bladder", "Breast", NA, "Cervix Uteri", "Colon", "Esophagus", NA, NA, "Lung", "Ovary", "Pancreas", "Prostate", "Salivary Gland", "Blood Vessel")
tissues_matching_METAPRISM <- data.frame(cancer_types, matching_normal_tissue)

METAPRISM_df_with_expr <- left_join(METAPRISM_df_with_expr, tissues_matching_METAPRISM, 
                                     join_by('Cancer.Type' == 'cancer_types'))

unique(METAPRISM_df_with_expr$Tumor.Site[!is.na(METAPRISM_df_with_expr$Cancer.Type)
                                                   & METAPRISM_df_with_expr$Cancer.Type == 'Head and Neck Cancer'])


METAPRISM_df_with_expr[!is.na(METAPRISM_df_with_expr$Cancer.Type) & 
                         METAPRISM_df_with_expr$Cancer.Type == 'Head and Neck Cancer' &
                         METAPRISM_df_with_expr$Tumor.Site == 'Parotid gland', c(length(colnames(METAPRISM_df_with_expr)))] <- 'Salivary Gland'

METAPRISM_df_with_expr[!is.na(METAPRISM_df_with_expr$Cancer.Type) & 
                         METAPRISM_df_with_expr$Cancer.Type == 'Head and Neck Cancer' &
                         METAPRISM_df_with_expr$Tumor.Site == 'Tongue, NOS', c(length(colnames(METAPRISM_df_with_expr)))] <- 'Esophagus'

METAPRISM_df_with_expr[!is.na(METAPRISM_df_with_expr$Cancer.Type) & 
                         METAPRISM_df_with_expr$Cancer.Type == 'Head and Neck Cancer' &
                         METAPRISM_df_with_expr$Tumor.Site == 'Oropharynx, NOS', c(length(colnames(METAPRISM_df_with_expr)))] <- 'Esophagus'


METAPRISM_df_with_expr <- left_join(METAPRISM_df_with_expr, GTEX, by = c('matching_normal_tissue' = 'X_primary_site', 'gene_name' = 'gene_name')) 


# Remove unnecessary columns
METAPRISM_df_with_expr$X <- NULL
METAPRISM_df_with_expr$median <- NULL
METAPRISM_df_with_expr$std <- NULL
METAPRISM_df_with_expr$n <- NULL

colnames(METAPRISM_df_with_expr)[length(colnames(METAPRISM_df_with_expr))] <- 'normal_tissue_mean_expr'

# In GTEX expression data are presented as log2(TPM + 0.001), in MET500 as raw counts => unify
METAPRISM_df_with_expr$normal_tissue_mean_expr_delog <- 2^(METAPRISM_df_with_expr$normal_tissue_mean_expr)
METAPRISM_df_with_expr$normal_tissue_mean_expr_delog_round <- round(METAPRISM_df_with_expr$normal_tissue_mean_expr_delog, digits = 6)

METAPRISM_df_with_expr$normal_tissue_mean_expr <- NULL
METAPRISM_df_with_expr$normal_tissue_mean_expr_delog <- NULL
colnames(METAPRISM_df_with_expr)[length(colnames(METAPRISM_df_with_expr))] <- 'normal_tissue_mean_expr'

write.csv(METAPRISM_df_with_expr, 'METAPRISM/METAPRISM_with_expr_and_normal_expr.csv', row.names = FALSE)
