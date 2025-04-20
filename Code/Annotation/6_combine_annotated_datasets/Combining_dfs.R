library(dplyr)

# COMBINE TCGA, MET500 AND METAPRISM, HARTWIG ANNOTATED DATA TOGETHER  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

TCGA_df <- read.csv('../4_annotation_with_DBs/TCGA_and_Hartwig/TCGA_with_GO.csv')
MET500_df <- read.csv('../5_expression/MET500/MET500_with_expr_and_normal_expr.csv')
METAPRISM_df <- read.csv('../5_expression/METAPRISM/METAPRISM_with_expr_and_normal_expr.csv')
HARTWIG <- read.csv('../4_annotation_with_DBs/TCGA_and_Hartwig/HARTWIG_with_GO.csv')

# TCGA
TCGA_df$dataset <- 'TCGA'

TCGA_df_min <- TCGA_df %>%
  dplyr::select(dataset, sample, chr, startpos, endpos, nMajor, gene_name, Role.in.Cancer.COSMIC, 
                Is.Oncogene.oncoKB, Is.Tumor.Suppressor.Gene.oncoKB, Is.TF., Family_kinase, median.CRISPR,
                mean.CRISPR, sd.CRISPR, min.CRISPR, pubmed_invasion, pubmed_cancer, pubmed_growth, pubmed_migration,
                pubmed_prolifiration, pubmed_mean, PPI_count, expr, GO_terms, norm_expr)

colnames(TCGA_df_min)[6] <- 'copy_number'

# MET500
MET500_df$sample <- paste(MET500_df$sample.id.tumor, '_Vs_', MET500_df$sample.id.normal, sep='')
MET500_df$pubmed_mean <- (MET500_df$pubmed_invasion + MET500_df$pubmed_cancer + MET500_df$pubmed_growth +
                              MET500_df$pubmed_migration + MET500_df$pubmed_prolifiration) / 5

MET500_df_min <- MET500_df %>%
  dplyr::select(dataset, sample, chrom, start, end, tcn.em, gene_name, Role.in.Cancer.COSMIC, 
                Is.Oncogene.oncoKB, Is.Tumor.Suppressor.Gene.oncoKB, Is.TF., Family_kinase, median.CRISPR,
                mean.CRISPR, sd.CRISPR, min.CRISPR, pubmed_invasion, pubmed_cancer, pubmed_growth, pubmed_migration,
                pubmed_prolifiration, pubmed_mean, PPI_count, expr, GO_terms, normal_tissue_mean_expr)


colnames(MET500_df_min)[3] <- 'chr'
colnames(MET500_df_min)[4] <- 'startpos'
colnames(MET500_df_min)[5] <- 'endpos'
colnames(MET500_df_min)[6] <- 'copy_number'
colnames(MET500_df_min)[26] <- 'norm_expr'

# METAPRISM
METAPRISM_df$sample <- paste(METAPRISM_df$sample.id.tumor, '_Vs_', METAPRISM_df$sample.id.normal, sep='')
METAPRISM_df$pubmed_mean <- (METAPRISM_df$pubmed_invasion + METAPRISM_df$pubmed_cancer + METAPRISM_df$pubmed_growth +
                               METAPRISM_df$pubmed_migration + METAPRISM_df$pubmed_prolifiration) / 5

METAPRISM_df_min <- METAPRISM_df %>%
  dplyr::select(dataset, sample, chrom, start, end, tcn.em, gene_name, Role.in.Cancer.COSMIC, 
                Is.Oncogene.oncoKB, Is.Tumor.Suppressor.Gene.oncoKB, Is.TF., Family_kinase, median.CRISPR,
                mean.CRISPR, sd.CRISPR, min.CRISPR, pubmed_invasion, pubmed_cancer, pubmed_growth, pubmed_migration,
                pubmed_prolifiration, pubmed_mean, PPI_count, expr, GO_terms, normal_tissue_mean_expr)


colnames(METAPRISM_df_min)[3] <- 'chr'
colnames(METAPRISM_df_min)[4] <- 'startpos'
colnames(METAPRISM_df_min)[5] <- 'endpos'
colnames(METAPRISM_df_min)[6] <- 'copy_number'
colnames(METAPRISM_df_min)[26] <- 'norm_expr'

# HARTWIG
HARTWIG$dataset <- 'HARTWIG'
HARTWIG$pubmed_mean <- (HARTWIG$pubmed_cancer + HARTWIG$pubmed_growth + HARTWIG$pubmed_invasion +
                          HARTWIG$pubmed_migration + HARTWIG$pubmed_prolifiration) /5
HARTWIG_df_min <- HARTWIG %>%
  dplyr::select(dataset, sample, chr, start, end, n_copy, gene_name, Role.in.Cancer.COSMIC, 
  Is.Oncogene.oncoKB, Is.Tumor.Suppressor.Gene.oncoKB, Is.TF., Family_kinase, median.CRISPR,
  mean.CRISPR, sd.CRISPR, min.CRISPR, pubmed_invasion, pubmed_cancer, pubmed_growth, pubmed_migration,
  pubmed_prolifiration, pubmed_mean, PPI_count, expr, GO_terms, norm_expr, expr_fold)

colnames(HARTWIG_df_min)[4] <- 'startpos'
colnames(HARTWIG_df_min)[5] <- 'endpos'
colnames(HARTWIG_df_min)[6] <- 'copy_number'

# MAKE ALL THE EXPRESSION BE IN RAW TPM VALUES + ADD FOLD CHANGE
# GTEX HAS EXPRESSION AS log2(tpm+0.001) (if tpm is 0, than substractig 0.001 leads to negative value -> leave substraction of 0.001 out)
# I already calculated raw tpm for normal expression from GTEX in MET500 and METAPRISM DF before
# TCGA also has log2(tpm+0.001) expression values
# MET500 and METAPRISM have raw tpm values
TCGA_df_min$expr <- round(2^TCGA_df_min$expr, digits = 6)
TCGA_df_min$norm_expr <- round(2^TCGA_df_min$norm_expr, digits = 6)
TCGA_df_min$expr_fold <- round(TCGA_df_min$expr / TCGA_df_min$norm_expr, digits = 6)

MET500_df_min$expr_fold <- round(MET500_df_min$expr / MET500_df_min$norm_expr, digits = 6)
METAPRISM_df_min$expr_fold <- round(METAPRISM_df_min$expr / METAPRISM_df_min$norm_expr, digits = 6)


# COMBINE DATASETS
TCGA_HARTWIG_metastatic <- rbind(TCGA_df_min, MET500_df_min, METAPRISM_df_min, HARTWIG_df_min)
write.csv(TCGA_HARTWIG_metastatic, 'TCGA_MET500_METAPRISM_HARTWIG.csv', row.names = FALSE)
