library(dplyr)
library(tidyr)
library(rvest)
library(tibble)
library(data.table)

# GTEX data summary by primary site


# Data are taken from here
# TPM:  https://xenabrowser.net/datapages/?dataset=gtex_RSEM_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# phenotype:   https://xenabrowser.net/datapages/?dataset=GTEX_phenotype&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Units: log2(tpm+0.001)

# Prepare GTEX data 
GTEX_pheno <- fread("GTEX_phenotype.gz")
GTEX_tpm <- fread("gtex_RSEM_gene_tpm.gz")
GTEX_tpm[, 2:7863] <-  round(GTEX_tpm[, 2:7863], digits = 2)
GTEX_gene_annot <- fread("probeMap_gencode.v23.annotation.gene.probemap")
GTEX_gene_annot <- GTEX_gene_annot[, c(1,2)]


GTEX_tpm_long <- GTEX_tpm %>% 
  pivot_longer(cols = -c(sample), names_to = 'sample_Id', values_to = 'expr')

GTEX_tpm_long_gene_names <- left_join(GTEX_tpm_long, GTEX_gene_annot, join_by('sample' == 'id'))
colnames(GTEX_tpm_long_gene_names)[1] <- 'gene_id'
colnames(GTEX_tpm_long_gene_names)[4] <- 'gene_name'

GTEX_tpm_long_gene_names <- left_join(GTEX_tpm_long_gene_names, GTEX_pheno, join_by('sample_Id' == 'Sample'))


# Summarize GTEX data: calculate mean, median and std for each primary site and each gene in GTEX
GTEX_summary_primary_site <- GTEX_tpm_long_gene_names %>%
  group_by(`_primary_site`, gene_name) %>%
  summarize(
    mean = mean(expr),
    median = median(expr),
    std = sd(expr),
    n = n()
  ) %>%
  ungroup()

GTEX_summary_primary_site <- GTEX_summary_primary_site[GTEX_summary_primary_site$`_primary_site` != '<not provided>', ]
GTEX_summary_primary_site[, 3:5] <- round(GTEX_summary_primary_site[3:5], digits = 2)

# Save resulting df
write.csv(GTEX_summary_primary_site, 'GTEX_normal_tissues_expr_summary.csv', row.names = FALSE)
