library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(Hmisc)


download.file('https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz', '../../Data/Raw_data/tcga_RSEM_gene_tpm.gz', method="auto", mode="wb")

# DATA LOADING
tpm_data <- fread("../../Data/Raw_data/tcga_RSEM_gene_tpm.gz")
round_tpm_data <- round(tpm_data[,2:10536], digits = 2)
round_tpm_data$sample <- tpm_data$sample
annotation_probeMap <- read.table('../../Data/Raw_data/probeMap_gencode.v23.annotation.gene.probemap', sep = '\t', header = TRUE)
table_with_gene_name <- read.csv('../../Data/Processed_data/Table_with_gene_names_from_gencode.csv', sep = ',')
phenotype_anot <- read.table('../../Data/Raw_data/TCGA_phenotype_denseDataOnlyDownload.tsv', sep = '\t', header = TRUE)


#DATA PREPARING
gene_list <- which(annotation_probeMap$gene %in% table_with_gene_name$gene_name)
annotation_gene_list <- annotation_probeMap[gene_list,]
ens_gene_list <- which(round_tpm_data$sample %in% annotation_gene_list$id)
table_gene_tpm_data <- round_tpm_data[ens_gene_list,]
annotation_gene_list <- annotation_gene_list %>%
  rename(sample = id)
table_gene_tpm_data <- table_gene_tpm_data %>%
  left_join(annotation_gene_list, by = 'sample') 
table_gene_tpm_data <- table_gene_tpm_data[,1:10537]
sample_anot <- phenotype_anot[,-2]
filtered_table_gene_tpm_data <- colnames(table_gene_tpm_data)[grep("-01$", colnames(table_gene_tpm_data))]
table_gene_tpm_data_filt <- table_gene_tpm_data[,..filtered_table_gene_tpm_data]
table_gene_tpm_data_filt$gene_name <- table_gene_tpm_data$gene
table_gene_tpm_data_filt$sample <- table_gene_tpm_data$sample
table_gene_tpm_data_long <- table_gene_tpm_data_filt %>%
  pivot_longer(cols = -c(gene_name, sample), names_to = 'sample_cut', values_to = 'expr')
table_gene_tpm_data_long <- left_join(table_gene_tpm_data_long, sample_anot, join_by('sample_cut' == 'sample'))
table_with_gene_name$sample_cut <- str_sub(table_with_gene_name$sample, 1, 12)
table_gene_tpm_data_long$sample_cut <- str_sub(table_gene_tpm_data_long$sample_cut, 1, 12)
result_df <- left_join(table_with_gene_name, table_gene_tpm_data_long, by = c('sample_cut', 'gene_name'))
colnames(result_df)[colnames(result_df) == "sample.y"] <- "ENS_id"
colnames(result_df)[colnames(result_df) == "sample.x"] <- "sample"


# DATA MERGING
table_gene_tpm_data <- round_tpm_data[ens_gene_list,]
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
result_df_with_norm <- result_df_with_norm[,-c(17,18,19)]
result_df_with_norm <- result_df_with_norm %>%
  group_by(sample_cut, gene_name) %>%
  filter(expr == max(expr) | is.na(expr))

for (i in 1:48860){
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

result_df_with_norm <- result_df_with_norm[,-c(17,18,19)]
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