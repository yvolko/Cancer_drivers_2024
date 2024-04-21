library(dplyr)
library(tidyr)

# IMPORT DATASETS

data <- read.csv("../result_files/Table_expr_merged_TF_kinase_CRISPR_min.csv")
pubmed <- read.csv("./pubmed_data_260224.txt", sep="\t", header = TRUE)

data <- select(data, -pubmed_cancer, -pubmed_growth, -pubmed_proliferation)

# MERGE PUBMED SERCH RESULTS INTO OUR DATA

data_pubmed <- left_join(data, pubmed, join_by("gene_name" == "gene"))

# EXPORT

write.csv(data_pubmed, file = "Table_version_240228.csv", row.names = FALSE)
