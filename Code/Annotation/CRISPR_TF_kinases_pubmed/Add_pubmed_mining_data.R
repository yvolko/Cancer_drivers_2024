library(dplyr)
library(tidyr)

# IMPORT DATASETS

data <- read.csv("../../../Data/Processed_data/Table_COSMIC_oncoKB_oncogenes_TSG.csv")
pubmed <- read.csv("../../../Data/Raw_data/pubmed_data_260224.txt", sep="\t", header = TRUE)

# MERGE PUBMED SERCH RESULTS INTO OUR DATA

data_pubmed <- left_join(data, pubmed, join_by("gene_name" == "gene"))

# EXPORT

write.csv(data_pubmed, file = "../../../Data/Processed_data/Table_with_pubmed.csv", row.names = FALSE)
