library(tidyr)
library(dplyr)
library(openxlsx)

setwd() #here write you information

# IMPORT DATASETS (keep information about oncogenes and TSG)

data <- read.csv("../../../Data/Processed_data/table_with_gene_expr.csv")
COSMIC <- read.csv("../../../Data/Raw_data/Cancer_Gene_census_COSMIC.tsv", sep = "\t")
oncoKB <- read.csv("../../../Data/Raw_data/cancerGeneList_OncoKB.tsv", sep = "\t")

# KEEP ONLY RELATED COLUMNS

COSMIC <- COSMIC %>%
  select("Gene.Symbol", "Tier", "Hallmark", "Role.in.Cancer", "Mutation.Types", 
         "Translocation.Partner")
colnames(COSMIC) <- paste(colnames(COSMIC), ".COSMIC", sep = "")
COSMIC <- COSMIC[grepl("TSG|oncogene", COSMIC$Role.in.Cancer), ]


oncoKB <- oncoKB %>%
  select("Hugo.Symbol", "Is.Oncogene", "Is.Tumor.Suppressor.Gene", "OncoKB.Annotated")
colnames(oncoKB) <- paste(colnames(oncoKB), ".oncoKB", sep = "")
oncoKB <- oncoKB[(oncoKB$Is.Tumor.Suppressor.Gene.oncoKB=="Yes"|oncoKB$Is.Oncogene=="Yes"), ]

# MERGE TABLES INTO OUR DATA

data_with_COSMIC <- left_join(data, COSMIC, join_by("gene_name" == "Gene.Symbol.COSMIC"))
data_with_COSMIC_oncoKB <- left_join(data_with_COSMIC, oncoKB, join_by("gene_name" == "Hugo.Symbol.oncoKB"))

data_with_COSMIC_oncoKB$Tier.COSMIC <- as.character(data_with_COSMIC_oncoKB$Tier.COSMIC)
colnames_to_replace_na <- c("Tier.COSMIC", "Hallmark.COSMIC", 
                            "Role.in.Cancer.COSMIC", "Mutation.Types.COSMIC",
                            "Translocation.Partner.COSMIC", "Is.Oncogene.oncoKB",
                            "Is.Tumor.Suppressor.Gene.oncoKB", "OncoKB.Annotated.oncoKB")   
data_with_COSMIC_oncoKB <- data_with_COSMIC_oncoKB %>%
  mutate_at(colnames_to_replace_na, ~replace_na(.x, "No_data"))

# EXPORT TABLE

write.csv(data_with_COSMIC_oncoKB, file = "../../../Data/Processed_data/Table_COSMIC_oncoKB_oncogenes_TSG.csv", row.names = FALSE)

