library(dplyr)
library(tidyr)
library(rvest)
library(tibble)

# IMPORT DATA

data <- read.csv("../../../Data/Processed_data/Table_with_pubmed.csv")

# ADD INFORMATION ABOUT TRANSCRIPTIONAL FACTORS

TF <- read.csv("../../../Data/Raw_data/Human_TFs.txt", sep="\t")
TF <- TF %>% 
  select("HGNC.symbol", "Is.TF.", "DBD")
TF <- TF[TF$Is.TF.=="Yes", ]

data_TF <- left_join(data, TF, join_by("gene_name" == "HGNC.symbol"))
data_TF <- data_TF %>%
  mutate_at(c("Is.TF.", "DBD"), ~replace_na(.x, "No_data"))

# ADD INFORMATION ABOUT KINASES

webpage <- read_html("http://www.kinhub.org/kinases.html")
kinase <- html_table(webpage)[[1]]
kinase <- kinase %>%
  select("xName", "HGNC Name", "Group", "Family")
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
# Please download file from this link https://depmap.org/portal/download/all/?releasename=DepMap+Public+21Q3&filename=CRISPR_gene_effect.csv
# File is too big for github
CRISPR <- read.csv("CRISPRGeneEffect.csv")
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

# EXPORT FINAL DATAFRAME

write.csv(data_TF_kinase_CRISPR, file = "../../../Data/Processed_data/Table_TF_kinase_CRISPR.csv", row.names = FALSE)
