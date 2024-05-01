#install.packages("openxlsx")
#install.packages("BiocManager")
#install.packages("RMariaDB")
#install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("org.Hs.eg.db")


library(openxlsx)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db) # Homo Sapiens annotation with gene names
library(rtracklayer) # Needed to read GTF file

setwd() #here write you information

# Raw data, page 2, csv file
filtered_data <- read.csv("../../Data/Processed_data/Table_filtering_by_nMajor_and_length.csv")

# Use Gencode V19 GTF file for annotation
# Download from here: https://www.gencodegenes.org/human/release_19.html
# Read GTF file
gtf <- readGFF("gencode.v19.annotation.gtf")
gtf <- gtf[gtf$gene_type == "protein_coding" & gtf$type == "gene", ]

# REMOVE DUPLICATES FROM ANNOTATION FILE
# check for duplicated names:
# length(names(which(table(gtf$gene_name) > 1))) -> 103
# They are very different in the sense of why they are duplicates:
# There are genes that are both on X and Y, or some actual duplication. 
duplicated <- names(which(table(gtf$gene_name) > 1))
potential_duplicates <- gtf[gtf$gene_name %in% duplicated, ]
# Let's remove chrX and chrY (there is actually no chrY in our data)
potential_duplicates <- potential_duplicates[!(potential_duplicates$seqid == "chrY" | potential_duplicates$seqid == "chrX"), ]

# Now we are going to choose the things we want to remove
# If one of the duplicate is "KNOWN" and another is "PUTATIVE", let's remove putative
# If one is "NOVEL" and another is "KNOWN", let's remove novel
selected_lines_1 <- potential_duplicates %>%
  group_by(gene_name) %>%
  filter(gene_status == "PUTATIVE") %>%
  ungroup()
selected_lines_2 <- potential_duplicates %>%
  group_by(gene_name) %>%
  filter(gene_status == "NOVEL") %>%
  ungroup()
selected_lines_3 <- potential_duplicates %>%
  group_by(gene_name) %>%
  filter(gene_status == "KNOWN") %>%
  ungroup()

# We can see that there are still some dulicates where both of them have "KNOWN" gene_status
# names(which(table(selected_lines_3$gene_name) > 1))
# Let remove the ones with smaller region
duplicated <- names(which(table(selected_lines_3$gene_name) > 1))
selected_lines_3 <- selected_lines_3[selected_lines_3$gene_name %in% duplicated, ]
selected_lines_3$length <- selected_lines_3$end-selected_lines_3$start
selected_lines_3 <- selected_lines_3 %>%
  group_by(gene_name) %>%
  slice_min(order_by = length) %>%
  ungroup()
selected_lines_3$length <- NULL

things_to_remove <- rbind(selected_lines_1, selected_lines_2, selected_lines_3)

# Now let's remove these duplicates from our gtf
# There is one more true duplicate on X chromosome (gene: IDS). Remove it too
gtf <- anti_join(gtf, things_to_remove, by = "gene_id")
gtf <- subset(gtf, subset = !(rownames(gtf) %in% c("20060")))

# GET OVERLAP BETWEEN OUR DATA AND ANNOTATION FILE
# Make genomic ranges from gtf
gr_gencode <- GRanges(
  seqnames = gtf$seqid,
  ranges = IRanges(start = gtf$start, end = gtf$end)
)

# Make genomic ranges from your data
gr_data <- GRanges(
  seqnames = filtered_data$chr,
  ranges = IRanges(start = filtered_data$startpos, end = filtered_data$endpos)
)

# Overlap gtf and your data
overlapping_genes_gencode <- findOverlaps(gr_data, gr_gencode, type = "any", ignore.strand=TRUE)

# Getting subset of gtf table by using indexces from findOverlaps Hits
gtf <- gtf %>% mutate(Index = row_number())
filtered_data <- filtered_data %>% mutate(Index = row_number())
gtf_subset <- gtf[subjectHits(overlapping_genes_gencode), ]
gtf_subset$data_indexces <- queryHits(overlapping_genes_gencode)
columns_to_keep <- c("gene_name", "Index", "data_indexces")
gtf_subset <- gtf_subset[, columns_to_keep]

# Merge data from gtf into your data
extended_filtered_data <- left_join(filtered_data, gtf_subset, join_by("Index" == "data_indexces"))
extended_filtered_data$Index <- NULL
extended_filtered_data$Index.y <- NULL

# Evaluate how many genes are without the name 
# sum(extended_filtered_data[["gene_name"]] == "" | is.na(extended_filtered_data[["gene_name"]])) -> 607 (~1.2%)
# check na per chromosome:
check_na <- extended_filtered_data %>%
  group_by(chr) %>%
  summarize(na_count = sum(is.na(gene_name)))

# EXPORT ANNOTATED DATA AND GENE LIST
# Prepare and export table
extended_filtered_data$chr <- sub("^chr", "", extended_filtered_data$chr)
extended_filtered_data$gene_name[is.na(extended_filtered_data$gene_name)] <- ""
write.csv(extended_filtered_data, file = "Table_with_gene_names_from_gencode.csv", row.names = FALSE)


# Extract gene names for esearch (search in PubMed in the contect of cancer etc.)
writeLines(unique(extended_filtered_data$gene_name), "selected_genes.txt") 
