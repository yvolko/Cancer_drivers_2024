library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# IMPORT DATA

data <- read.csv("../../Data/Processed_data/Merge_table_with_gene_names_from_gencode.csv")

# SELECT ALL UNIQUE REGIONS FROM OUR DATA

data <- data %>% 
  dplyr::select("chr", "startpos", "endpos")
data$chr <- paste("chr", data$chr, sep="")
data <- unique(data)

# EXPORT SELECTED REGIONS TO BED FILE

write.table(data, file = "regions_from_data.bed", sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# MAKE FILE WITH CHROMOSOME LENGTHS (can be taken from here https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&chromInfoPage=)

chromosome_lengths <- data.frame(
  chr = paste0("chr", 1:22),  
  length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 
             159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 
             115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 
             59128983, 63025520, 48129895, 51304566)  
)

chromosome_lengths <- rbind(chromosome_lengths, c("chrX", 155270560))
chromosome_lengths <- rbind(chromosome_lengths, c("chrY", 59373566))

write.table(chromosome_lengths, file = "chromosome_lengths.bed", sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# RUN THIS COMMAND IN TERMINAL: bedtools shuffle -i regions_from_data.bed -g chromosome_lengths.bed > random_full_genome_output.bed 
# THIS WILL GET RANDOM INTERVALS OF THE SAME LENGTH ACROSS HUMAN GENOME
# IMPORT RESULT FILE

random_intervals <- read.csv("../../Data/Raw_data/random_full_genome_output.bed", sep="\t", header = FALSE, col.names = c("chr", "start", "end"))

# CHECK IF INTERVAL LENGTH IN RANDOM INTERVALS IS THE SAME AS IN OUR DATA
data$length <- data$endpos-data$startpos
random_intervals$length <- random_intervals$end-random_intervals$start

all(sort(random_intervals$length) == sort(data$length))

# ANNOTATE RANDOM INTERVALS WITH GENE NAMES
# Fist prepared genome annotation in a way it was done for initial data annotation

# Use Gencode V19 GTF file for annotation
# Download from here: https://www.gencodegenes.org/human/release_19.html
# Read GTF file

gtf_new <- readGFF("gencode.v19.annotation.gtf")
gtf <- gtf_new[gtf_new$gene_type == "protein_coding" & gtf_new$type == "gene", ]


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

# MAKE GENOMIC RANGES FROM GTF
gr_gencode <- GRanges(
  seqnames = gtf$seqid,
  ranges = IRanges(start = gtf$start, end = gtf$end)
)

# MAKE GENOMIC RANGES FROM RANDOM INTERVALS
gr_random_intervals <- GRanges(
  seqnames = random_intervals$chr,
  ranges = IRanges(start = random_intervals$start, end = random_intervals$end)
)

# OVERLAP GTF AND RANDOM INTERVALS
overlapping_genes_gencode <- findOverlaps(gr_random_intervals, gr_gencode, type = "any", ignore.strand=TRUE)

# Getting subset of gtf table by using indexces from findOverlaps Hits
gtf <- gtf %>% mutate(Index = row_number())
random_intervals <- random_intervals %>% mutate(Index = row_number())
gtf_subset <- gtf[subjectHits(overlapping_genes_gencode), ]
gtf_subset$data_indexces <- queryHits(overlapping_genes_gencode)
columns_to_keep <- c("gene_name", "Index", "data_indexces")
gtf_subset <- gtf_subset[, columns_to_keep]

# Merge data from gtf into random intervals
random_intervals <- left_join(random_intervals, gtf_subset, join_by("Index" == "data_indexces"))
random_intervals$Index <- NULL
random_intervals$Index.y <- NULL

# ANALYSE
# WHAT PERCENTAGE OF INTERVALS IN OUR DATA HAVE KNOWN ONCOGENES

COSMIC <- read.csv("../../Data/Raw_data/Cancer_Gene_census_COSMIC.tsv", sep = "\t")
oncoKB <- read.csv("../../Data/Raw_data/cancerGeneList_OncoKB.tsv", sep = "\t")

COSMIC <- COSMIC %>%
  dplyr::select("Gene.Symbol", "Role.in.Cancer")
colnames(COSMIC) <- paste(colnames(COSMIC), ".COSMIC", sep = "")

oncoKB <- oncoKB %>%
  dplyr::select("Hugo.Symbol", "Is.Oncogene")
colnames(oncoKB) <- paste(colnames(oncoKB), ".oncoKB", sep = "")

our_data <- read.csv("../../Data/Processed_data/Merge_table_with_gene_names_from_gencode.csv")

our_data <- left_join(our_data, COSMIC, by=c("gene_name" = "Gene.Symbol.COSMIC"))
our_data <- left_join(our_data, oncoKB, by=c("gene_name" = "Hugo.Symbol.oncoKB"))


our_data <- our_data %>% dplyr::select("chr", "startpos", "endpos", "Role.in.Cancer.COSMIC", "Is.Oncogene.oncoKB")
our_data$coordinates <- paste(as.character(our_data$chr), ":", our_data$startpos, "-", our_data$endpos)
our_data$chr <- NULL
our_data$startpos <- NULL
our_data$endpos <- NULL

our_data <- our_data %>%
  mutate_at(c("Role.in.Cancer.COSMIC", "Is.Oncogene.oncoKB"), ~replace_na(.x, "No_data"))

regions_with_without_oncogenes <- our_data %>%
  group_by(coordinates) %>%
  mutate(
    oncogene_presence = ifelse(
      Role.in.Cancer.COSMIC %in% c("oncogene, fusion", "oncogene") | 
        Is.Oncogene.oncoKB == "Yes", 
      1, 
      0
    )
  )  %>%
  ungroup()

regions_with_without_oncogenes$Role.in.Cancer.COSMIC <- NULL
regions_with_without_oncogenes$Is.Oncogene.oncoKB <- NULL

regions_with_without_oncogenes <- unique(regions_with_without_oncogenes)

length(regions_with_without_oncogenes$coordinates) # 7777 unique regions
sum(regions_with_without_oncogenes$oncogene_presence) # 2126 regions have oncogenes (27.3% regions have oncogene)


# WHAT PERCENTAGE OF INTERVALS IN RANDOM INTERVALS HAVE KNOWN ONCOGENES

random_intervals <- left_join(random_intervals, COSMIC, by=c("gene_name" = "Gene.Symbol.COSMIC"))
random_intervals <- left_join(random_intervals, oncoKB, by=c("gene_name" = "Hugo.Symbol.oncoKB"))

random_intervals <- random_intervals %>%
  mutate_at(c("Role.in.Cancer.COSMIC", "Is.Oncogene.oncoKB"), ~replace_na(.x, "No_data"))
random_intervals$coordinates <- paste(as.character(random_intervals$chr), ":", random_intervals$start, "-", random_intervals$end)
random_intervals <- random_intervals %>% dplyr::select("Role.in.Cancer.COSMIC", "Is.Oncogene.oncoKB", "coordinates")

random_intervals_summary <- random_intervals %>%
  group_by(coordinates) %>%
  mutate(
    oncogene_presence = ifelse(
      Role.in.Cancer.COSMIC %in% c("oncogene, fusion", "oncogene") | 
        Is.Oncogene.oncoKB == "Yes", 
      1, 
      0
    )
  )  %>%
  ungroup()

random_intervals_summary$Role.in.Cancer.COSMIC <- NULL
random_intervals_summary$Is.Oncogene.oncoKB <- NULL
random_intervals_summary <- unique(random_intervals_summary)

length(random_intervals_summary$coordinates) # 6465 unique intervals
sum(random_intervals_summary$oncogene_presence) # 761 regions have oncogenes (11.8%)

# So, 27.3% in our data vs 11.8% in random intervals - % of regions with oncogenes

onco_data <- 2126
total_minus_onco_data <- 7777 - 2126

onco_random <- 761
total_minus_onco_random <- 6465 - 761

contingency_table <- matrix(c(onco_data, onco_random, total_minus_onco_data,
                              total_minus_onco_random), 
                            nrow = 2)

colnames(contingency_table) <- c("With known oncogenes", "Without known oncogenes")
rownames(contingency_table) <- c("High-level amplification regions", "Random intervals")

fisher.test(contingency_table)

# VISUALIZE
slice_colors <- c("#ED9564", "#6495ED")
percentages_data <- contingency_table[1,] / sum(contingency_table[1,]) * 100
percentages_random <- contingency_table[2,] / sum(contingency_table[2,]) * 100


barplot(t(contingency_table), beside = TRUE, col = c("#fed9a5", "#b3cee3"), legend = TRUE, ylab = "Regions count")

# Conclusion: There is an enrichment by oncogenes in our data regions vs random intervals 
