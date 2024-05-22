library(org.Hs.eg.db)
library(GO.db)
library(dplyr)

# DATA PREPARE

a <- read.csv("../../Data/Processed_data/Table_TF_kinase_CRISPR.csv")
df <- a %>%
  mutate(region_length = endpos - startpos,
         pubmed_mean = rowMeans(dplyr::select(., c("pubmed_growth", "pubmed_cancer", "pubmed_prolifiration", 
                                                   "pubmed_invasion", "pubmed_migration"))))
df$chr_conectied=paste(df$sample.x,df$chr,df$startpos,df$endpos,sep="_")
df_TSG <- df %>%
  group_by(chr_conectied) %>%
  dplyr::filter('TSG' %in% Role.in.Cancer.COSMIC & 'oncogene' %in% Role.in.Cancer.COSMIC)
TSG_oncogene_region_count <- length(unique(df_TSG$chr_conectied))

b=df[,colnames(df) %in% c("sample.x","chr","startpos","endpos","nMajor","gene_name","expr","norm_expr","Role.in.Cancer.COSMIC","Is.Oncogene.oncoKB","pubmed_mean","Is.TF.","Is.Tumor.Suppressor.Gene.oncoKB","Family_kinase","mean.CRISPR","min.CRISPR")]
b=b[complete.cases(b$norm_expr),] 

# REMOVE TUMOR SUPPRESSOR GENES

b <- filter(b, Role.in.Cancer.COSMIC != 'TSG' &
              Role.in.Cancer.COSMIC != 'TSG, fusion' &
              Is.Tumor.Suppressor.Gene.oncoKB != 'Yes')

b$chr=paste(b$sample.x,b$chr,b$startpos,b$endpos,sep="_")
b$fold=round((2^b$expr)/(2^b$norm_expr),digits = 2)

# CORRECTION OF DATA ON ONCOGENES

for (i in 1:length(b$chr)){
  if(b[i,"Role.in.Cancer.COSMIC"] == 'oncogene, TSG, fusion' & b[i, "Is.Oncogene.oncoKB"] != 'Yes'){
    b[i,"Role.in.Cancer.COSMIC"] = 'oncogene'
  }
}

for (i in 1:length(b$chr)){
  if(b[i,"Role.in.Cancer.COSMIC"] == 'oncogene, TSG' & b[i, "Is.Oncogene.oncoKB"] != 'Yes'){
    b[i,"Role.in.Cancer.COSMIC"] = 'oncogene'
  }
}
duplicated_oncogene <- dplyr::filter(b, Role.in.Cancer.COSMIC == 'oncogene' | Is.Oncogene.oncoKB == 'Yes')$gene_name

# ADD GO RESULTS

cluster_summary_onco <- read.table('../../Data/Processed_data/cluster_summary_oncogene.txt', header = TRUE, sep = ' ')
new_feature <- cluster_summary_onco[,c(3,9)]
for (i in cluster_summary_onco$Description){
  go_id = GOID( GOTERM[ Term(GOTERM) == i])
  get(go_id, org.Hs.egGO2ALLEGS)
  allegs = get(go_id, org.Hs.egGO2ALLEGS)
  genes = unlist(mget(allegs,org.Hs.egSYMBOL), use.names = FALSE)
  new_feature[new_feature$Description == i, 'geneID'] <- paste(unique(genes), sep = ',', collapse = '/')
}

new_feature <- new_feature %>%
  mutate(gene_list = strsplit(geneID, "/"))

GO_df <- data.frame()
for (i in unique(df$gene_name)){
  gene = i
  GO_terms = sum(sapply(new_feature$gene_list, function(x) i %in% x))
  GO_df <- rbind(GO_df, c(gene, GO_terms))
}

colnames(GO_df) <- c('gene', 'GO_terms')
b <- left_join(b, GO_df, join_by(gene_name == gene))

# ADD PPI RESULTS

summarized_gene_list <- read.table('../../Data/Processed_data/PPI_list.txt', header = TRUE, sep = ' ')
b <- left_join(b, summarized_gene_list, join_by(gene_name == gene))
b$PPI_count <- ifelse(is.na(b$PPI_count), 0, b$PPI_count)


summary_data <- data.frame()
for (i in unique(b$chr)){
  c=b[b$chr == i,]
  c$is_kinase=ifelse(c$Family_kinase != "No_data","Yes","No")
  c$is_TS=ifelse(c$Is.Tumor.Suppressor.Gene.oncoKB == "Yes","Da","No")
  c$expr2=rank(c$expr,na.last = F)
  c$pbm=rank(c$pubmed_mean,na.last = F)
  c$tf=rank(c$Is.TF.,na.last = F)
  c$kin=rank(c$is_kinase,na.last = F)
  c$crisp1=rank(1-c$mean.CRISPR,na.last = F)
  c$crisp2=rank(1-c$min.CRISPR,na.last = F)
  c$fold2=rank(c$fold,na.last = F)
  c$TS = rank(c$is_TS, na.last = F)
  c$GO = rank(c$GO_terms, na.last = F)
  c$PPI = rank(c$PPI_count, na.last = F)
  c$rank=c$expr2+c$pbm+c$tf+c$kin+c$crisp1+c$crisp2+c$fold2 + c$TS + c$GO + c$PPI 
  summary_data <- rbind(summary_data, c)
}


# FILTERING DUPLICATE ONCOGENES

new_summary_data <- data.frame()
for (i in unique(summary_data$chr)){
  sample_df <- dplyr::filter(summary_data, chr == i)
  
  count_1 <- sum(sample_df$Role.in.Cancer.COSMIC == 'oncogene')
  count_2 <- sum(sample_df$Is.Oncogene.oncoKB == 'Yes')
  
  if (count_1 <= 1 && count_2 <= 1){
    new_summary_data <- rbind(new_summary_data, sample_df)
  }
}
set.seed(123)
duplicated_oncogene <- dplyr::filter(new_summary_data, Role.in.Cancer.COSMIC == 'oncogene' | Is.Oncogene.oncoKB == 'Yes')$gene_name
duplicates <- duplicated_oncogene[duplicated(duplicated_oncogene)]

for (gene in duplicates) {
  duplicated_regions <- dplyr::filter(new_summary_data, gene_name == gene)$chr
  if (length(duplicated_regions) > 1) {
    sample_data <- dplyr::filter(new_summary_data, chr %in% duplicated_regions)
    random_region_to_keep <- sample(duplicated_regions, 1)
    new_summary_data <- new_summary_data[!(new_summary_data$chr %in% duplicated_regions), ]
    new_summary_data <- rbind(new_summary_data, sample_data[sample_data$chr == random_region_to_keep,])
  }
}

write.csv(new_summary_data, file = "../../Data/Processed_data/Table_with_rank_filtering_by_oncogene.csv", row.names = FALSE)

# SELECT GENES WITH A RANK HIGHER THAN DECILE 9

summary_data_not_oncogene <- data.frame()
for (i in unique(new_summary_data$chr)){
  sample_data <- dplyr::filter(new_summary_data, chr == i)
  sample_data <- sample_data %>% arrange(desc(rank))
  if (!any(grepl("oncogene", sample_data$Role.in.Cancer.COSMIC))){
    if(!('Yes' %in% sample_data$Is.Oncogene.oncoKB)){
      top_rank <- quantile(sample_data$rank, probs = 0.9)
      new_data <- dplyr::filter(sample_data, rank >= top_rank)
      summary_data_not_oncogene <- rbind(summary_data_not_oncogene, new_data)
    }
  }
}  

write.table(unique(summary_data_not_oncogene$gene_name), '../../Data/Result_data/Oncogenes_from_ranks.txt', quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

