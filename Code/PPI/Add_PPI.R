#install.packages("rbioapi")
library(rbioapi)
library(dplyr)


df <- read.csv('../../Data/Processed_data/Table_expression_COSMIC_oncoKB_merged_intervals.csv')
int_net_1 <- rba_string_interactions_network(ids = unique(df$gene_name[1:2000]),
                                           species = 9606,
                                           required_score = 700)
int_net_2 <- rba_string_interactions_network(ids = unique(df$gene_name[2000:4000]),
                                             species = 9606,
                                             required_score = 700)
int_net_3 <- rba_string_interactions_network(ids = unique(df$gene_name[4000:6000]),
                                             species = 9606,
                                             required_score = 700)
int_net_4 <- rba_string_interactions_network(ids = unique(df$gene_name[6000:8000]),
                                             species = 9606,
                                             required_score = 700)
int_net_5 <- rba_string_interactions_network(ids = unique(df$gene_name[8000:10000]),
                                             species = 9606,
                                             required_score = 700)
int_net_6 <- rba_string_interactions_network(ids = unique(df$gene_name[10000:11475]),
                                             species = 9606,
                                             required_score = 700)
summary_int_net <- rbind(int_net_1, int_net_2, int_net_3, int_net_4, int_net_5, int_net_6)

my_gene_list <- data.frame()
for (i in unique(df$gene_name)){
  if (i %in% summary_int_net$preferredName_B){
    my_sum <- 0
    my_gene_name <- i
    for (j in unique(filter(summary_int_net, preferredName_B == i)$preferredName_A)){
      if (j %in% unique(df$gene_name)){
        new_list <- filter(df, gene_name == j)
        if (('oncogene' %in% new_list$Role.in.Cancer.COSMIC) | ('Yes' %in% new_list$Is.Oncogene.oncoKB)){
          my_sum <- my_sum + 1
        }
      }
    }
    my_gene_list <- rbind(my_gene_list, c(my_gene_name, my_sum))
  }
}


for (i in unique(df$gene_name)){
  if (i %in% summary_int_net$preferredName_A){
    my_sum <- 0
    my_gene_name <- i
    for (j in unique(filter(summary_int_net, preferredName_A == i)$preferredName_B)){
      if (j %in% unique(df$gene_name)){
        new_list <- filter(df, gene_name == j)
        if (('oncogene' %in% new_list$Role.in.Cancer.COSMIC) | ('Yes' %in% new_list$Is.Oncogene.oncoKB)){
          my_sum <- my_sum + 1
        }
      }
    }
    my_gene_list <- rbind(my_gene_list, c(my_gene_name, my_sum))
  }
}

colnames(my_gene_list) <- c('gene', 'PPI_count')
my_gene_list$PPI_count <- as.numeric(my_gene_list$PPI_count)
summarized_gene_list <- aggregate(. ~ gene, data = my_gene_list, FUN = sum)
write.table(summarized_gene_list, '../../Data/Processed_data/PPI_list.txt', quote = FALSE, row.names = FALSE)
