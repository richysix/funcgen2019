# This script takes an expression file like the example ones for the 
# Zebrafish Functional Genomics Course that contains
# GeneIds, Adjusted pvalues, Log2 Fold Changes and Normalised Counts for samples

# load the tidyverse
library(tidyverse)
# load data
# you will need to change this if you are using a different sig file
sig_data_file <- 'uninf_3dpf_hom_vs_sib.sig.tsv'
sig_data <- read_tsv(sig_data_file)

# select normalised count columns
norm_counts <- 
  select(sig_data, contains("normalised"))

# calculate correlation matrix
# the cor function works on columns so we have to transpose the matrix with the t() function
cor_matrix <- cor(t(norm_counts))

# add gene ids and make a data frame
cor_df <- data.frame(
  GeneID = sig_data$GeneID,
  as.data.frame(cor_matrix)
  )
# add column and row names
colnames(cor_matrix) <- sig_data$GeneID
rownames(cor_matrix) <- sig_data$GeneID

# This matrix has reciprocal correlations for each gene pair
# we only want one edge per gene pair so
# get one half of matrix
num_edges <- nrow(sig_data) * (nrow(sig_data) - 1)/2
gene1  <- character(length = num_edges)
gene2  <- character(length = num_edges)
cor_values  <- numeric(length = num_edges)
index <- 1
for(i in seq_len(nrow(sig_data))) {
  for(j in seq_len(nrow(sig_data))) {
    if(i > j){
      gene1[index] <- rownames(cor_matrix)[i]
      gene2[index] <- rownames(cor_matrix)[j]
      cor_values[index] <- cor_matrix[i, j]
      index <- index + 1
     }
  }
}

# this is now the edges in the source -> target format
cor_long <- data.frame(
  Gene1ID = gene1,
  type = "geneExprCor",
  Gene2ID = gene2,
  PearsonCor = cor_values
)

# reorder columns and filter for Pearson > 0.8
# Change this number to change the filtering
cor_long <- filter(cor_long, PearsonCor > 0.8)

# reorder the columns for the sif format
network_df <- select(cor_long, Gene1ID, type, Gene2ID)

# write out network file
# change this filename if necessary
network_file <- 'uninf_3dpf_hom_vs_uninf_3dpf_sib.sif'
write.table(network_df, file = network_file, 
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# get node info
node_info <- select(sig_data, Name = GeneID, GeneName = Name, 
                    adjp = ends_with('adjp'), 
                    log2fc = ends_with('log2fc') )

# write out nodes file
# change this filename if necessary
nodes_file <- 'uninf_3dpf_hom_vs_uninf_3dpf_sib.nodes.txt'
write.table(node_info, file = nodes_file, 
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

# make edges df
edges_df <- data.frame(
  edge = paste0(cor_long$Gene1ID, ' (', cor_long$type, ') ', cor_long$Gene2ID),
  PearsonCor = cor_long$PearsonCor
)

# write out edges file
# change this filename if necessary
edges_file <- 'uninf_3dpf_hom_vs_uninf_3dpf_sib.edges.txt'
write.table(edges_df, file = edges_file, 
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

