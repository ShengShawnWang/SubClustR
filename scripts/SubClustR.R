#!/usr/bin/env Rscript
library(optparse)
library(stringr)

# ------------------------------
# Command-line options
# ------------------------------
option_list <- list(
  make_option(c("-d", "--expression_data"), type="character",
              help="Input CSV file with normalized gene expression values (e.g., microarray RMA or log2-transformed RNA-seq TPM/RPKM)", 
              metavar="character"),
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="Optional comma-separated list of genes of interest, or 'all' for all genes. Default: all genes", 
              metavar="character"),
  make_option(c("-p", "--patience"), type="numeric", default=25,
              help="Termination patience for FDR cutoff estimation [default = 25]"),
  make_option(c("-s", "--stabilization_threshold"), type="numeric", default=0.01,
              help="Threshold to stop auto-shuffle convergence [default = 0.01]"),
  make_option(c("-f", "--fdr"), type="numeric", default=0.01,
              help="False discovery rate for cutoff [default = 0.01]"),
  make_option(c("-o", "--out"), type="character",
              help="Output directory to save cluster CSVs", metavar="character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ------------------------------
# Helper Functions
# ------------------------------

# Flatten distance matrix excluding zeros
mat2dist <- function(input_matrix){  
  d_mat <- dist(input_matrix, method = "euclidean")
  d_list <- as.vector(as.matrix(d_mat))
  zero_count <- sum(d_list == 0)
  d_list <- d_list[d_list != 0]
  c(d_list, rep(0, zero_count - nrow(input_matrix)))
}

# Hierarchical clustering cut
cut_dend_tree <- function(data_frame, cutoff){
  hc <- hclust(dist(data_frame, method = "euclidean"), method = "complete")
  cutree(hc, h = cutoff)
}

# Cut into 2 clusters
cuttree_k2 <- function(data_frame){
  hc <- hclust(dist(data_frame, method = "euclidean"), method = "complete")
  cutree(hc, k = 2)
}

# Percentile of correlation matrix excluding diagonal
cor_matrix_percentile <- function(cor_matrix, percentile){
  cor_values <- as.vector(as.matrix(cor_matrix))
  cor_values <- cor_values[cor_values < 1]  # remove self-correlation
  quantile(cor_values, percentile)
}

# Auto-shuffle to estimate FDR cutoff
auto_shuffle <- function(input_matrix, shuffle_num = 1000, fdr = 0.1, 
                         break_threshold = 0.01, patience = 10) {
  qile <- numeric(shuffle_num)
  
  for (i in seq_len(shuffle_num)) {
    # Shuffle rows only
    shuffled <- input_matrix[sample(nrow(input_matrix)), , drop = FALSE]
    rand_cor <- cor(t(shuffled), method = "spearman")
    rand_dist <- mat2dist(rand_cor)
    qile[i] <- quantile(rand_dist, 1 - fdr)
    
    # Early stopping
    if (i > patience) {
      recent_range <- max(qile[(i - patience + 1):i]) - min(qile[(i - patience + 1):i])
      if (recent_range < break_threshold) {
        message("Auto-shuffle converged at iteration ", i)
        return(qile[i])
      }
    }
  }
  
  qile[shuffle_num]
}

# Recursive cluster extraction
extract_clusters <- function(mean_data, gene_list, threshold1 = 0.017, threshold2 = 0.024){
  temp_df <- mean_data[gene_list, , drop = FALSE]
  temp_cor <- cor(t(temp_df), method = "spearman")
  p1 <- cor_matrix_percentile(temp_cor, 0.1)
  p2 <- cor_matrix_percentile(temp_cor, 0.5)
  
  clusters <- list()
  removed <- character()
  
  if (p1 > threshold1 & p2 > threshold2){
    clusters <- list(gene_list)
  } else {
    new_cut <- cuttree_k2(temp_cor)
    for (k in 1:2){
      sub_genes <- names(new_cut[new_cut == k])
      if (length(sub_genes) >= 10){
        sub_result <- extract_clusters(mean_data, sub_genes, threshold1, threshold2)
        clusters <- c(clusters, sub_result$clusters)
        removed <- c(removed, sub_result$removed)
      } else {
        removed <- c(removed, sub_genes)
      }
    }
  }
  
  list(clusters = clusters, removed = removed)
}

# ------------------------------
# Main Workflow
# ------------------------------

# Load data
data <- read.csv(opt$expression_data, row.names = 1, check.names = FALSE)

# Mean-center genes (optional for Spearman)
mean_data <- data - rowMeans(data)

# Prepare genes of interest
if (!is.null(opt$genes)) {
  if (tolower(opt$genes) == "all") {
    genes <- rownames(mean_data)
  } else {
    genes <- str_split(opt$genes, ",")[[1]]
    genes <- genes[genes %in% rownames(mean_data)]
    if (length(genes) == 0) stop("No input genes match rows in the expression matrix.")
  }
} else {
  genes <- rownames(mean_data)  # default to all genes if not provided
}

# Calculate correlation
cor_mat <- cor(t(mean_data[genes, , drop = FALSE]), method = "spearman")

# Estimate FDR cutoff
fdr_cutoff <- auto_shuffle(cor_mat, fdr = opt$fdr, 
                           break_threshold = opt$stabilization_threshold, 
                           patience = opt$patience)

# Cut dendrogram using FDR cutoff
subtree <- cut_dend_tree(cor_mat, fdr_cutoff)

# Create output directory
dir.create(opt$out, showWarnings = FALSE)

# Extract clusters recursively and save CSVs
for (cluster_id in unique(subtree)){
  genes_in_cluster <- names(subtree[subtree == cluster_id])
  result <- extract_clusters(mean_data, genes_in_cluster)
  
  # Assign cluster IDs
  cluster_records <- unlist(lapply(seq_along(result$clusters), function(i){
    setNames(result$clusters[[i]], rep(i, length(result$clusters[[i]])))
  }))
  removed_records <- setNames(result$removed, rep("Removed", length(result$removed)))
  all_records <- c(cluster_records, removed_records)
  
  write.csv(data.frame(id = names(all_records), cluster = all_records),
            file.path(opt$out, paste0("tree_", cluster_id, ".csv")),
            row.names = FALSE)
}
