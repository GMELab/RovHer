#!/usr/bin/env Rscript
# Last updated: June 2025
# https://github.com/Keonapang/RovHer

# For a given protein-coding gene, retrieve all rare variants with RovHer scores
################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
cat("\n")

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
setDTthreads(threads = 15)

args <- commandArgs(trailingOnly = TRUE)
DIR_OUT <- args[1] 
GENE <- args[2] # GENE="LDLR APOB BRCA1"

if (length(args) < 2) {
  stop("Usage: Rscript ./Scripts/get_scores_per_gene.r <GENE> <DIR_OUT>")
}

if (!dir.exists(DIR_OUT)) {
  dir.create(DIR_OUT, recursive = TRUE)
  cat(paste("Directory created:", DIR_OUT, "\n"))
}

# Parse single or multiple genes
gene_list <- unlist(strsplit(GENE, " "))  
cat("Selected gene(s):", gene_list, "\n")
cat("Output directory:", DIR_OUT, "\n\n")

# Load master score file
cat("Loading master score file...\n")
score_file <- fread("All_RovHer_Scores.txt.gz")

# Iterate over each gene
combined_gene_rows <- data.table()
skipped_genes <- c()
for (gene in gene_list) {
  cat("Processing gene:", gene, "...\n")
  
  # Check if the gene is valid
  if (!gene %in% unique(score_file$Gene)) {
    cat(paste("Warning:", gene, "not found in the score file. Skipping...\n"))
    skipped_genes <- c(skipped_genes, gene)  # Add to skipped genes list
    next
  }

  gene_rows <- score_file[Gene == gene]
  combined_gene_rows <- rbind(combined_gene_rows, gene_rows)
  cat(dim(gene_rows)[1], "variants found for", gene, "\n")
}

if (length(skipped_genes) > 0) {
  cat("The following genes were skipped because they were not found in the score file:\n")
  cat(paste(skipped_genes, collapse = ", "), "\n")
} else {
  cat("All gene variants were found.\n")
}
if (nrow(combined_gene_rows) == 0) {
  stop("No rare variants found for the provided gene(s).")
}

# Save 
if (length(gene_list) == 1) {
  outfile_gene <- paste0(DIR_OUT, "/output_", gene_list[1], ".txt")
} else {
  outfile_gene <- paste0(DIR_OUT, "/output_geneset.txt")
}
outfile_gene <- paste0(DIR_OUT, "/output_geneset.txt")
fwrite(combined_gene_rows, file = outfile_gene, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

