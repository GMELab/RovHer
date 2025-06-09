#!/usr/bin/env Rscript
# Last updated: June 2025
# https://github.com/Keonapang/RovHer

# Get Rovher scores for a list of input rare variants. 
# Inputs a .txt file with a list of variant PLINK IDs (formatted as "chr:pos:ref:alt")
# Outputs a .txt file with three tab-delimited columns: PLINK_SNP_NAME, Gene, RovHer_score

#########################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
install.packages('R.utils')

cat("\n")
library(data.table)
suppressMessages(library(R.utils))
setDTthreads(threads = 15)

args <- commandArgs(trailingOnly = TRUE)
INFILE <- args[1]
DIR_OUT <- args[2]

if (length(args) < 2) {
  stop("Usage: Rscript ./Scripts/get_scores.r <INFILE> <DIR_OUT>")
}

if (!dir.exists(DIR_OUT)) {
  stop(paste("Directory does not exist:", DIR_OUT))
}
if (!file.exists(INFILE)) {
  stop(paste("Input file does not exist:", INFILE))
} 

# Check if "All_RovHer_Scores.txt.gz" is in the same directory as this script

cat("Input:", INFILE, "\n\n")
cat("Output directory:", DIR_OUT, "\n\n")

# Load
cat("Loading input variant list...\n")
data <- tryCatch({
  fread(paste0(INFILE))
}, error = function(e) {
  stop(paste("Error reading input file:", e$message))
})

if (nrow(data) == 0) {
  stop(paste("Input file is empty:", INFILE))
}
cat(nrow(data), "variants found\n\n")

# add column header "PLINK_SNP_NAME"
colnames(data) <- "PLINK_SNP_NAME"

# Load master score file
cat("Loading master score file...\n")
SCORE_FILE <- "All_RovHer_Scores.txt.gz"
if (!file.exists(SCORE_FILE)) {
  stop(paste("RovHer scoring file (", SCORE_FILE, ") does not exist. Instructions: download from Zenodo and place it in the /RovHer directory."))
}
scores <- fread(SCORE_FILE)

# Check if data is formatted correctly
if (ncol(data) != 1) {
  stop("The input data does not consist of exactly one column.")
}
is_valid_variant <- grepl("^\\d+:\\d+:[A-Za-z]+:[A-Za-z]+$", data[[1]])
if (!all(is_valid_variant)) {
  stop("Not all values in the column are valid variant IDs.")
}

# Merge
cat("Retrieving RovHer scores...\n")
setkey(scores, PLINK_SNP_NAME)
setkey(data, PLINK_SNP_NAME)
data <- merge(data, scores, by = "PLINK_SNP_NAME", all.x = TRUE)

# Count number of missing values or NA in the "RovHer_Score" column
missing_scores <- sum(is.na(data$RovHer_Score))
if (missing_scores > 0) {
  cat("Done!", missing_scores, "variants do not have RovHer scores.\n\n")
} else {
  cat("Done! All variants have RovHer scores.\n\n")
}

# Save
outfile <- paste0(DIR_OUT, "/output_", sub("\\.txt$", "", basename(INFILE)), ".txt")
fwrite(data, file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Results saved to:", outfile, "\n")
