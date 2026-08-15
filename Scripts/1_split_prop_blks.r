#!/usr/bin/Rscript
#
# SCRIPT 1:  Split a clumped variant list into top-proportion subsets
#
# DESCRIPTION
#   Extracts the top X% of variants from a clumped variant list, ranked by an
#   annotation score. Variants are aligned to their scores, sorted (ascending,
#   descending, or random), and written out as one ".clumped" file per requested
#   proportion bin. These per-bin files are the input to SCRIPT 2 (2_geno_blks.R).
#
#   Optionally, if a gene-aggregation directory ('/9_GENO_AGG_RDATA') exists, the
#   script first removes low-MAC variants that were aggregated into single SNPs
#   before ranking/splitting.
#
# INPUT REQUIREMENTS
#   clumped file : first column holds SNP IDs (any column name is fine); may
#                  contain additional columns. SNP IDs must be ':'-delimited.
#   score file   : must contain a column named <anno_name> to align to the SNP
#                  IDs. SNP IDs must be ':'-delimited.
#
# COMMAND-LINE ARGUMENTS
#   1. input_file (path)   Clumped variant list.
#   2. anno_name  (string) Annotation/score column name to rank by.
#   3. prop_blks  (string) Comma-separated proportion bins, e.g. "1,5,10,100".
#   4. sort       (string) "ascending", "descending", or "random".
#   5. score_file (path)   Table of variant scores (must contain <anno_name>).
#   6. DIR_WORK   (path)   Working directory root.
#
# EXAMPLE
#   Rscript 1_split_prop_blk.R "$input_file" "$anno_name" "1,5,10" "$sort" "$score_file" "$DIR_WORK"
# =============================================================================

suppressWarnings(suppressMessages(base::library(data.table)))
suppressWarnings(suppressMessages(base::library(doParallel)))
suppressWarnings(suppressMessages(base::library(parallel)))

# NOTE: a cluster is registered but no foreach/%dopar% is used in this script, and it is never stopped via stopCluster().
cores <- 2
registerDoParallel(makeCluster(cores))

# -----------------------------------------------------------------------------
# Command-line arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1])   # clumped variant list
anno_name  <- as.character(args[2])   # score column to rank by
prop_blks  <- as.character(args[3])   # e.g. "1,5,10,15,20,25,30,40,50,60,70,80,90,100"
sort       <- as.character(args[4])   # "ascending" / "descending" / "random"
score_file <- as.character(args[5])   # table of variant scores (used below)
DIR_WORK   <- as.character(args[6])   # working directory root

cat("\ninput_file (clump):", input_file, "\n")
cat("\nanno_name:", anno_name, "\n")
cat("\nprop_blks:", prop_blks, "\n")
cat("\nsort:", sort, "\n")
cat("\nScore file:", score_file, "\n")
cat("\nDIR_WORK:", DIR_WORK, "\n")

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------
if (!file.exists(score_file)) {
  stop(paste("\n\nscore_file", score_file, "does not exist.\n\n"))
}
if (!file.exists(input_file)) {
  stop(paste("\n\ninput_file", input_file, "does not exist.\n\n"))
}

# =============================================================================
# 1. Load and normalize the clumped variant list
# =============================================================================
clumped_RVs <- fread(input_file, header = FALSE)
setnames(clumped_RVs, "PLINK_SNP_NAME")

validate_format <- function(snp_names) {
  pattern <- "^\\d+:[0-9]+:[ACGT]+:[ACGT]+$"
  all(grepl(pattern, snp_names))
}

# Validate the PLINK_SNP_NAME column
if (!validate_format(clumped_RVs$PLINK_SNP_NAME)) {
  stop("ERROR: The PLINK_SNP_NAME column has values that do not match the 'chr:pos:ref:alt' format.")
} else {
  cat("All values in PLINK_SNP_NAME are in the correct 'chr:pos:ref:alt' format.\n")
}
# Example resulting IDs:
#   4:110575688:G:A

# Drop blank/NA SNP IDs and standardize the column name
clumped_RVs <- clumped_RVs[clumped_RVs$PLINK_SNP_NAME != "" & !is.na(clumped_RVs$PLINK_SNP_NAME), ]
original_RV_count <- nrow(clumped_RVs)
cat("\nNumber of clumped RVs:", original_RV_count, "\n")


# =============================================================================
# 2. Load the score file and normalize its SNP-ID column
# =============================================================================
score_df <- fread(score_file)
cat("Loading", nrow(score_df), "variants with aligned scores...\n")
head(score_df, 2)
cat("\n")

# First column must contain SNP IDs delimited by ':' or '_'
first_col <- score_df[[1]]
cat("Class of col1:", class(score_df[[1]]), "\n")
print(head(first_col, 3))

if (any(grepl("[:]", first_col))) {
  setnames(score_df, old = colnames(score_df)[1], new = "PLINK_SNP_NAME")
} else {
  stop(paste("\n", score_file, "must have first column containing SNP ids that are partitioned by ':'\n"))
}

if (!anno_name %in% colnames(score_df)) {
  head(score_df, 2)
  stop(paste("\nColumn", anno_name, "not found in score_df."))
}

# =============================================================================
# 3. Align annotation scores to the (retained) clumped RVs, then sort
# =============================================================================

new_df <- merge(clumped_RVs, score_df[, .(PLINK_SNP_NAME, get(anno_name))],
                by = "PLINK_SNP_NAME", all.x = TRUE)
setnames(new_df, old = "V2", new = anno_name)

# Every retained RV must have a score
if (any(is.na(new_df[[anno_name]]))) {
  stop(paste("\nTheres ", sum(is.na(new_df[[anno_name]])), "NA in", anno_name,
             "column.\nPlease check if all PLINK_SNP_NAME in clumped_RVs have corresponding entries in score_df.\n"))
}

# Sort by the annotation column in the requested order
if (sort == "ascending") {
  new_df <- new_df[order(new_df[[anno_name]]), ]
} else if (sort == "descending") {
  new_df <- new_df[order(-new_df[[anno_name]]), ]
} else if (sort == "random") {
  new_df <- new_df[sample(nrow(new_df)), ]
} else {
  stop("Error: sort must be either 'ascending', 'descending', or 'random'.")
}
cat("\nAligned a total of", nrow(new_df), "variants in", sort, "order...\n\n")

# Keep only PLINK_SNP_NAME and the annotation column
new_df <- new_df[, c("PLINK_SNP_NAME", anno_name), with = FALSE]
print(head(new_df, 2))
cat("\n")
print(tail(new_df, 2))
cat("\n")
total_rows <- nrow(new_df)

outfile <- paste0(DIR_WORK, "/", anno_name, "_01range_all_CHR_agg_aligned.txt")
fwrite(new_df, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# =============================================================================
# 4. Split the ranked list into top-proportion ".clumped" files
# =============================================================================
if (is.character(prop_blks)) {
  prop_blks <- as.numeric(unlist(strsplit(prop_blks, ",")))
}


for (prop_blk in prop_blks) {

  num_rows <- ceiling(total_rows * prop_blk / 100)
  top_rows <- new_df[1:num_rows, ]

  export_file <- gsub("\\.txt$", paste0("_top", prop_blk, ".clumped"), input_file)

  fwrite(top_rows, file = export_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  cat(paste("Top", prop_blk, "%  |  ", nrow(top_rows), "variants \n"))
}

cat(paste("\nResults in:\n", input_file, "\n\n", outfile, "\n\n"))

# Example row format:
#   PLINK_SNP_NAME     RovHer
#   16:71715915:GA:G   12.344
