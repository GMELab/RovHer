#!/usr/bin/Rscript
# DESCRIPTION
#   Aligns and processes clumped genotype and phenotype matrices in preparation
#   for heritability estimation. Supports two layouts:
#     - all chromosomes processed as ONE RARity block ("CHR_all" input), or
#     - each chromosome processed as separate block(s)  <-- default.
#
# INPUTS
#   - Unaligned/unprocessed clumped genotype and phenotype matrices
#   - Covariate/PC file, e.g.
#     "COVARS_UKB_MERGED_PHENO_..._FINAL_INCLUDE_PCs.RData"
#
# STEPS
#   1) Align GENOTYPE and PHENOTYPE matrices (by 'eid')
#   2) Continuous PHENOTYPE processing: quantile-normalize, adjust for
#      age/sex/20 PCs, then standardize (mean 0, sd 1)
#   3) GENOTYPE processing: NA-impute, drop MAC < 2 columns, normalize
#
# OUTPUTS
#   - Filtered/processed/aligned GENO in '/3_GENO_aligned_<trait>'
#   - Filtered/processed/aligned PHENO in '/3_PHENO_aligned_<trait>'
#
# COMMAND-LINE ARGUMENTS  (as actually parsed below; 9 total)
#   1. anno       (character) Annotation label, e.g. "yhat", "RovHer".
#   2. trait      (character) Trait name, e.g. "LDL_direct".
#   3. DIR_WORK   (character) Working directory root.
#   4. top        (character) Proportion bin, e.g. "100" or "80".
#   5. cores      (numeric)   Worker count for the parallel cluster.
#   6. chr_start  (numeric)   First chromosome (per-chromosome mode).
#   7. chr_end    (numeric)   Last chromosome (per-chromosome mode).
#   8. pheno_file (character) Phenotype table (eid + trait).
#   9. PC_file    (character) Age/sex/PC covariate table.
#   NOTE: the original usage line listed only 6 args and omitted top/pheno/PC;
#         the parsing below is the source of truth.
#
# EXAMPLE
#   Rscript 3_align_geno_pheno.R $anno $trait $DIR_WORK $top $cores \
#           $chr_start $chr_end $pheno_file $PC_file

############################################################################################
args <- commandArgs(trailingOnly = TRUE)

anno <- as.character(args[1])   # annotation label, e.g. "yhat"
trait <- as.character(args[2])   # trait name
DIR_WORK <- args[3]                 # working directory root
top  <- as.character(args[4])   # proportion bin, e.g. "100" or "80"
cores <- as.numeric(args[5])     # parallel worker count
pheno_file <- as.character(args[6])   # phenotype table (eid + trait)
PC_file  <- as.character(args[7])   # age/sex/PC covariate table

chr_start  <- 1
chr_end <- 22

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("MBESS"))
suppressMessages(library("corpcor"))    # provides pseudoinverse()
suppressMessages(library("doParallel"))
suppressMessages(library("parallel"))


registerDoParallel(makeCluster(cores))

total_start_time <- Sys.time()
cat("\n")
cat("=========== SCRIPT 3: Align geno and pheno (chr ", chr_start, "-", chr_end, ") ===========\n\n", sep = "")

# if pheno_file or PC_file do not exist, stop
if (!file.exists(pheno_file)) {
  stop("Error: pheno_file does not exist: ", pheno_file)
}
if (!file.exists(PC_file)) {
  stop("Error: PC_file does not exist: ", PC_file)
}

################################################################################
# OUTPUT DIRECTORIES
################################################################################

root1 <- paste0(DIR_WORK, "/h2_result")
if (!dir.exists(root1)) { dir.create(root1) }

root <- paste0(root1, "/top", top)
if (!dir.exists(root)) { dir.create(root) }
cat("3. Output root :", root, "\n\n")

# Output GENO and PHENO directories
GENO_OUTDIR <- paste0(root, "/3_GENO_aligned_", trait)
if (!dir.exists(GENO_OUTDIR)) { dir.create(GENO_OUTDIR) }

PHENO_OUTDIR <- paste0(root, "/3_PHENO_aligned_", trait)
if (!dir.exists(PHENO_OUTDIR)) { dir.create(PHENO_OUTDIR) }

# Output files
out_file1           <- paste0(PHENO_OUTDIR, "/PHENOS_", trait, "_top", top, "_IMPUTED_aligned.RData")  
norm_df_path        <- paste0(PHENO_OUTDIR, "/PHENOS_", trait, "_IMPUTED_aligned_QNORM_RESID_COVAR_normz.RData")
norm_df_path_no_eid <- paste0(PHENO_OUTDIR, "/PHENOS_", trait, "_IMPUTED_aligned_QNORM_RESID_COVAR_normz_no_eid.RData") 

################################################################################
# INPUT DIRECTORIES
################################################################################

# Unaligned genotypes from the clumping pipeline (SCRIPT 2 output)
GENO_DIR <- paste0(DIR_WORK, "/h2_result/top", top, "/2_GENO_RDATA")

cat("1. Genotypes :", GENO_DIR, "\n\n")
if (!dir.exists(GENO_DIR) || length(list.files(GENO_DIR)) == 0) {
  stop("Error: GENO_DIR does not exist or is empty")
}
# Remove a stray nohup.out if present
if (file.exists(paste0(GENO_DIR, "/nohup.out"))) {
  file.remove(paste0(GENO_DIR, "/nohup.out"))
}

##########################################################################################################################################################################################
# Mean-impute NAs (expects a data.frame column, not a data.table column)
NA2mean  <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
# Standardize to mean 0, sd 1  (z = (x - mean) / sd)
normFunc <- function(x) { (x - (mean(x, na.rm = TRUE))) / sd(x, na.rm = TRUE) }
# Rank-based inverse-normal (quantile) transform
quantNorm <- function(x) { qnorm(rank(x, ties.method = "average") / (length(x) + 1)) }

# ---------- Start of pheno processing  ----------

if (!file.exists(norm_df_path)) {  

  # ---------------------------------------------------------------------------
  # STEP 1) Align GENOTYPE and PHENOTYPE matrices
  # ---------------------------------------------------------------------------

  # Reference genotype file (first file in GENO_DIR)
  geno_files <- list.files(GENO_DIR, full.names = TRUE)
  if (length(geno_files) > 0) {
    ref_geno <- geno_files[1]
  } else {
    stop("No files found in GENO_DIR:", GENO_DIR, "\n")
  }
  base::load(ref_geno)   # loads object: GENE_DF
  cat("GENE_DF:", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], "\n")

  # Load phenotype matrix ('eid' + trait)
  PHENO <- fread(pheno_file)
  cat(trait, ":", dim(PHENO)[1], "x", dim(PHENO)[2], "\n")
  print(head(PHENO, 2))
  cat("\n")

  # Restrict PHENO to the eids present in GENE_DF, in GENE_DF order
  colnames(PHENO)[1]   <- "eid"
  colnames(GENE_DF)[1] <- "eid"
  index <- match(GENE_DF$eid, PHENO$eid)
  PHENO <- PHENO[index, ]
  PHENO <- PHENO[order(match(PHENO$eid, GENE_DF$eid)), ]
  PHENO <- PHENO[1:nrow(GENE_DF), ]

  # Drop participants with NA values
  na_count <- sum(is.na(PHENO[, 2]))
  cat(trait, "PHENO before NA-removal:", nrow(PHENO), " (Missing values:", na_count, ")\n")

  PHENO <- na.omit(PHENO) 
  cat(trait, "PHENO AFTER NA-removal:", nrow(PHENO), "\n\n")

  # Re-align GENE_DF eids to the surviving PHENO eid
  index   <- match(PHENO$eid, GENE_DF$eid)
  GENE_DF <- GENE_DF[index, ]
  GENE_DF <- GENE_DF[1:nrow(PHENO), ]

  # Confirm the eid columns match in order and length
  identical_ordered <- identical(GENE_DF$eid, PHENO$eid)
  identical_nrow    <- nrow(GENE_DF) == nrow(PHENO)

  if (identical_ordered && identical_nrow) {
    cat("'eid' in GENO and PHENO are aligned\n")
  } else {
    cat("GENO and PHENO are not aligned or the number of rows don't match.\n")
    if (!identical_ordered) {
      stop("'eid' columns are mismatched\n")
    }
    if (!identical_nrow) {
      cat("Rows in GENE_DF:", nrow(GENE_DF), "\nPHENO:", nrow(PHENO), "\n")
    }
  }
  cat("Aligned GENE_DF :", dim(GENE_DF)[1], "x", dim(GENE_DF)[2],
      " PHENO:", dim(PHENO)[1], "x", dim(PHENO)[2], "\n\n")
  base::save(PHENO, file = out_file1)

  # ---------------------------------------------------------------------------
  # STEP 2) Continuous phenotype processing
  # ---------------------------------------------------------------------------

  # (2a) Quantile-normalize the phenotype column
  cat("PHENO step 1: Quantile normalization....\n")
  PHENO <- as.matrix(PHENO)

  PHENO <- PHENO[, c(1, 2)] # keep eid + trait
  PHENO <- PHENO[complete.cases(PHENO), ] 

  PHENO[, 2] <- quantNorm(PHENO[, 2])

  colnames(PHENO)[2] <- trait
  colnames(PHENO)[2] <- paste(colnames(PHENO)[2], "Q_NORM", sep = "_")

  # (2b) Adjust for age, sex, and 20 PCs
  cat("PHENO step 2: Adj for age, sex, 20 PCs...\n")
  covar_traits <- c("eid", "AGE", "SEX",
                    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                    "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
  AgeSexPCA <- fread(PC_file)

  df <- merge(AgeSexPCA, PHENO[, 1:2], by = "eid")
  lm_f <- function(x) {
    x <- residuals(lm(
      data = df,
      formula = x ~ AGE + as.factor(SEX) +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
        PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20))
  }

  # Regress out covariates from the trait column
  trait_column <- names(df)[ncol(df)]
  resid        <- lm_f(df[[trait_column]])
  resid_final  <- data.frame(eid = df$eid)
  resid_final[[trait_column]] <- resid

  # (2c) Standardize residuals to mean 0, variance 1
  cat("PHENO step 3: Standardize mean and variance ...\n\n")
  norm_trait_values <- normFunc(resid_final[, 2])
  second_col_name   <- colnames(resid_final)[2]

  norm_df <- data.frame(eid = resid_final[, 1])
  norm_df[[second_col_name]] <- norm_trait_values
  cat("Final 'norm_df':", dim(norm_df)[1], "x", dim(norm_df)[2], "\n")
  print(head(norm_df, 2))
  cat("\n")

  # Save with and without the 'eid' column
  base::save(norm_df, file = norm_df_path)
  norm_df <- as.matrix(norm_df[, -1])            # 173599 x 1 (drop eid)
  base::save(norm_df, file = norm_df_path_no_eid)
  rm(df, resid, AgeSexPCA, covar_traits)
  invisible(gc()) 
} else {
  cat("Pheno files complete, skip to GENO processing...\n\n")
}


# End of all phenotype processing... skips to here
###################################################################################################
# STEP 3) GENOTYPE PROCESSING
#       (a)  NA-impute 
#       (b)  MINMAC: removes RV columns with a MAC < 2.
#       (c)  NORMALIZE to mean = 0 and sd = 1
#################################################################################
setwd(GENO_OUTDIR)

# Load the aligned, standardized phenotype (with eid) as a data.table
base::load(norm_df_path)
setDT(norm_df)

geno_files <- list.files(GENO_DIR, full.names = TRUE)

# "CHR_all" single-block files (exclude any "old" files)
chr_all_files <- geno_files[grepl("CHR_all", geno_files) & !grepl("old", geno_files)]

# MODE A: all chromosomes combined in a single "CHR_all" block
# ---------------------------------------------------------------------------
if (length(chr_all_files) > 0) {
  counter    <- 0
  start_time <- Sys.time()
  cat("======== SCRIPT 3:", anno, trait, "all chr in a single block  =======\n\n")

  for (GENO_path in chr_all_files) {
    outfile <- paste0(GENO_OUTDIR, "/NORM_MAC2_BRIT_NONBRIT_IMPUTE_CHR_all_", counter, ".RData")
    cat("Processing GENO_path:", GENO_path, "\n")
    cat("Output file:", outfile, "\n")

    if (!file.exists(GENO_path) || file.exists(outfile)) {
      next
    }

    base::load(GENO_path)   # loads object: GENE_DF
    # Defensive: some upstream blocks may have been saved as 'temp_GENO'
    if (exists("temp_GENO")) {
      GENE_DF <- temp_GENO
      rm(temp_GENO)
    }
    setDT(GENE_DF)
    colnames(GENE_DF)[1] <- "eid"

    # Keep only GENO rows whose eid is present in norm_df (inner join on eid)
    GENO_aligned <- merge(GENE_DF, norm_df[, .(eid)], by = "eid", allow.cartesian = TRUE)
    GENE_DF <- GENO_aligned[, names(GENE_DF), with = FALSE]
    if (names(GENE_DF)[ncol(GENE_DF)] == paste0(trait, "_Q_NORM")) {
      GENE_DF <- GENE_DF[, -ncol(GENE_DF), with = FALSE]
    }
    cat(paste("GENO: ", paste(dim(GENE_DF), collapse = " x "),
              "  PHENO: ", paste(dim(norm_df), collapse = " x "), "\n"))

    # STEP 1: NA-impute (columns 2..ncol; column 1 is eid)
    GENE_DF <- as.data.frame(GENE_DF)
    for (i in 2:ncol(GENE_DF)) {
      GENE_DF[, i] <- NA2mean(GENE_DF[, i])
    }

    # STEP 2: drop columns with MAC < 2
    id_col <- NULL
    if ("eid" %in% names(GENE_DF)) {
        id_col  <- GENE_DF[, "eid", drop = FALSE]
        GENE_DF <- GENE_DF[, setdiff(names(GENE_DF), "eid"), drop = FALSE]
    }
    non_num <- names(GENE_DF)[!sapply(GENE_DF, is.numeric)]
    if (length(non_num) > 0) {
        cat("WARNING: dropping non-numeric columns:", paste(non_num, collapse = ", "), "\n")
        GENE_DF <- GENE_DF[, sapply(GENE_DF, is.numeric), drop = FALSE]
    }
    cs      <- colSums(GENE_DF, na.rm = TRUE)
    keep    <- cs >= 2 & !is.nan(cs)
    GENE_DF <- as.data.frame(GENE_DF[, keep, drop = FALSE])

    if (!is.null(id_col)) {
        GENE_DF <- cbind(id_col, GENE_DF)
    }
    cat("MAC=2 filtered:", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], "\n")

    # STEP 3: normalize (columns 2..ncol; column 1 is eid)
    for (i in 2:ncol(GENE_DF)) {
      GENE_DF[, i] <- normFunc(GENE_DF[, i])
    }

    # Save with the 'eid' column (eid is only removed downstream in SCRIPT 4)
    base::save(GENE_DF, file = outfile)
    cat("Duration:", round(as.numeric(Sys.time() - start_time, units = "mins"), 2), "mins\n\n")
    rm(GENE_DF)
    invisible(gc())
    counter <- counter + 1
  }

# MODE B (default): each chromosome/block processed separately
# ---------------------------------------------------------------------------
} else {

  for (chr in chr_start:chr_end) {
    for (blk in 0:1000) {

      GENO_path <- paste0(GENO_DIR, "/CHR_", chr, ".clumped_", sprintf("%02d", blk), ".RData")
      outfile   <- sprintf(paste0(GENO_OUTDIR, "/NORM_MAC2_IMPUTE_CHR_", chr, "_", sprintf("%02d", blk), ".RData"))

      if (!file.exists(GENO_path) || file.exists(outfile)) {
        next
      }

      if (file.exists(GENO_path)) {
        cat("======== SCRIPT 3:", anno, trait, "Chr", chr, "blk", sprintf("%02d", blk), "========\n")
        start_time <- Sys.time()
        base::load(GENO_path)   # loads object: GENE_DF
        # Defensive: some upstream blocks may have been saved as 'temp_GENO'
        if (exists("temp_GENO")) {
          GENE_DF <- temp_GENO
          rm(temp_GENO)
        }
        setDT(GENE_DF)
        colnames(GENE_DF)[1] <- "eid"

        # Keep only GENO rows whose eid is present in norm_df (inner join)
        GENO_aligned <- merge(GENE_DF, norm_df[, .(eid)], by = "eid", allow.cartesian = TRUE)
        GENE_DF <- GENO_aligned[, names(GENE_DF), with = FALSE]
        if (names(GENE_DF)[ncol(GENE_DF)] == paste0(trait, "_Q_NORM")) {
          GENE_DF <- GENE_DF[, -ncol(GENE_DF), with = FALSE]
        }
        cat(paste("GENO: ", paste(dim(GENE_DF), collapse = " x "),
                  "  PHENO: ", paste(dim(norm_df), collapse = " x "), "\n"))

        # STEP 1: NA-impute (columns 2..ncol; column 1 is eid)
        GENE_DF <- as.data.frame(GENE_DF)
        for (i in 2:ncol(GENE_DF)) {
          GENE_DF[, i] <- NA2mean(GENE_DF[, i])
        }

        # STEP 2: drop columns with MAC < 2 (eid moved to rownames first)
        rownames(GENE_DF) <- GENE_DF$eid
        GENE_DF$eid <- NULL

        cs <- colSums(GENE_DF, na.rm = TRUE)   # every column is now numeric
        GENE_DF <- as.matrix(GENE_DF[, cs >= 2 & !is.nan(cs), drop = FALSE])

        GENE_DF <- as.data.frame(GENE_DF)
        cat("MAC=2 filtered:", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], "\n")

        # STEP 3: normalize (all columns are variants at this point)
        for (i in seq_len(ncol(GENE_DF))) {
          GENE_DF[, i] <- normFunc(GENE_DF[, i])
        }

        # Restore 'eid' (from rownames) as the first column
        GENE_DF <- cbind(eid = rownames(GENE_DF), GENE_DF)
        rownames(GENE_DF) <- NULL

        # Save with the 'eid' column
        base::save(GENE_DF, file = outfile)
        ncols <- ncol(GENE_DF)
        cat("Duration:", round(as.numeric(Sys.time() - start_time, units = "mins"), 2), "mins\n\n")

      }  # end: one block
    }  # end: all blocks for one chromosome
  }   # end: all chromosomes
}

end_time   <- Sys.time()
time_taken <- as.numeric(end_time - total_start_time, units = "mins")
cat("Time:", round(time_taken, 2), "mins\n")
cat("Start:", format(total_start_time, "%B %d %H:%M:%S"),
    "   End:", format(end_time, "%B %d %H:%M:%S"), "\n")

cat("============ SCRIPT 3: Aligned geno and pheno data ==============\n\n")

