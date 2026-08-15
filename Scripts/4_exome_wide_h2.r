#!/usr/bin/Rscript

# SCRIPT DESCRIPTION
#   Estimates exome-wide RV heritability (h2) for a single trait from a set of
#   pre-computed genotype "blocks". For each block, an adjusted R^2 is estimated
#   by regressing the normalized phenotype on the block's genotype matrix; block
#   estimates and their variances are then aggregated into a single trait-level
#   h2 with a 95% confidence interval. The heritability estimation method is RARity
#   by Pathan et. al 2024 (https://www.nature.com/articles/s41467-024-45407-8).
#
#   Two computation modes are supported:
#     - "regular" : base-R implementation (use if the lmutils package is
#                   unavailable or not working).
#     - "lmutils" : uses lmutils::calculate_r2() for the per-block R^2 step.
#
# USAGE
#   Rscript 4_exome_wide_h2.r <anno> <trait> <threads> <DIR_WORK> <top> <cores> <mode>
#
# ARGUMENTS
#   1. anno     (chr) Annotation label, e.g. "RovHer", "yhat", "delta".
#   2. trait    (chr) Trait name, e.g. "LDL_direct".
#   3. threads  (num)   Worker threads (physical vCPUs/cores available).
#   4. DIR_WORK (chr) Working directory root.
#   5. top      (num)   Top-bin identifier (selects the "top<top>" folder).
#   6. cores    (num)   Core parallelism (how many files held in memory).
#   7. mode     (chr) "regular" or "lmutils".
#
# INPUT FILES
#   1. norm_pheno_file : mean-imputed, NORMALIZED phenotype (.RData)
#   2. geno_file       : one genotype "block" (NA-imputed, MAC > 2, normalized)
#
# OUTPUT
#   Appended to TOTAL_H2_<trait>_exome.txt in the "/4_<trait>_H2_RESULTS" folder.
#
####################################################################################
rm(list = ls())
gc()

# Input arguments
args     <- commandArgs(trailingOnly = TRUE)
anno     <- as.character(args[1])   # annotation label (e.g. yhat, delta, refvar)
trait    <- as.character(args[2])   # trait name
threads  <- as.numeric(args[3])     # worker threads (vCPUs/cores available)
DIR_WORK <- args[4]                 # working directory root
top      <- as.numeric(args[5])     # top-bin identifier
cores    <- as.numeric(args[6])     # core parallelism (files held in memory)
mode     <- as.character(args[7])   # "regular" or "lmutils"

# top <- 80
# DIR_WORK <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/RovHer/revision_round4"
# trait <- "LDL_direct"
# threads <- 1 
# anno <- "RovHer"
# cores <- 2
# mode <- "regular"

# Validate arguments and configure lmutils parallelism
if (is.na(threads) || threads <= 0) {
  stop("Error: 'threads' must be a valid positive number.")
}
if (is.na(cores) || cores <= 0) {
  stop("Error: 'cores' must be a valid positive number.")
}

# Read more at https://github.com/GMELab/lmutils.r
lmutils::set_core_parallelism(cores)     # how many files to hold in memory
lmutils::set_num_worker_threads(threads) # how many vCPUs/cores available
lmutils::disable_predicted()             # toggle: enable_predicted() / disable_predicted()

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("MBESS"))
suppressMessages(library("corpcor"))     # provides pseudoinverse()

cat("\n=========== Estimating exome-wide h2 RV for", trait, "===========\n")
start_time <- Sys.time()
cat("START:", format(start_time, "%B %d %H:%M:%S"), "\n\n")


################################################################################
# INPUT: locate and validate directories and phenotype files
################################################################################
root1 <- paste0(DIR_WORK, "/h2_result")
root  <- paste0(root1, "/top", top)
cat("1. Output root :", root, "\n\n")
if (!file.exists(root)) {
  cat("1. root:", root, "\n\n")
  stop("root missing\n")
}

# Genotype and phenotype directories
GENO_DIR  <- paste0(root, "/3_GENO_aligned_",  trait)
PHENO_DIR <- paste0(root, "/3_PHENO_aligned_", trait)
if (!file.exists(GENO_DIR) || length(list.files(GENO_DIR)) == 0) {
  cat("2. GENO_DIR:", GENO_DIR, "\n\n")
  stop("Geno dir missing\n")
}
if (!file.exists(PHENO_DIR) || length(list.files(PHENO_DIR)) == 0) {
  cat("3. PHENO_DIR:", PHENO_DIR, "\n\n")
  stop("Pheno dir missing\n")
}

# Phenotype input files
norm_df_path  <- paste0(PHENO_DIR, "/PHENOS_", trait, "_IMPUTED_aligned_QNORM_RESID_COVAR_normz.RData")
norm_df_path_no_eid <- paste0(PHENO_DIR, "/PHENOS_", trait, "_IMPUTED_aligned_QNORM_RESID_COVAR_normz_no_eid.RData")

################################################################################
# OUTPUT: create results directory and define output file paths
################################################################################
H2_RESULTS <- paste0(root, "/4_", trait, "_H2_RESULTS/")
if (!dir.exists(H2_RESULTS)) { dir.create(H2_RESULTS) }
setwd(H2_RESULTS)

outfile1 <- paste0(H2_RESULTS, trait, "_H2_exome_raw.txt")            # per-block raw R^2
outfile2 <- paste0(H2_RESULTS, trait, "_H2_exome_raw_processed.txt")  # per-block, CIs added
outfile3 <- paste0(H2_RESULTS, "TOTAL_H2_", trait, "_exome.txt")      # trait-level total h2
if (file.exists(outfile3)) {
  stop(paste("\n\nOutput already exists: ", outfile3, "\n\n"))
}

################################################################################
#  1.  Load Phenotype and configure the file
################################################################################

if (!file.exists(norm_df_path_no_eid)) {
  base::load(norm_df_path)               # loads object: norm_df
  norm_df <- as.matrix(norm_df[, -1])    # 173599 x 1 (drop the eid column)
  base::save(norm_df, file = norm_df_path_no_eid)
} else {
  base::load(norm_df_path_no_eid)        # loads object: norm_df
  cat("\nnorm_df :", norm_df_path_no_eid, "\n")
}


################################################################################
# 2. Count available genotype blocks
################################################################################
list_of_rdata <- list.files(GENO_DIR, pattern = "\\.RData$", full.names = TRUE)
rdata_count   <- length(list_of_rdata)
cat("Number of .RData files in GENO_DIR:", rdata_count, "\n\n")
cat("=========== [ SCRIPT 4 :", anno, trait, "h2 RV] ==========", "\n\n")

get_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}
start_time <- Sys.time()

################################################################################
# 3. Estimate per-block heritability (R^2) using RARity method 
#    Read more at Pathan et. al 2024 (https://www.nature.com/articles/s41467-024-45407-8)
################################################################################

# MODE: "regular"  (base-R version) - recommended for small-scale use
# get_R2(): estimate adjusted R^2 (and its variance/CI) for one block.
#   chr        block identifier (used for logging/output)
#   Datas      genotype matrix for the block (n individuals x m variants)
#   pheno_name trait name (written to output)
#   pheno_norm normalized phenotype vector/matrix (n x 1)
#   outfile    file to append results to

if (mode == "regular") {
  get_R2 = function(chr,Datas, pheno_name, pheno_norm, outfile){ 
    n = dim(Datas)[1]
    m = dim(Datas)[2]
    C_matrix = crossprod((Datas)) 
    C_all = t(Datas) %*% pheno_norm		
    possibleError <- tryCatch( 
      inv_matrix<-solve(C_matrix),	
      error=function(e) e
      )
    if (inherits(possibleError, "error")) {
      inv_matrix <- pseudoinverse(C_matrix)
      write.table(paste("ERROR : on = ", chr, "for trait", pheno_name),
                  paste0(out_dir, "/Warning_H2.txt"),
                  append = TRUE, quote = FALSE, row.names = FALSE,
                  col.names = FALSE, sep = "\t")
    }
    rm(C_matrix)
    # Regression coefficients, predictions, and R^2 / adjusted R^2
    Betas     <- inv_matrix %*% C_all
    Predicted <- Datas %*% Betas
    R2        <- var(Predicted)
    adj_R2    <- 1 - (1 - R2) * (n - 1) / (n - m - 1)

    # 95% CI on R^2
    out <- ci.R2(R2, df.1 = m, df.2 = n - m, conf.level = 0.95)[
      c("Lower.Conf.Limit.R2", "Upper.Conf.Limit.R2")]

    # Variance of R^2 and of adjusted R^2 (per Olkin & Finn approximation)
    rho2             <- R2
    block_VarR2      <- ((4 * rho2) * (1 - rho2)^2 * (n - m - 1)^2) / (((n^2) - 1) * (n + 3))
    block_Var_adj_R2 <- ((n - 1) / (n - m - 1))^2 * block_VarR2
    adj_R2_perVar    <- adj_R2 / m

    write.table(cbind(pheno_name, R2, adj_R2, n, m, chr,
                      out$Lower.Conf.Limit.R2, out$Upper.Conf.Limit.R2,
                      adj_R2_perVar, block_VarR2, block_Var_adj_R2),
                outfile, append = TRUE, quote = FALSE,
                row.names = FALSE, col.names = FALSE, sep = "\t")
    cat(paste("Processed block:", chr, "for trait:", pheno_name, "\n"))
    rm(C_all)
    rm(Datas)
    rm(inv_matrix)
    gc()
  }

  # initialized output file
  write.table(
    data.frame(trait = character(), r2 = numeric(), adj_r2 = numeric(),
               n = integer(), m = integer(), data = character(),
               block_LCL_R2 = numeric(), block_UCL_R2 = numeric(),
               adj_R2_perVar = numeric(), block_VarR2 = numeric(),
               block_Var_adj_R2 = numeric()),
    file = outfile1, append = FALSE, quote = FALSE,
    row.names = FALSE, col.names = TRUE, sep = "\t")
         
  # Estimate R^2 for each genotype block
  for (rdata_file in list_of_rdata) {
    tryCatch({
      base::load(rdata_file)    # loads object: GENE_DF
      Datas <- as.matrix(GENE_DF[, -1])
      chr   <- tools::file_path_sans_ext(basename(rdata_file))
      get_R2(chr, Datas, trait, norm_df, outfile1)
    }, error = function(e) {
      cat("SKIPPED", rdata_file, ":", conditionMessage(e), "\n")
    })
  }

  # Recompute CIs and variances for every block
  results_subset <- fread(outfile1)

  for (i in seq_len(nrow(results_subset))) {
    r2     <- as.numeric(results_subset$r2[i])
    adj_r2 <- as.numeric(results_subset$adj_r2[i])
    n      <- as.numeric(results_subset$n[i])
    m      <- as.numeric(results_subset$m[i])

    # Skip rows with NA or out-of-range values (adjust range checks as needed)
    if (is.na(r2) || is.na(adj_r2) || is.na(n) || is.na(m) || n <= m || m <= 0) {
      next
    }

    ci_output <- ci.R2(r2, df.1 = m, df.2 = n - m, conf.level = 0.95)
    results_subset$block_LCL_R2[i]  <- ci_output$Lower.Conf.Limit.R2
    results_subset$block_UCL_R2[i]  <- ci_output$Upper.Conf.Limit.R2
    results_subset$adj_R2_perVar[i] <- adj_r2 / m

    rho2                             <- r2
    block_VarR2                      <- ((4 * rho2) * (1 - rho2)^2 * (n - m - 1)^2) / ((n^2 - 1) * (n + 3))
    results_subset$block_VarR2[i]       <- block_VarR2
    results_subset$block_Var_adj_R2[i]  <- ((n - 1) / (n - m - 1))^2 * block_VarR2
  }

  results_subset$trait <- trait
  results_subset <- results_subset[, c("trait", setdiff(names(results_subset), "trait")), with = FALSE]

  write.table(results_subset, file = outfile2, row.names = FALSE, quote = FALSE, sep = "\t")
}

# ---------------------------------------------------------------------------
# MODE: "lmutils"  (uses lmutils::calculate_r2 for the per-block step)
# see: https://github.com/GMELab/lmutils.r?tab=readme-ov-file#lmutilscalculate_r2
# ---------------------------------------------------------------------------
if (mode == "lmutils") {

  if (rdata_count > 0) {
    cat("\nCALCULATING h2 for", rdata_count, " RData files... \n")
    results <- lmutils::calculate_r2(list_of_rdata, norm_df)
    cat("Results: ", paste(dim(results), collapse = " x "), "\n")

    results_subset <- results[, c("r2", "adj_r2", "n", "m", "data")]

    valid_data <- grepl("CHR_[0-9]+_[0-9]+", results_subset$data)

    if (all(valid_data)) {
      cat("All file paths are valid. Truncating and rearranging...\n")

      results_subset$data <- sapply(strsplit(as.character(results_subset$data), "/"), function(x) tail(x, 1))
      results_subset$data <- sub(".*_(CHR_[0-9]+_[0-9]+).*", "\\1", results_subset$data)

      # Sort into chronological chromosome block order
      results_subset <- results_subset[order(
        as.numeric(sub("CHR_([0-9]+)_.*", "\\1", results_subset$data)),
        as.numeric(sub("CHR_[0-9]+_([0-9]+)", "\\1", results_subset$data))
      ), ]
    } else {
      cat("Done processing a single rarity block.\n")
    }
  }
  write.table(results_subset, file = outfile1, append = FALSE, quote = FALSE,
              row.names = FALSE, col.names = TRUE, sep = "\t")

  # -------------------------------------------------------------------------
  # 4. Compute CIs and variances:
  #    block_LCL_R2, block_UCL_R2, adj_R2_perVar, block_VarR2, block_Var_adj_R2
  #
  # NOTE: results_subset here has only c(r2, adj_r2, n, m, data); the block_*
  #       columns below are created on first assignment. If any row is skipped
  #       via `next`, the created column can end up shorter than nrow() and the
  #       final `results_subset$col <- ...` reassignment may error with a
  #       "replacement has N rows" mismatch. Present in original logic.
  # -------------------------------------------------------------------------
  for (i in seq_len(nrow(results_subset))) {
    r2     <- as.numeric(results_subset$r2[i])
    adj_r2 <- as.numeric(results_subset$adj_r2[i])
    n      <- as.numeric(results_subset$n[i])
    m      <- as.numeric(results_subset$m[i])

    # Skip rows with NA or out-of-range values (adjust range checks as needed)
    if (is.na(r2) || is.na(adj_r2) || is.na(n) || is.na(m) || n <= m || m <= 0) {
      next
    }

    ci_output <- ci.R2(r2, df.1 = m, df.2 = n - m, conf.level = 0.95)
    results_subset$block_LCL_R2[i]  <- ci_output$Lower.Conf.Limit.R2
    results_subset$block_UCL_R2[i]  <- ci_output$Upper.Conf.Limit.R2
    results_subset$adj_R2_perVar[i] <- adj_r2 / m

    rho2                             <- r2
    block_VarR2                      <- ((4 * rho2) * (1 - rho2)^2 * (n - m - 1)^2) / ((n^2 - 1) * (n + 3))
    results_subset$block_VarR2[i]       <- block_VarR2
    results_subset$block_Var_adj_R2[i]  <- ((n - 1) / (n - m - 1))^2 * block_VarR2
  }

  # Move 'trait' to the first column and save
  results_subset$trait <- trait
  results_subset <- results_subset[, c("trait", setdiff(names(results_subset), "trait"))]

  write.table(results_subset, file = outfile2, row.names = FALSE, quote = FALSE, sep = "\t")
  cat(paste("\nOutput blocks: ", outfile2, "\n\n"))
}

#################################################################################
# 5. Aggregate per-block estimates into a single trait-level h2
#################################################################################

results_subset_tot <- fread(outfile2)

results_subset_tot[, adj_r2 := as.numeric(adj_r2)]
results_subset_tot[, duplicate := r2 == adj_r2]
results_subset_tot <- results_subset_tot[!duplicate | !duplicated(results_subset_tot[, .(r2, adj_r2)]), ]
results_subset_tot[, duplicate := NULL]   # remove temporary flag column

# Trait-level totals
N_RVs   <- base::sum(results_subset_tot$m, na.rm = TRUE)                  # total RVs
N       <- base::mean(results_subset_tot$n, na.rm = TRUE)                 # participants
ADJ_R2  <- base::sum(as.numeric(results_subset_tot$adj_r2), na.rm = TRUE) # total h2
STD2    <- base::sqrt(base::sum(results_subset_tot$block_Var_adj_R2))     # aggregated SD
LCL_adj <- ADJ_R2 - 1.96 * STD2                                          # 95% lower CI
UCL_adj <- ADJ_R2 + 1.96 * STD2                                          # 95% upper CI

# Write the trait-level result (create with header, else append)
output_headers <- c("TRAIT", "N", "N_RVs", "ADJ_R2", "STD2", "LCL_adj", "UCL_adj")
data_to_write  <- cbind(trait, N, N_RVs, ADJ_R2, STD2, LCL_adj, UCL_adj)

if (!file.exists(outfile3)) {
  cat(paste("Creating: ", outfile3, "\n"))
  suppressWarnings(write.table(data_to_write, file = outfile3, append = TRUE, quote = FALSE,
                               row.names = FALSE, col.names = output_headers, sep = "\t"))
} else {
  cat(paste("Appending: ", outfile3, "\n\n"))
  suppressWarnings(write.table(data_to_write, file = outfile3, append = TRUE, quote = FALSE,
                               row.names = FALSE, col.names = FALSE, sep = "\t"))
}

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n", trait, "h2: ", ADJ_R2, "\n")
cat("RVs (m):", N_RVs, "     Participants (n):", N, "\n")
cat("Duration:", round(duration), "mins\nEnd:", format(end_time, "%B %d %H:%M:%S"), "\n")
cat("========= SCRIPT 5:", anno, trait, " =======\n\n")  
cat("DIR_WORK:", DIR_WORK, "\n")

#################################################################################
# 6. CLEAN UP: remove intermediate .RData directories once results are written
#################################################################################

final_result <- fread(outfile3, header = TRUE, sep = "\t")
if (!is.null(final_result) && ncol(final_result) > 0 && nrow(final_result) > 0 && all(!is.na(final_result))) {
  unlink(GENO_DIR, recursive = TRUE)
  unlink(PHENO_DIR, recursive = TRUE)
  file.remove(outfile1)
} else {
  cat("Number of .RData:", rdata_count, "\n\n")
  cat("GENO_DIR:", GENO_DIR, "\n\n")
}