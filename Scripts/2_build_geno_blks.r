#!/usr/bin/Rscript
#
# OBJECTIVE
#   For each requested proportion bin, scan the full set of clumped genotype
#   matrices and keep only the variants present in that bin's variant list. Retained
#   columns are accumulated and written out as a smaller set of genotype blocks
#   (~5,000 variants each). Small subsets (< 1,000 variants) are written as a single
#   combined block spanning all chromosomes.
# 
# INPUT
#   - PLINK_list  : clumped variant IDs (one ".clumped" file per proportion bin,
#                   created by 1_split_prop_blk.r)
#   - GENOTYPE_DIR: all clumped ~250 exome-wide (chr1-22) genotype blocks
#                   (e.g. from '/8_GENO_0.1LD50kWIN_RDATA')
#
# OUTPUT
#   - A smaller/filtered set of ~5k-variant genotype blocks (RVs only), written
#     per proportion bin to '<DIR>/h2_result/top<top>/2_GENO_RDATA'.
#
# COMMAND-LINE ARGUMENTS
#   1. DIR          (character) Working directory root.
#   2. PLINK_list   (character) Base clumped-variant list (".txt"); the per-bin
#                               file "<base>_top<top>.clumped" is loaded per bin.
#   3. prop_blks    (character) Comma-separated proportion bins,
#                               e.g. "1,5,10,15,20,25,30,35,40,50,60,70,80,90".
#   4. GENOTYPE_DIR (character) Directory of clumped genotype matrices.
#   NOTE: the original usage line listed a 4th arg "cores", but the code reads
#         args[4] as GENOTYPE_DIR and never reads a 'cores' argument; the
#         parallel cluster below is hard-coded to 4 workers.
#
# EXAMPLE
#   DIR="/your/repo/dir"
#   PLINK_list="${DIR}/4_CLUMP_RESULT/yhat_01range_CHR_17.clumped"
#   prop_blks="1,5,10"
#   Rscript 2_geno_blks.R $DIR $PLINK_list $prop_blks $GENOTYPE_DIR

############################################################################################################
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
PLINK_list <- args[2] 
prop_blks <- as.character(args[3]) # c(1,5,10,15,20,25,30,35,40,50,60,70,80,90)
GENOTYPE_DIR <- as.character(args[4]) 


suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(doParallel)
  library(parallel)
}))
registerDoParallel(makeCluster(4))

cat("\n================== [ SCRIPT 2 ] ==================\n")

if (is.character(prop_blks)) {
  prop_blks <- as.numeric(unlist(strsplit(prop_blks, ",")))
}
cat("\nprop_blks:", prop_blks, "\n\n")

PLINK_list_base <- sub("\\.txt$", "", PLINK_list)
cat("PLINK_list_base:", PLINK_list_base, "\n\n")

######################################################################
# # Directories
######################################################################
# Input directory (clumped genotype matrices)
if (!dir.exists(GENOTYPE_DIR)) {
  stop("ERROR: No genotype matrix directory found")
}

# Root output directory
root <- paste0(DIR, "/h2_result")
if (!dir.exists(root)) { dir.create(root) }

# List all genotype matrices (includes files with "MAXMAC" in the filename)
geno_files <- list.files(path = GENOTYPE_DIR, full.names = TRUE)
geno_files <- sort(geno_files)
cat("Number of genotype files:", length(geno_files), "\n\n")

######################################################################
# MAIN LOOP: iterate over each proportion bin
######################################################################

start_time <- Sys.time()

for (top in prop_blks) {

  # ---------------------------------------------------------------------------
  # Set up per-bin output director
  # ---------------------------------------------------------------------------
  root2 <- paste0(root, "/top", top)
  if (!dir.exists(root2)) { dir.create(root2) }

  GENO_DIR <- paste0(root2, "/2_GENO_RDATA")
  if (!dir.exists(GENO_DIR)) { dir.create(GENO_DIR, recursive = TRUE)}
  if (length(list.files(GENO_DIR)) > 0) {
    file.remove(list.files(GENO_DIR, full.names = TRUE))
  }

  # ---------------------------------------------------------------------------
  # Load the per-bin clumped variant list (from 1_split_prop_blk.r)
  # ---------------------------------------------------------------------------
  PLINK_list <- paste0(PLINK_list_base, "_top", top, ".clumped")
  if (!file.exists(PLINK_list)) {
    stop(PLINK_list, "\n\n")
  }
  cat("\nLoading:", PLINK_list, "\n\n")
  var_list        <- fread(PLINK_list)   # input file
  num_top_variant <- nrow(var_list)

  # Vector of SNP names to keep (plus the "IID" column)
  PLINK_IDs        <- var_list$PLINK_SNP_NAME
  required_RV_cols <- c("IID", PLINK_IDs)

  # ---------------------------------------------------------------------------
  # Initialize
  # ---------------------------------------------------------------------------
  temp_GENO <- NULL  # data.table of retained columns
  col_count <- 0   
  new_blks  <- 0  

  # Chromosomes represented in this RV subset
  chromosomes <- sapply(strsplit(as.character(var_list[[1]]), ":"), "[", 1)
  unique_chr  <- as.numeric(sort(unique(chromosomes)))
  cat(paste0("\nFinding ", num_top_variant, " RVS in top ", top, "%\n"))

  # ===========================================================================
  # Search each genotype matrix for the requested variants
  # ===========================================================================
  cat("Searching for variants in chromosomes:", paste(unique_chr, collapse = ", "), "\n\n")

  SINGLE_BLK <- num_top_variant < 1000

  for (chr in unique_chr) {
    count <- 0  
    # Genotype files relevant to this chromosome
    target_files <- grep(paste0("_", chr, "\\."), geno_files, value = TRUE)
    if (length(target_files) == 0) {
      cat("No genotype matrices for chromosome", chr, "\n")
      next
    } else {
      cat("Found", length(target_files), "clumped GENO files for chr", chr, "\n")
    }

    for (geno_file in target_files) {
      original_basename <- basename(geno_file)
      new_basename <- sub(
        "clumped_\\d{2}",
        sprintf("clumped_%02d", as.numeric(count)),
        original_basename
      )
      out_path <- file.path(GENO_DIR, new_basename)
      cat("Output:", new_basename, "\n\n")

      cat("[SCRIPT 2: Top", top, "%  ", basename(geno_file), "]\n")
      base::load(geno_file)   # loads object: GENE_DF

      # Normalize PLINK column names (strip a trailing "_<ALT>" allele suffix)
      colnames(GENE_DF) <- ifelse(
        colnames(GENE_DF) == "IID",
        "IID",
        gsub("_[A-Z]$", "", colnames(GENE_DF))
      )

      # If non-IID columns use "_" delimiters (chr_pos_ref_alt), convert to ":"
      variant_cols <- setdiff(colnames(GENE_DF), "IID")
      if (length(variant_cols) > 0 && all(grepl("_", variant_cols)) && !any(grepl(":", variant_cols))) {
        cat("Renaming SNPs from chr_pos_ref_alt --> chr:pos:ref:alt\n")
        cols    <- colnames(GENE_DF)
        non_iid <- cols != "IID"
        # Replace only the first 3 underscores (chr_pos_ref_alt -> chr:pos:ref:alt)
        cols[non_iid] <- sub("^([^_]+)_([^_]+)_([^_]+)_(.+)$", "\\1:\\2:\\3:\\4", cols[non_iid])
        colnames(GENE_DF) <- cols
      }

      if (!all(grepl(":", PLINK_IDs)) || !all(grepl(":", colnames(GENE_DF)[colnames(GENE_DF) != "IID"]))) {
        stop("ERROR: PLINK_IDs and GENE_DF columns do not have the same format (script needs chr:pos:ref:alt)\n")
      }

      # -----------------------------------------------------------------------
      # Keep only the PLINK_list variants that belong to this chromosome
      # -----------------------------------------------------------------------
      pattern              <- paste0("^", chr, ":")
      required_RV_cols_chr <- required_RV_cols[grepl(pattern, required_RV_cols)]
      required_RV_cols_chr <- c("IID", required_RV_cols_chr)          # keep IID
      common_columns       <- intersect(names(GENE_DF), required_RV_cols_chr)

      if (length(common_columns) == 1) {   # only the IID column matched
        rm(GENE_DF)
        next
      }
      cat("Found", length(common_columns) - 1, "RVs\n")

      GENE_DF <- GENE_DF[, ..common_columns, with = FALSE]

      # Rename IID -> eid
      if ("IID" %in% names(GENE_DF)) {
        GENE_DF <- rename(GENE_DF, eid = IID)
      }

      # Append the filtered GENE_DF to the running temp_GENO
      if (is.null(temp_GENO)) {
        # First block: keep all columns
        temp_GENO <- GENE_DF
      } else {
        # Subsequent blocks: bind all columns except the shared "eid"
        stopifnot(identical(temp_GENO$eid, GENE_DF$eid))
        temp_GENO <- bind_cols(temp_GENO, GENE_DF[, -1])
        # temp_GENO <- bind_cols(temp_GENO, GENE_DF[, !colnames(GENE_DF) %in% "eid"])
      }
      num_participants <- nrow(temp_GENO)   # always valid here
      rm(GENE_DF); invisible(gc())

      # Running variant count (excludes the eid column)
      col_count <- col_count + ncol(temp_GENO) - 1
      cat("Variants found so far:", col_count, "\n\n")
      # print(head(temp_GENO[, 1:5]))

      # -----------------------------------------------------------------------
      # Flush the block to disk once it reaches ~5k columns (multi-block mode)
      # -----------------------------------------------------------------------
      if (!SINGLE_BLK && ncol(temp_GENO) >= 4999) {
        GENE_DF <- temp_GENO
        # FIXME: saves the object 'temp_GENO', but every other save (and SCRIPT 4)
        #        expects the object to be named 'GENE_DF'. Loading this block
        #        downstream will not find GENE_DF. Likely intended: save(GENE_DF, ...).
        save(temp_GENO, file = out_path)
        cat("Blk has reached 5k size:", dim(temp_GENO)[1], "x", dim(temp_GENO)[2], "\n")
        cat("Saved:", new_basename, "\n\n")

        temp_GENO <- NULL
        new_blks  <- new_blks + 1
        count     <- count + 1
      }
    }  # ---- end loop over genotype files for this chromosome ----

    # Save any remaining accumulated data (multi-block mode, per-chromosome)
    if (!SINGLE_BLK && !is.null(temp_GENO) && ncol(temp_GENO) > 0) {
      cat("Final temp_GENO:", dim(temp_GENO)[1], "x", dim(temp_GENO)[2], "\n")
      GENE_DF <- temp_GENO
      cat("\n")
      base::save(GENE_DF, file = out_path)
      cat("Saved :", out_path, "\n\n")

      temp_GENO <- NULL
      GENE_DF   <- NULL
      new_blks  <- new_blks + 1
    } else if (SINGLE_BLK && !is.null(temp_GENO)) {
      cat("Accumulated after chr", chr, ":", dim(temp_GENO)[1], "x", dim(temp_GENO)[2], "\n\n")
    }
  }  # ============ end loop over chromosomes ============

  # NOTE: in multi-block mode temp_GENO is NULL here, so num_participants
  #       becomes NULL; it is not used meaningfully afterward.
  num_participants <- nrow(temp_GENO)
  cat("Total RVS found:", col_count, "\n")

  # ===========================================================================
  # SINGLE_BLK: write one combined block spanning every chromosome
  # ===========================================================================
  if (SINGLE_BLK) {
    if (is.null(temp_GENO) || ncol(temp_GENO) <= 1) {
      cat("WARNING: No variants collected across any chromosome. Nothing to save.\n")
    } else {
      SINGLE_GENE_DF <- temp_GENO

      # Sanity check: did we recover every requested variant?
      found_vars <- setdiff(colnames(SINGLE_GENE_DF), "eid")
      cat("Chromosomes combined:", paste(unique_chr, collapse = ", "), "\n")
      cat("Dimensions:", nrow(SINGLE_GENE_DF), "x", ncol(SINGLE_GENE_DF), "\n")
      cat("Variants found:", length(found_vars), "/", num_top_variant, "\n")

      missing_vars <- setdiff(PLINK_IDs, found_vars)
      if (length(missing_vars) > 0) {
        cat("Variants NOT found (", length(missing_vars), "):\n", sep = "")
        print(missing_vars)
      }

      # Per-chromosome breakdown of the retained columns
      found_chr <- sapply(strsplit(found_vars, ":"), "[", 1)
      cat("\nVariants per chromosome:\n")
      print(table(found_chr))
      cat("\n================ SINGLE_GENE_DF ================\n")

      # Keep the GENE_DF alias so downstream code expecting `GENE_DF` still works
      GENE_DF <- SINGLE_GENE_DF

      new_basename <- sub(
        "CHR_\\d+\\.clumped_\\d+",
        "CHR_all",
        original_basename
      )
      single_out_path <- file.path(GENO_DIR, new_basename)
      base::save(GENE_DF, file = single_out_path)
      cat("Saved", nrow(GENE_DF), " participants x", ncol(GENE_DF) - 1, "variants\n\n")
      cat("Saved SINGLE_GENE_DF:", single_out_path, "\n\n")

      temp_GENO <- NULL
      new_blks  <- 1
    }
  }

}  # ============ end loop over proportion bins ============

end_time   <- Sys.time()
time_taken <- as.numeric(end_time - start_time, units = "mins")
cat("Time:", round(time_taken, 2), "mins\n\n")

cat("Results:", GENO_DIR, "\n\n")
cat("============ SCRIPT 2:Top", top, "% RVs============", "\n")