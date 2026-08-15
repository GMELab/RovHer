#!/bin/bash
# Master run script — RovHer rare-variant heritability (RARity) pipeline
#
# Sets all user parameters (directories, settings, and input files) in one
# place, then calls the 4 scripts sequentially to estimate the
# heritability captured by top-ranked variant subsets:

#   1_split_prop_blks.r   Split clumped RVs into top-proportion bins (top 1%..95%)
#   2_build_geno_blks.r   Extract/split genotype blocks for each subset
#   3_align_geno_pheno.r  Align + process genotype and phenotype matrices by IID
#   4_exome_wide_h2.r     Estimate per-block h2 and aggregate to a trait-level value
#
# Steps:
#       Step 1: Split clumped RVs into lists of top 1%, 5%, 10%, ..., 100% 
#       Step 2: Split geno blocks into smaller number of RVs (i.e. 5000 RVs)
#       Script 3. Create genotype matrix and align with phenotype matrix in /3_NORM_MAC2_GENO_0.1LD50kWIN_aligned_trait 
#       Script 4. Calculate adj_R2 (h2) in /4_H2_RESULTS

# Directories
DIR_WORK="/your/repo/dir"
cd $DIR_WORK

# Download github
git clone https://github.com/GMELab/RovHer.git
cd RovHer

DIR_SCRIPT="./Scripts"
DIR_WORK="./sample"
mkdir -p $DIR_WORK

# Variables
anno_name="RovHer" 
prop_blks=(1,5,10,15,20,25,30,35,40,50,60,70,75,80,85,90,95)
sort="descending" 
traits=("LDL_direct")
mode="regular" # regular is the R implementation; "lmutils" is the faster Rust-based implementation 

# Input files
input_file="${DIR_WORK}/sample/clumped_plinkids.txt" #clumped_variants_list.txt"
# Example: 
    # 1:92837557:A:G
    # 1:247424041:G:A

score_file="${DIR_WORK}/sample/master_score_file.txt"
# Example: 
    # PLINK_SNP_NAME  GENEBASS_AF     LOF_DISRUPTIVE  Gene.refGene    RovHer
    # 1:925948:A:C    6.2487e-05      0       SAMD11  0.89136924763911
    # 1:925950:G:C    1.3295e-06      0       SAMD11  0.892169307714771

pheno_file="${DIR_WORK}/sample/pheno_file.txt"
# Example:
    # IID	LDL_direct
    # p1	11.919
    # p2	6.948

PC_file="${DIR_WORK}/sample/PCs_1_20.txt"
# Example:
    # eid	AGE	SEX	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15	PC16	PC17	PC18	PC19	PC20
    # p1	62	1	-13.9042	6.6688	-0.982632	-1.46579	-0.0511454	2.64369	2.53286	1.23257	-1.51514	-1.2893	-3.17964	1.36194	-2.0918	-4.34219	-0.210294	0.37539	0.20819	0.406243	1.42201	-2.02739


########################################################################
# Step 1: Split clumped RVs into lists of top 1%, 5%, 10%, ..., 100% 
#   - output:  list of clumped RVs split by proportion of blocks
########################################################################

Rscript "${DIR_SCRIPT}/1_split_prop_blks.r" $input_file $anno_name $prop_blks $sort $score_file $DIR_WORK


########################################################################
# 2: Split geno blocks into smaller number of RVs (i.e. 5000 RVs)
########################################################################
GENOTYPE_DIR="${DIR_WORK}/sample/clumped_geno_matrices"
Rscript ${DIR_SCRIPT}/2_build_geno_blks.r $DIR_WORK $input_file $prop_blks $GENOTYPE_DIR

############################################################################
# SCRIPT 3 & 4: exome-wide per-block heritability calculation 
############################################################################

for top in 1; do
  for trait in "${traits[@]}"; do
    cores=2
    threads=2
    Rscript "${DIR_SCRIPT}/3_align_geno_pheno.r" ${anno_name} ${trait} ${DIR_WORK} ${top} ${cores} ${pheno_file} ${PC_file}
    export TMPDIR=/tmp 
    Rscript -e 'Sys.setenv(TMPDIR = "/tmp"); source(paste0("'$DIR_SCRIPT'", "/4_exome_wide_h2.r"))' ${anno_name} ${trait} ${threads} ${DIR_WORK} ${top} ${cores} ${mode}
  done
done

