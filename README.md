# RovHer: heritability-optimized scores for the functional prioritization of rare missense variants
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15596103.svg)](https://doi.org/10.5281/zenodo.15596103)

<!-- TABLE OF CONTENTS -->
<a name="readme-top"></a>
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about">About</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#hardware">Hardware</a></li>
        <li><a href="#software-dependencies">Software Dependencies</a></li>
      </ul>
    </li>
    <li><a href="#usage-retrieve-pre-computed-scores">Usage: Retrieve Pre-computed Scores</a></li>
    <li><a href="#resources">Resources</a></li>
    <li><a href="#acknowledgements">Acknowledgments</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#citing-rovher">Citing RovHer</a></li>
    <li><a href="#references">References</a></li>
  </ol>
</details>

<!-- ABOUT -->
## About

Predicting **rare variants** (RVs; MAF < 1%) that influence complex disease risk is a major challenge in human genetics. We introduce **RovHer** (RV heritability-optimized scores), an unbiased, scalable method that scores missense RVs based on their probability of functional effect. RovHer employs the [Multivariate Adaptive Regression Splines](https://CRAN.R-project.org/package=earth) [1] model to integrate feature annotations with [Genebass](https://app.genebass.org/) [2] (N=394,841) exome-wide association study (ExWAS) summary statistics of height, which serves as the training trait. For the dependent variable, we used the **false discovery rate** across 4,927,152 RVs as a surrogate measure for the likelihood of variant functionality.

The raw outputs of RovHer therefore reflects the predicted FDR of each variant. We subsequently inverted them to obtain a "true positive" scale, where higher RovHer scores (closer to 1) represent a greater probability of a variant being functional, and lower scores (closer to 0) represent likely neutral variants.


**Performance**:
* prioritizes missense RVs that explain the greatest proportion of genome-wide trait variance, which are variants that more likely functional [3]

**Usage**:
* scores are trait-agnostic (i.e. not tied to a specific disease or phenotype)
* RovHer scores were generated for all possible rare SNVs deposited in the dbNSFP database; however, they are optimized for missense variants.

![Workflow Overview](RovHer%20workflow.png)
*An overview of the development and application of RovHer.*

# Getting started
### Hardware
RovHer can generate predictions on major operating systems, including GNU/Linux, macOS, and Windows. For biobank-scale analyses, we recommend Unix-based hardware with a minimum of 100GB RAM.

### Software dependencies

The package development version is tested on Linux operating systems. However, the CRAN packages required should be compatible with Windows, Mac, and Linux operating systems. Type the following into your R console:
   ```R
  install.packages("data.table")
  install.packages("R.utils")
  install.packages("dplyr")
  install.packages("tidyverse")
  install.packages("corpcor") # pseudoinverse fallback
  install.packages("MBESS") # confidence intervals on R²

  ```
To install `lmutils ` (optional), which is only required for `mode="lmutils"` in Step 4. Type the following into your R console:

    ```R
    install.packages(
      "https://github.com/mrvillage/lmutils.r/archive/refs/heads/master.tar.gz",
      repos = NULL) 
    ```
    
---

<!-- Usage: Retrieve pre-computed scores -->
# Usage: Retrieve pre-computed scores 

1. Clone the repo into your working directory. This step should take only a few seconds on any computer.
   ```sh
   mydir="/my/working/dir" # modify 
   cd $mydir

   git clone https://github.com/Keonapang/RovHer.git
   cd RovHer
   ```
2. Download `All_RovHer_Scores.txt.gz` from [Zenodo](https://zenodo.org/records/15596103?preview=1) and save it to your local repo directory as `/RovHer/All_RovHer_Scores.txt.gz`. This step should take up to several minutes.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15596103.svg)](https://doi.org/10.5281/zenodo.15596103)

3. Run automated scripts to obtain a subset of RovHer scores for downstream analyses in two ways:

| Option | Description | Output | Script name |
|--:|-----------|-----------|-----------|
|  A| Retrieve scores for a list of RVs in `variants.txt` | A `output_variants.txt` with tab-delimited columns: PLINK_SNP_NAME, Gene, RovHer_score | `get_scores.r` |
|  B| For a given protein-coding `gene` or set of genes, retrieve a list of scored RVs | A `output_gene.txt` or `output_geneset.txt` with tab-delimited columns: PLINK_SNP_NAME, Gene, RovHer_score | `get_scores_per_gene.r` |


For **option A**, the required input is a single text file `variants.txt` consisting of a column of PLINK IDs (no headers). Variants do not have to be sorted. For example `$INFILE` may look like:
|             |
|-------------|
|  1:10030:A:T| 
|  8:203440:G:C| 

   ```sh
    INFILE="$mydir/RovHer/Demo/variants.txt" # input list 
    DIR_OUT="$mydir/RovHer/Demo" # output directory 

    cd $mydir/RovHer
    Rscript Scripts/get_scores.r $INFILE $DIR_OUT
  ```

For **option B**, no input file is required. Simply specify a list of gene(s) in all capital letters:
  ```sh
    DIR_OUT="$mydir/RovHer/Demo" # output directory 

    GENE="LDLR BRCA1 APOB" # example of a gene set  
    # or
    GENE="LDLR" # example of one gene 

    cd $mydir/RovHer
    Rscript Scripts/get_scores_per_gene.r "$DIR_OUT" "$GENE"
  ``` 

<p align="right">(<a href="#readme-top">back to top</a>)</p>

The runtime for these scripts ranges from a few seconds up to a few minutes (depending on input variant file size) on a desktop computer (RAM: 16+ GB, CPU: 4+ cores, 3.3+ GHz/core).

---

# Benchmark: Rare variant heritability estimation with RARity 

Given a ranked list of variants (e.g. prioritized by the RovHer annotation score), this benchmarking pipeline estimates how much trait heritability is captured by a chosen subset of variants. This enables a comparison of **heritability-enrichment** across different variant prioritization methods. This is a four-step pipeline for estimating **exome-wide rare-variant (RV) heritability** from genotype and phenotype data, using the RARity approach. These scripts run sequentially:

| Script | Purpose | Outputs |
|--------|---------|---------|
|`1_split_prop_blks.R` | Splits a large variant list into top proportion bins | |
|`2_geno_blks.R` | Extract genotype blocks for a target variant subset | Filtered genotype blocks (`GENE_DF` objects) written to `<DIR>/h2_result/top<top>/2_GENO_RDATA/`, one set of blocks per proportion bin| 
|`3_align_geno_pheno.R` | Align both genotype and phenotype matrices by participant `eid` and processes them (mean-impute, quantile normalize, standardize...) | genotype blocks in `3_GENO_aligned_<trait>/` and phenotype (`norm_df`) in `3_PHENO_aligned_<trait>/|
|`4_exome_h2_RV.R` | Estimate per-block heritability then aggregate into one trait-level estimate| `TOTAL_H2_<trait>_exome.txt` in `4_<trait>_H2_RESULTS/`|


 **RARity** estimates heritability from *blocks* of approximately LD-independent variants. Within each block, the normalized phenotype is regressed on the block's genotype matrix, and the adjusted R² of that regression is an estimate of the heritability explained by the block. Block estimates (and their variances) are then summed into a single trait-level heritability with a 95% confidence interval.
 
By running the pipeline on different **variant subsets** — for example, the top *X%* of variants ranked by a functional annotation — you can measure how heritability concentrates among prioritized variants. This makes the pipeline useful both as a standalone heritability estimator and as an **evaluation metric for variant-prioritization methods** (does ranking variants by a given score concentrate heritability at the top?). 

## Upstream requirements

The pipeline assumes the following have already been produced (not included in this repo):

1. **Genotype matrices** per chromosome, saved as `.RData` files each containing a
   single object named `GENE_DF` — a table with an `IID`/`eid` column plus one column
   per variant. Variant IDs are in `chr:pos:ref:alt` (or `chr_pos_ref_alt`) format.
2. **A clumped variant list** (LD-pruned index variants), typically from PLINK
   `--clump`, split into **per-proportion files** by the upstream `1_split_prop_blk.r`
   step. Step 2 expects files named `<base>_top<top>.clumped`.
3. **A phenotype table** (`pheno_file`): first column `eid`, second column the trait.
4. **A covariate/PC table** (`PC_file`): `eid`, `AGE`, `SEX`, `PC1`…`PC20`.

---

## Directory structure

All outputs are written under `<DIR_WORK>/h2_result/top<top>/`:

```
<DIR_WORK>/h2_result/
└── top<top>/
    ├── 2_GENO_RDATA/                    # Step 2: filtered genotype blocks
    ├── 3_GENO_aligned_<trait>/          # Step 3: processed, aligned genotype blocks
    ├── 3_PHENO_aligned_<trait>/         # Step 3: processed phenotype (norm_df)
    └── 4_<trait>_H2_RESULTS/            # Step 4: heritability results
        ├── <trait>_H2_exome_raw.txt            # per-block raw R² (intermediate)
        ├── <trait>_H2_exome_raw_processed.txt  # per-block, with CIs/variances
        └── TOTAL_H2_<trait>_exome.txt          # ← final trait-level heritability
```

Here `<top>` is the proportion bin (e.g. `1` for top 1% of RVs, or `100` for all RVs) and
`<trait>` is the trait name (e.g. `LDL_direct`).


```bash

# Directories
DIR_SCRIPT="/your/repo/dir"
DIR_WORK="${DIR_SCRIPT}"
mkdir -p $DIR_WORK

# Variables to modify 
anno_name="RovHer" # score column name
prop_blks=(1,5,10,15,20,25,30,35,40,50,60,70,75,80,85,90,95)
traits=("LDL_direct")
mode="regular" # regular is the R implementation; "lmutils" is the Rust-based implementation
sort="descending" # "descending" (high scores = more functional) or "ascending" (lower scores = more functional) 

# Input files/directories
input_file="${DIR_WORK}/sample/clumped_plinkids.txt" #clumped_variants_list.txt"
score_file="${DIR_WORK}/sample/master_score_file.txt"
GENOTYPE_DIR="${DIR_WORK}/sample/clumped_geno_matrices"
pheno_file="${DIR_WORK}/sample/pheno_file.txt"
PC_file="${DIR_WORK}/sample/PCs_1_20.txt"

# Run
Rscript "${DIR_SCRIPT}/1_split_prop_blks.r" $input_file $anno_name $prop_blks $sort $score_file $DIR_WORK
Rscript ${DIR_SCRIPT}/2_build_geno_blks.r $DIR_WORK $input_file $prop_blks $GENOTYPE_DIR

for top in 1 5 10 15 20 25 30 35 40 50 60 70 75 80 85 90 95; do
  for trait in "${traits[@]}"; do
    cores=2
    threads=2
    Rscript "${DIR_SCRIPT}/3_align_geno_pheno.r" ${anno_name} ${trait} ${DIR_WORK} ${top} ${cores} ${pheno_file} ${PC_file}
    export TMPDIR=/tmp 
    Rscript -e 'Sys.setenv(TMPDIR = "/tmp"); source(paste0("'$DIR_SCRIPT'", "/4_exome_wide_h2.r"))' ${anno_name} ${trait} ${threads} ${DIR_WORK} ${top} ${cores} ${mode}
  done
done

```

### Interpreting the output

The final file `TOTAL_H2_<trait>_exome.txt` is tab-delimited with one row per run:

| Column | Meaning |
|--------|---------|
| `TRAIT` | Trait name |
| `N` | Mean number of participants across blocks |
| `N_RVs` | Total number of rare variants included |
| `ADJ_R2` | **Total heritability estimate** (sum of block adjusted R²) |
| `STD2` | Aggregated standard deviation (√ of summed block variances) |
| `LCL_adj` | Lower 95% CI bound (`ADJ_R2 − 1.96 × STD2`) |
| `UCL_adj` | Upper 95% CI bound (`ADJ_R2 + 1.96 × STD2`) |


**Heritability enrichment.** To assess a prioritization method, run the pipeline for
several proportion bins (e.g. `top1`, `top5`, `top100`) and compare `ADJ_R2`. If a
small top bin captures a large share of the `top100` heritability, the ranking
concentrates heritability among its highest-scored variants.

---

<!-- Resources -->

# Resources

RovHer scores for all 4,927,152 rare variants analyzed in this study from the UK Biobank, as well as pre-computed scores for 79,971,228 possible autosomal rare variants in humans, will be made publicly available for download on [Zenodo](https://zenodo.org/records/15596103?preview=1) upon publication. 

A separate set of 79,971,228 scores, trained without prediction features subject to commercial licensing restrictions, will also be provided.

<!-- Acknowledgements -->
## Acknowledgements

We gratefully acknowledge and thank the authors of various in silico tools that we utilized in our study for making their pre-computed scores and training data readily available.

## References

[1] Milborrow, S. (2023, January 26). earth: Multivariate Adaptive Regression Splines

[2] Karczewski, K. J., Solomonson, M., Chao, K. R., Goodrich, J. K., Tiao, G., Lu, W., Riley-Gillis, B. M., Tsai, E. A., Kim, H. I., Zheng, X., Rahimov, F., Esmaeeli, S., Grundstad, A. J., Reppell, M., Waring, J., Jacob, H., Sexton, D., Bronson, P. G., Chen, X., … Neale, B. M. (2022). Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. Cell Genomics, 2(9), 100168. https://doi.org/10.1016/j.xgen.2022.100168

[3] Pathan, N., Deng, W. Q., Khan, M., Scipio, M. D., Mao, S., Morton, R. W., Lali, R., Pigeyre, M., Chong, M. R., & Paré, G. (2022). A method to estimate the contribution of rare coding variants to complex trait heritability. Nature Communications, 15(1), 1245. https://doi.org/10.1038/s41467-024-45407-8

<!-- Citing -->
## Citing RovHer

If you use RovHer in your research, please cite our paper (citation details will be added upon publication).

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

## Contact

- Maintainer: Keona Pang (email: pangk3@mcmaster.ca) 
- Lab: Genetic and Molecular Epidemiology Lab, McMaster University
- License: MIT

<p align="right">(<a href="#readme-top">back to top</a>)</p>

