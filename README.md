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

Predicting **rare variants** (RVs; MAF < 1%) that influence complex disease risk is a significant challenge. We introduce **RovHer** (RV heritability-optimized scores), an unbiased, scalable method that scores missense RVs based on their probability of functional effect. RovHer employs the [Multivariate Adaptive Regression Splines](https://CRAN.R-project.org/package=earth) [1] model to integrate feature annotations with [Genebass](https://app.genebass.org/) [2] (N=394,841) exome-wide association study (ExWAS) summary statistics of height, which serves as the training trait. For the dependent variable, we used the **false discovery rate** across 4,927,152 RVs as a surrogate measure for the likelihood of variant functionality.

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

The package development version is tested on Linux operating systems. However, the CRAN packages required should be compatible with Windows, Mac, and Linux operating systems. Ensure you have the following:
   ```R
  install.packages("data.table")
  install.packages("R.utils")
  install.packages("dplyr")
  install.packages("tidyverse")
  install.packages("corpcor") # pseudoinverse fallback
  install.packages("MBESS") # confidence intervals on R²

  ```

To install **`lmutils`** (optional), which is only required for `mode="lmutils"` in Step 4, type in R:
    ```r
    install.packages(
      "https://github.com/mrvillage/lmutils.r/archive/refs/heads/master.tar.gz",
      repos = NULL)   # use the .zip archive on Windows
    ```


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


# Rare variant heritability estimation with RARity 

A four-step pipeline for estimating **exome-wide rare-variant (RV) heritability** from UK Biobank–style genotype and phenotype data, using the RARity approach. **RARity** estimates heritability from *blocks* of approximately LD-independent variants. Within each block, the normalized phenotype is regressed on the block's genotype matrix, and the adjusted R² of that regression is an estimate of the heritability explained by the block. Block estimates (and their variances) are then summed into a single trait-level heritability with a 95% confidence interval.

Given a ranked list of variants (e.g. prioritized by the RovHer annotation score), the pipeline estimates how much trait heritability is captured by a chosen subset of variants. This enables a comparison of **heritability-enrichment** across different variant prioritization methods. The pipeline runs as sequential scripts:

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `1_split_prop_blks.R` | Splits a large variant list into top proportion bins |
| 2 | `2_geno_blks.R` | Extract genotype blocks for a target variant subset |
| 3 | `3_align_geno_pheno.R` | Align and process genotype + phenotype matrices |
| 4 | `4_exome_h2_RV.R` | Estimate per-block and total trait heritability |


By running the pipeline on different **variant subsets** — for example, the top *X%*
of variants ranked by a functional annotation — you can measure how heritability
concentrates among prioritized variants. This makes the pipeline useful both as a
standalone heritability estimator and as an **evaluation metric for variant-prioritization
methods** (does ranking variants by a given score concentrate heritability at the top?).

## Upstream dataset requirements

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

Here `<top>` is the proportion bin (e.g. `80`, or `100` for all variants) and
`<trait>` is the trait name (e.g. `LDL_direct`).

---

## Step 2 — Extract genotype blocks

**Script:** `2_geno_blks.R`

**Aim.** For each requested proportion bin, scan the full set of clumped genotype
matrices and keep only the variants present in that bin's variant list. Retained
columns are accumulated and written out as a smaller set of genotype blocks
(~5,000 variants each). Small subsets (< 1,000 variants) are written as a single
combined block spanning all chromosomes.

**Arguments** (positional):

| # | Name | Type | Description |
|---|------|------|-------------|
| 1 | `DIR` | path | Working directory root |
| 2 | `PLINK_list` | path | Base clumped-variant list (`.txt`); the per-bin file `<base>_top<top>.clumped` is loaded |
| 3 | `prop_blks` | string | Comma-separated proportion bins, e.g. `"1,5,10,80,100"` |
| 4 | `GENOTYPE_DIR` | path | Directory of clumped genotype matrices |

**Run:**

```bash
Rscript 2_geno_blks.R "$DIR" "$PLINK_list" "1,5,10" "$GENOTYPE_DIR"
```

**Output.** Filtered genotype blocks (`GENE_DF` objects) written to
`<DIR>/h2_result/top<top>/2_GENO_RDATA/`, one set per proportion bin.

---

## Step 3 — Align and process matrices

**Script:** `3_align_geno_pheno.R`

**Aim.** Prepare matched, model-ready genotype and phenotype matrices:

1. **Align** genotype and phenotype rows by participant `eid`.
2. **Process the phenotype:** quantile-normalize the trait, regress out `AGE`, `SEX`,
   and 20 principal components (keep residuals), then standardize to mean 0, sd 1.
3. **Process the genotype:** mean-impute missing values, drop columns with a minor
   allele count (MAC) < 2, and standardize each variant column to mean 0, sd 1.

**Arguments** (positional):

| # | Name | Type | Description |
|---|------|------|-------------|
| 1 | `anno` | string | Annotation label, e.g. `yhat`, `RovHer` |
| 2 | `trait` | string | Trait name, e.g. `LDL_direct` |
| 3 | `DIR_WORK` | path | Working directory root |
| 4 | `top` | string | Proportion bin, e.g. `80` or `100` |
| 5 | `cores` | int | Worker count for the parallel cluster |
| 6 | `chr_start` | int | First chromosome (per-chromosome mode) |
| 7 | `chr_end` | int | Last chromosome (per-chromosome mode) |
| 8 | `pheno_file` | path | Phenotype table (`eid` + trait) |
| 9 | `PC_file` | path | Age/sex/PC covariate table |

**Run:**

```bash
Rscript 3_align_geno_pheno.R RovHer LDL_direct "$DIR_WORK" 80 20 1 22 \
        "$pheno_file" "$PC_file"
```

**Output.**
- Processed, aligned genotype blocks → `3_GENO_aligned_<trait>/`
- Processed phenotype (`norm_df`, with and without `eid`) → `3_PHENO_aligned_<trait>/`

The phenotype step runs once and is skipped automatically on re-runs if its output
already exists. The script auto-detects whether to process a single combined
"`CHR_all`" block or per-chromosome blocks based on the input filenames.

---

## Step 4 — Estimate heritability

**Script:** `4_exome_h2_RV.R`

**Aim.** Estimate heritability block by block, then aggregate into one trait-level
estimate. For each block, the normalized phenotype is regressed on the block's
genotypes and the **adjusted R²** is taken as the block's heritability contribution;
block variances are combined to give a standard error and 95% confidence interval.

Two computation modes:
- `regular` — base-R implementation (no extra dependencies).
- `lmutils` — uses `lmutils::calculate_r2()` for the per-block step (faster; requires
  the `lmutils` package).

**Arguments** (positional):

| # | Name | Type | Description |
|---|------|------|-------------|
| 1 | `anno` | string | Annotation label |
| 2 | `trait` | string | Trait name |
| 3 | `threads` | int | Worker threads (physical cores available) |
| 4 | `DIR_WORK` | path | Working directory root |
| 5 | `top` | int | Proportion bin |
| 6 | `cores` | int | Core parallelism (files held in memory) |
| 7 | `mode` | string | `regular` or `lmutils` |

**Run:**

```bash
Rscript 4_exome_h2_RV.R RovHer LDL_direct 1 "$DIR_WORK" 80 2 regular
```

**Output.** `TOTAL_H2_<trait>_exome.txt` in `4_<trait>_H2_RESULTS/` (see next section).
On successful completion, the script removes the intermediate aligned genotype and
phenotype directories to save space.

---

## Interpreting the output

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

The per-block intermediate file `<trait>_H2_exome_raw_processed.txt` reports each
block's `r2`, `adj_r2`, sample size `n`, variant count `m`, per-variant adjusted R²,
and block-level variance terms — useful for diagnosing individual blocks.

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
- Lab: Genetic and Molecular Epidemiology Lab, McMaster University, Canada
- License: MIT

<p align="right">(<a href="#readme-top">back to top</a>)</p>


#####################################################################


## Pipeline flow

```
        clumped variant list (per proportion bin)         genotype matrices (chr 1–22)
                        │                                            │
                        └──────────────┬─────────────────────────────┘
                                       ▼
                     ┌─────────────────────────────────┐
   Step 2            │  2_geno_blks.R                   │
   Extract blocks    │  keep only target-subset RVs;    │
                     │  write ~5k-variant geno blocks   │
                     └─────────────────────────────────┘
                                       │  2_GENO_RDATA/
                                       ▼
   phenotype table ──▶┌─────────────────────────────────┐
   PC/covariate file │  3_align_geno_pheno.R            │
   Step 3            │  align by eid; quantile-norm +    │
   Process matrices  │  covariate-adjust + standardize   │
                     │  pheno; NA-impute + MAC≥2 + z-norm │
                     │  geno                             │
                     └─────────────────────────────────┘
                                       │  3_GENO_aligned_<trait>/
                                       │  3_PHENO_aligned_<trait>/
                                       ▼
                     ┌─────────────────────────────────┐
   Step 4            │  4_exome_h2_RV.R                 │
   Estimate h²       │  per-block adjusted R²;          │
                     │  aggregate → total h² + 95% CI    │
                     └─────────────────────────────────┘
                                       │
                                       ▼
                        TOTAL_H2_<trait>_exome.txt
```


