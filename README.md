# RovHer: heritability-optimized scores for the functional prioritization of rare missense variants
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21908110.svg)](https://doi.org/10.5281/zenodo.21908110)

<!-- TABLE OF CONTENTS -->
<a name="readme-top"></a>
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about">About</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#hardware">Hardware</a></li>
        <li><a href="#software-dependencies">Software Dependencies</a></li>
      </ul>
    </li>
    <li><a href="#usage-retrieve-pre-computed-scores">Usage: Retrieve Pre-computed Scores</a></li>
    <li>
      <a href="#benchmark-rare-variant-heritability-with-rarity">Benchmark: RV Heritability with RARity</a>
      <ul>
        <li><a href="#prepare-input-files">Prepare Input Files</a></li>
        <li><a href="#run-the-pipeline">Run the Pipeline</a></li>
      </ul>
    </li>
    <li><a href="#acknowledgements">Acknowledgments</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#citing-rovher">Citing RovHer</a></li>
    <li><a href="#references">References</a></li>
  </ol>
</details>


<!-- ABOUT -->
## About

Predicting **rare variants** (RVs; MAF < 1%) that influence complex disease risk is a major challenge in human genetics. We introduce **RovHer** (RV heritability-optimized scores), an unbiased, scalable method that scores missense RVs by their probability of functional effect. RovHer uses a [Multivariate Adaptive Regression Splines](https://CRAN.R-project.org/package=earth) [1] model to integrate feature annotations with [Genebass](https://app.genebass.org/) [2] (N = 394,841) exome-wide association study (ExWAS) summary statistics of height, the training trait. As the dependent variable, we used the **false discovery rate** across 4,927,152 RVs as a surrogate for the likelihood of variant functionality.

The raw RovHer output reflects each variant's predicted FDR. We invert these to a "true positive" scale, so higher RovHer scores (closer to 1) indicate a greater probability of a variant being functional, and lower scores (closer to 0) indicate likely neutral variants.

**Performance**
* Prioritizes missense RVs that explain the greatest proportion of genome-wide trait variance — the variants most likely to be functional [3].

**Usage**:
* Scores are trait-agnostic (not tied to a specific disease or phenotype).
* Scores were generated for all possible rare SNVs in the dbNSFP database, but are optimized for missense variants.

![Workflow Overview](RovHer%20workflow.png)
*An overview of the development and application of RovHer.*

# Getting started
### Hardware
RovHer runs on all major operating systems (GNU/Linux, macOS, Windows). For biobank-scale analyses, we recommend Unix-based hardware with at least 100 GB RAM.

### Software dependencies

The development version is tested on Linux, but the required CRAN packages are compatible with Windows, macOS, and Linux. In your R console:
   ```R
  install.packages("data.table")
  install.packages("R.utils")
  install.packages("dplyr")
  install.packages("tidyverse")
  install.packages("corpcor") # pseudoinverse fallback
  install.packages("MBESS") # confidence intervals on R²

  ```
`lmutils` is **optional** and only needed for `mode = "lmutils"` in Step 4 (a faster, Rust-based implementation). See: https://github.com/GMELab/lmutils.r

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

<!-- Usage: Retrieve pre-computed scores -->
# Usage: Retrieve pre-computed scores

1. Clone the repo into your working directory (a few seconds on any computer):
   ```sh
   DIR_WORK="/your/repo/dir"
   cd $DIR_WORK

   git clone https://github.com/GMELab/RovHer.git
   cd RovHer
   ```

2. Download `All_RovHer_Scores.txt.gz` from [Zenodo](https://doi.org/10.5281/zenodo.21908110) and save it to your local repo as `RovHer/All_RovHer_Scores.txt.gz` (up to several minutes).

3. Run one of two helper scripts to extract a subset of scores:

   | Option | Description | Output columns | Script |
   |:--:|-------------|----------------|--------|
   | A | Retrieve scores for a list of RVs in `variants.txt` | `PLINK_SNP_NAME`, `Gene`, `RovHer_score` | `get_scores.R` |
   | B | Retrieve scored RVs for one or more protein-coding genes | `PLINK_SNP_NAME`, `Gene`, `RovHer_score` | `get_scores_per_gene.R` |

**Option A** takes a single text file of PLINK IDs (no header, unsorted is fine):

```
1:10030:A:T
8:203440:G:C
```

```sh
INFILE="Demo/variants.txt"   # input list
DIR_OUT="./Demo"             # output directory

cd $DIR_WORK/RovHer
Rscript Scripts/get_scores.R $INFILE $DIR_OUT
```

**Option B** takes no input file — just specify one or more genes in uppercase:

```sh
DIR_OUT="./Demo"             # output directory

GENE="LDLR BRCA1 APOB"       # a set of genes
# or
GENE="LDLR"                  # a single gene

Rscript Scripts/get_scores_per_gene.R "$DIR_OUT" "$GENE"
```

Runtime ranges from seconds to a few minutes depending on input size, on a desktop (16+ GB RAM, 4+ cores).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

# Benchmark: Rare variant heritability estimation with RARity 

Given a ranked list of variants (e.g. prioritized by the RovHer annotation score), this pipeline estimates how much trait heritability is captured by a chosen subset of variants.  **RARity** estimates heritability from *blocks* of approximately LD-independent variants. Within each block, the normalized phenotype is regressed on the block's genotype matrix, and the adjusted R² of that regression is an estimate of the heritability explained by the block. Block estimates (and their variances) are then summed into a single trait-level heritability with a 95% confidence interval. Block estimates are summed into a single trait-level heritability with a 95% confidence interval.

The pipeline runs as four sequential scripts:

| Script | Purpose | Outputs |
|--------|---------|---------|
|`1_split_prop_blks.R` | Split a clumped variant list into top-proportion subsets | example: [sample/clumped_plinkids_top1.clumped](sample/clumped_plinkids_top1.clumped) |
|`2_geno_blks.R` | Extract genotype blocks for a target variant subset | Filtered genotype blocks (`GENE_DF` objects) written to `<DIR>/h2_result/top<top>/2_GENO_RDATA/`, one set of blocks per proportion bin. example: [sample/h2_result/top1/2_GENO_RDATA](sample/h2_result/top1/2_GENO_RDATA) | 
|`3_align_geno_pheno.R` | Align both genotype and phenotype matrices by participant `eid` and processes them (i.e. mean-impute, standardize) | directories of genotype blocks (example: [3_GENO_aligned_LDL_direct](sample/h2_result/top1/3_GENO_aligned_LDL_direct)) and phenotype files (example: [3_PHENO_aligned_LDL_direct](sample/h2_result/top1/3_PHENO_aligned_LDL_direct))|
|`4_exome_h2_RV.R` | Estimate per-block heritability then aggregate into one trait-level estimate| `TOTAL_H2_<trait>_exome.txt` in `4_<trait>_H2_RESULTS/` example: [sample/h2_result/top1/4_LDL_direct_H2_RESULTS](sample/h2_result/top1/4_LDL_direct_H2_RESULTS)|

 
By running the pipeline on different **variant subsets** (as defined by variable `prop_blks` below), for example the top *X%* of variants ranked by a functional annotation — you can measure how heritability concentrates among prioritized variants. RARity is both as a standalone RV heritability estimator and as an **evaluation metric for variant-prioritization methods** (i.e. does ranking variants by a given score concentrate heritability at the top?). 

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Prepare input files

Example versions of every input are in the [`sample/`](sample/) directory. All text files are **tab-delimited**. In the tables below, **Required** marks the columns the pipeline looks up *by name* — any additional columns are ignored. Variant IDs use the `chr:pos:ref:alt` format (underscore-delimited `chr_pos_ref_alt` is also accepted).

**1. Clumped variant list** — `input_file` &nbsp; ([sample/clumped_plinkids.txt](sample/clumped_plinkids.txt))
Approximately LD-independent variants, typically from PLINK `--clump`.
* Plain text, **no header**, one variant ID per line.

```
1:10030:A:T
8:203440:G:C
...
```

**2. Master score file** — `score_file` &nbsp; ([sample/master_score_file.txt](sample/master_score_file.txt))
The annotation scores used to rank variants.

| PLINK_SNP_NAME | `<anno_name>` |
|--------|:--------:|
| `chr:pos:ref:alt` | 0.80|
| ... | ...|

**3. Covariates / principal components** — `PC_file` &nbsp; ([sample/PCs_1_20.txt](sample/PCs_1_20.txt))

| IID | AGE | SEX | PC1 | ... | PC20 | 
|--------|:--------:|:--------:|:--------:|:--------:|:--------:|
| p1 | 68| 0 | 0.00043 | ... | 0.003|
| ... | ...| ... | ...| ... | ...|

**4. Phenotype file** — `pheno_file` &nbsp; ([sample/pheno_file.txt](sample/pheno_file.txt))

| IID | `<trait>` |
|--------|:--------:|
| p1 | 6.50|
| ... | ...|

**5. Genotype directory** — `GENOTYPE_DIR` &nbsp; ([sample/clumped_geno_matrices](sample/clumped_geno_matrices))
A directory of per-chromosome genotype matrices, each an `.RData` file containing a **single object named `GENE_DF`**.
* File naming: `CHR_<chr>.clumped_<blk>.RData` — `<chr>` = 1–22, `<blk>` = 00, 01, 02, …
* `GENE_DF` layout: **first column** = participant ID (`IID`/`eid`); **remaining columns** = variant genotypes, named `chr:pos:ref:alt`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Run the pipeline

A copy of the wrapper script (as described below) to execute the four sequential scripts is found in: [Scripts/heritability_wrapper.sh](Scripts/heritability_wrapper.sh)

**Directory**:
* `DIR_WORK` is the full path to the repository

**Variables**:
* `anno_name` is the score column to rank by (default: `RovHer`) 
* `prop_blks` is a list of comma-separated proportion bins (e.g. `1` for top 1% of RVs, or `100` for all RVs)
* `trait` is the trait name (e.g. `LDL_direct`)
* `sort` is how to sort list of variants; `descending` (higher score = more functional) or `ascending`
* `mode` is the step 4 backend: `regular` (R) or `lmutils` (Rust) — default `regular`

**Prepared input files**:
* `input_file` path to a list of approximately LD-independent variants (formatted as chr:pos:ref:alt) typically generated from the PLINK
   `--clump` function; no column header
* `score_file` path to the master annotation file. Retrieve pre-computed RovHer annotations from: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21908110.svg)](https://doi.org/10.5281/zenodo.15596103)
* `GENOTYPE_DIR` directory containing per-chromosome genotype matrices for these LD-independent variants saved as `.RData` files each containing a
   single object named `GENE_DF`. Here, variant IDs are in `chr:pos:ref:alt_alt` format.
* `pheno_file` path to a participant trait file (first column IID, second column is the `trait`)
* `PC_file` path to the first 20 genetic principle components for each participant (columns: `eid`, `AGE`, `SEX`, `PC1`…`PC20`)
  
```bash

### Directories

DIR_WORK="/your/repo/dir"
cd ${DIR_WORK}/RovHer

### Settings

anno_name="RovHer"
prop_blks=(1,5,10,15,20,25,30,35,40,50,60,70,75,80,85,90,95)
trait="LDL_direct"
mode="regular" 
sort="descending"

### Input files / directories

input_file="./sample/clumped_plinkids.txt"
score_file="./sample/master_score_file.txt"
GENOTYPE_DIR="./sample/clumped_geno_matrices" # directory
pheno_file="./sample/pheno_file.txt"
PC_file="./sample/PCs_1_20.txt"

### Run

Rscript "./Scripts/1_split_prop_blks.r" ${input_file} ${anno_name} ${prop_blks} ${sort} ${score_file} ${DIR_WORK}
Rscript "./Scripts/2_build_geno_blks.r" ${DIR_WORK} ${input_file} ${prop_blks} ${GENOTYPE_DIR}

for top in 1 5 10 15 20 25 30 35 40 50 60 70 75 80 85 90 95; do
  cores=2
  threads=2
  Rscript "./Scripts/3_align_geno_pheno.r" ${anno_name} ${trait} ${DIR_WORK} ${top} ${cores} ${pheno_file} ${PC_file}
  Rscript "./Scripts/4_exome_wide_h2.r" ${anno_name} ${trait} ${threads} ${DIR_WORK} ${top} ${cores} ${mode}
done

```

### Interpreting the output

The final heritability estimate for each proportion bin is written to `h2_result/top<top>/4_<trait>_H2_RESULTS/TOTAL_H2_<trait>_exome.txt`
`h2_result/top<top>/4_<trait>_H2_RESULTS/TOTAL_H2_<trait>_exome.txt`. Example: [sample/h2_result/top1/4_LDL_direct_H2_RESULTS/TOTAL_H2_LDL_direct_exome.txt](sample/h2_result/top1/4_LDL_direct_H2_RESULTS/TOTAL_H2_LDL_direct_exome.txt)

(columns: `TRAIT`, `N`, `N_RVs`, `ADJ_R2` = heritability, `STD2`, `LCL_adj`, `UCL_adj` = 95% CI).


<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Directory structure

All outputs are written under `<DIR_WORK>/h2_result/top<top>/`:

```
<DIR_WORK>/h2_result/
└── top<top>/
    ├── 2_GENO_RDATA/                           # Step 2: filtered genotype blocks
    ├── 3_GENO_aligned_<trait>/                 # Step 3: processed, aligned genotype blocks
    ├── 3_PHENO_aligned_<trait>/                # Step 3: processed, aligned phenotype (norm_df)
    └── 4_<trait>_H2_RESULTS/                   # Step 4: heritability results
        ├── <trait>_H2_exome_raw.txt            # per-block raw R² (intermediate)
        ├── <trait>_H2_exome_raw_processed.txt  # per-block, with CIs/variances
        └── TOTAL_H2_<trait>_exome.txt          # ← final trait-level heritability
```

### Assess heritability enrichment

To assess a prioritization method, run the pipeline for several proportion bins (e.g. `top1`, `top5`, `top100`) and compare `ADJ_R2`. If a small top bin captures a large share of the `top100` heritability, the ranking concentrates heritability among its highest-scored variants.

---

<!-- Resources -->

# Resources

RovHer scores for all 4,927,152 rare variants analyzed in this study from the UK Biobank, as well as pre-computed scores for 79,971,228 possible autosomal rare variants in humans, will be made publicly available for download on [Zenodo](https://zenodo.org/records/15596103?preview=1) upon publication. A separate set of 79,971,228 scores, trained without prediction features subject to commercial licensing restrictions, will also be provided.

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

