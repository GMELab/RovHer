# RovHer: heritability-optimized scores for the functional prioritization of rare missense variants

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

Predicting **rare variants** (RVs; MAF < 1%) that influence complex disease risk is a significant challenge. We introduce **RovHer** (RV heritability-optimized scores), an unbiased, scalable method that scores missense RVs based on their probability of functional effect. RovHer employs the [Multivariate Adaptive Regression Splines](https://CRAN.R-project.org/package=earth) [1] model to integrate feature annotations with [Genebass](https://app.genebass.org/) [2] (N=394,841) exome-wide association study (ExWAS) summary statistics of height, which serves as the training trait. The dependent variable for training is the per-variant **false discovery rate**, a surrogate measure for the likelihood of variant functionality.

**Performance**:
* prioritizes RVs that explain the largest proportion of genome-wide phenotypic variance, which are variants that more likely functional and disease-relevant [3]

**Usage**:
* scores are trait-agnostic (i.e. not tied to a specific disease or phenotype)
* RovHer scores were generated for all RVs; however, they are optimized for missense variants, which demonstrates significant enrichment in RV heritability across 21 traits tested in the UK Biobank.

![Workflow Overview](RovHer%20workflow.png)
*An overview of the development and application of RovHer.*

# Getting started
### Hardware
RovHer can generate predictions on major operating systems, including GNU/Linux, macOS, and Windows. For biobank-scale analyses, we recommend Unix-based hardware with a minimum of 100GB RAM.

### Software dependencies
R packages data.table and tidyverse (latest versions).
  ```R
  install.packages("data.table")
  install.packages("tidyverse")
  ```
<!-- Usage: Retrieve pre-computed scores -->
# Usage: Retrieve pre-computed scores 

1. Clone the repo into your working directory
   ```sh
   mydir="/my/working/dir" # modify 
   cd $mydir

   git clone https://github.com/Keonapang/RovHer.git
   cd RovHer
   ```
2. Download `All_RovHer_Scores.txt.gz` from [Zenodo](https://zenodo.org/records/15596103?preview=1) and place it in the same directory as this script.

3. Run script to obtain RovHer scores in two ways:

| Option | Description | Output | Script name |
|--:|-----------|-----------|-----------|
|  A| Retrieve scores for a list of RVs in `variants.txt` | A `output_variants.txt` with tab-delimited columns: PLINK_SNP_NAME, Gene, RovHer_score | `get_scores.r` |
|  B| For a given protein-coding `gene` or set of genes, retrieve a list of scored RVs | A `output_gene.txt` or `output_geneset.txt` with tab-delimited columns: PLINK_SNP_NAME, Gene, RovHer_score | `get_scores.r` | `get_scores_per_gene.r` |

For **option A**, the required input is a single text file `variants.txt` consisting of a column of PLINK IDs (no headers). Variants do not have to be sorted. For example:
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

For **option B**, specify gene(s) in all capital letters:
  ```sh
    DIR_OUT="$mydir/RovHer/Demo" # output directory 
    GENE="LDLR BRCA1 APOB" # example of a gene set  
    # or
    GENE="LDLR" # example of one gene 

    cd $mydir/RovHer
    Rscript Scripts/get_scores_per_gene.r "$DIR_OUT" "$GENE"
  ``` 

<p align="right">(<a href="#readme-top">back to top</a>)</p>

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

<p align="right">(<a href="#readme-top">back to top</a>)</p>
