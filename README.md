# TisCoMM
TisCoMM leverages the co-regulation of genetic variations across different tissues explicitly via a unified probabilistic model. TisCoMM not only performs hypothesis testing to prioritize gene-trait associations black, but also detects the tissue-specific role of candidate target genes in complex traits. To make use of widely available GWAS summary statistics, we extend TisCoMM to use summary-level data, namely, TisCoMM-S$^2$. 

**TisCoMM** R package implements the TisCoMM method described in [Shi 2020](https://www.biorxiv.org/content/10.1101/789396v1). Get started with the [User Manual](https://github.com/XingjieShi/TisCoMM/blob/master/vignettes/TisCoMM.pdf). 

## Installation
To install the development version of **TisCoMM**, it's easiest to use the 'devtools' package. Note that **TisCoMM** depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```{r, fig.show='hold', eval=FALSE}
library(devtools)
install_github("XingjieShi/TisCoMM")
```
## Usage

```{r, fig.show='hold', eval=FALSE}
library(TisCoMM)
?TisCoMM
```
## Replicate simulation results in Shi et al. (2019)
All the simulation results can be reproduced by using the code in folder [simulation](https://github.com/XingjieShi/TisCoMM/tree/master/simulation). Before running simulation to reproduce the results, please familiarize yourself with **TisCoMM** using 'TisCoMM' vignette. 

1. Simulation results for multi-tissue joint test can be reproduced by following steps:

    - ExampleOne.R: This function can be run in a HPC cluster (with minor revisions, it could be run on a PC), it will output files, named pvalue_hz0.1_hc0.25_rhoX5_s5_batch-6.txt, which contain inference results of each replicate, for all multi-tissue TWAS  methods: TisCoMM, TisCoMM-S$^2$, MultiXcan, S-MultiXcan and UTMOST. 


    - ExampleOnePlot.R: This function produces simulation figures of joint test in Shi et al. (2019).

2. Simulation results for tissue-specific test can be reproduced by following steps:

     - PartCoMMCOR.R: This function can be run in a HPC cluster (with minor revisions, it could be run on a PC), it will output files, named part_hc4_rhoX8_rhoW8nz_ti2_batch-2.rds, which contain inference results of each replicate, for all single-tissue TWAS methods: CoMM, PrediXcan, and TWAS. 
     
    - SummaryCOR.R: This function produces simulation figures of tissue-specific test in Shi et al. (2019).
 

## Reference
[A tissue-specific collaborative mixed model for jointly analyzing multiple tissues in transcriptome-wide association studies](https://www.biorxiv.org/content/10.1101/789396v1) 
