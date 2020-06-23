[![](https://travis-ci.com/js2264/periodicDNA.svg?branch=master)](https://travis-ci.com/js2264/periodicDNA)
[![](https://codecov.io/gh/js2264/periodicDNA/branch/master/graph/badge.svg)](https://codecov.io/github/js2264/periodicDNA?branch=master)
[![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/github/languages/code-size/js2264/periodicDNA.svg)](https://github.com/js2264/periodicDNA)
[![](https://img.shields.io/badge/license-GPL--3-orange.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# periodicDNA <img src="man/figures/logo.png" align="right" alt="" />

![](https://raw.githubusercontent.com/js2264/periodicDNA/master/man/figures/TT_tissue-specific-classes.png)
![](https://raw.githubusercontent.com/js2264/periodicDNA/master/man/figures/WW-TT-AA-10bp-periodicity_tissue-spe-TSSs.png)

## Introduction

This R package helps the user identify k-mers (e.g. di- or 
tri-nucleotides) present periodically in a set of genomic loci (typically 
regulatory elements). It is not aimed at identifying motifs separated by a 
conserved distance; for this type of analysis, please visit 
[MEME](http://meme-suite.org) website.

## Installation

periodicDNA is available in Bioconductor. To install the current release use:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```

For advanced users, the most recent periodicDNA can be installed 
from Github as follows (might be buggy):

```r
install.packages("devtools")
devtools::install_github("js2264/periodicDNA")
library(periodicDNA)
```

## How to use periodicDNA

periodicDNA includes a vignette where its usage is 
illustrated. To access the vignette, please use:

```r
vignette('periodicDNA')
```

A description of the internal steps performed by periodicDNA is accessible 
in the following vignette:

```r
vignette('internal-steps', package = 'periodicDNA')
```

## Main functions 

The three main user-level functions of periodicDNA are `getPeriodicity()`,
`getPeriodicityTrack()` and `getFPI()`.

* `getPeriodicity()` is used to compute the power spectral density 
  (PSD) of a chosen k-mer (i.e. `TT`) in a set of sequences. The PSD 
  score at a given period indicates the strength of the k-mer at 
  this period. 
* `getPeriodicityTrack()` can be used to generate linear tracks representing 
  the periodicity strength of a given k-mer at a chosen period, over genomic
  loci of interest. 
* `getFPI()` is used to compute the Fold Power Increase, 
  a more sophisticated metric derived from the PSD. 
  It was initially developed by Pich et al., Cell 2018. 
  It estimates the background periodicity of the k-mer 
  of interest in shuffled sequences. 

### `getPeriodicity()` function

```r
data(ce11_proms)
PSDs <- getPeriodicity(
    ce11_proms,
    genome = 'ce11',
    motif = 'TT', 
    BPPARAM = BiocParallel::SnowParam(workers = 4)
)
plotPeriodicityResults(PSDs)
```

### `getPeriodicityTrack()` function

```r
WW_10bp <- getPeriodicityTrack(
    genome = 'ce11',
    granges = ce11_proms, 
    motif = 'WW',
    period = 10,
    bw_file = 'WW-10-bp-periodicity_over-proms.bw', 
    BPPARAM = BiocParallel::SnowParam(workers = 12)
)
```

### `getFPI()` function

```{r}
ce_seq <- Biostrings::getSeq(
    BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
)
FPI <- getFPI(
    ce_seq[ce11_proms], 
    period = 10,
    motif = 'TT', 
    cores_shuffling = 10, 
    n_shuffling = 10
)
plotFPI(FPI)
```

**Warning**: It is recommended to run this command across many processors 
using BiocParallel. This command typically takes one day to produce 
a periodicity track over 15,000 GRanges of 150 bp (with default parameters) 
using `BPPARAM = BiocParallel::SnowParam(workers = 12)`. 
It is highly recommended to run this command in a new `screen` session.

## Contributions
Code contributions, bug reports, fixes and feature requests are most welcome.
Please make any pull requests against the master branch at 
https://github.com/js2264/periodicityDNA
and file issues at https://github.com/js2264/periodicityDNA/issues

## License 
**periodicDNA** is licensed under the GPL-3 license.
