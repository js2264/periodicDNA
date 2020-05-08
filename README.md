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

periodicDNA can be installed from Github as follows:

```r
install.packages("devtools")
devtools::install_github("js2264/periodicDNA")
library(periodicDNA)
```

## Main functions 

The three main user-level functions of periodicDNA are `getPeriodicity()`, 
`getFPI()` and `getPeriodicityTrack()`. 

* `getPeriodicity()` is used to compute the power spectral density 
  (PSD) of a chosen k-mer (i.e. `TT`) in a set of sequences. The PSD 
  score at a given period indicates the strength of the k-mer at 
  this period. 
* `getPeriodicityTrack()` can be used to generate linear tracks representing 
  the periodicity strength of a given k-mer at a chosen period, over genomic
  loci of interest. 
* `getFPI()` is used to compute the Fold Power Increase, a more sophisticated 
  metric derived from the PSD. It was initially developed by Pich et al., 
  Cell 2018. It takes into account the background periodicity
  of the k-mer of interest in the provided sequences. 

### `getPeriodicity()` function

```r
data(proms)
PSDs <- getPeriodicity(
    proms,
    genome = 'ce11',
    motif = 'TT', 
    cores = 4
)
plotPeriodicityResults(PSDs)
```

### `getPeriodicityTrack()` function

```r
WW_10bp <- getPeriodicityTrack(
    genome = 'ce11',
    granges = proms, 
    motif = 'WW',
    period = 10,
    cores = 12, 
    bw_file = 'WW-10-bp-periodicity_over-proms.bw'
)
```

### `getFPI()` function

```r
ce_seq <- Biostrings::getSeq(
    BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
)
FPI <- getFPI(
    ce_seq[proms], 
    motif = 'TT', 
    cores_shuffling = 10, 
    n_shuffling = 10, 
    FUN = stats::spectrum
)
plotFPI(FPI)
```

## Contributions
Code contributions, bug reports, fixes and feature requests are most welcome.
Please make any pull requests against the master branch at 
https://github.com/js2264/periodicityDNA
and file issues at https://github.com/js2264/periodicityDNA/issues

## License 
**periodicDNA** is licensed under the GPL-3 license.
