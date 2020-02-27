# periodicDNA

![](man/images/TT_tissue-specific-classes.png)
![](man/images/WW-TT-AA-10bp-periodicity_tissue-spe-TSSs.png)

## Introduction

This R package helps the user to identify very short sequences (e.g. di- or 
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

## Overview

To begin, a genome sequence and a set of genomic loci must be defined. Let's 
focus on the C. elegans genome for now, and more specifically around its TSSs. 

```r
require(magrittr)
require(GenomicRanges)
require(ggplot2)
ce_seq <- Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11)
ce_proms <- readRDS(url('http://ahringerlab.com/VplotR/promoters_Granges.rds'))
```

## Dinucleotide periodicity over a set of Genomic Ranges

Let's look at the TT 10-bp periodicity strength around ubiquitous promoters:

```r
MOTIF <- 'TT'
ubiq_TT <- getPeriodicity(
    alignToTSS(ce_proms[ce_proms$which.tissues == 'Ubiq.'], 30, 500), 
    genome = ce_seq, 
    motif = MOTIF, 
    cores = 4
)
list_plots <- plotPeriodicityResults(ubiq_TT)
``` 

![TT-periodicity](man/images/ubiquitous-promoters_TT-periodicity.png)

## Make track of periodicity over a set of Genomic Ranges

The other major use of this package is to generate specific tracks 
over a set of loci, e.g. the strength of WW 10-pb periodicity over promoters.  
**Important note:** We recommand to run this command across at least a dozen of
processors (use the `PROCS` argument). This command will take several hours and
possibly up to a day to run. It typically takes one day to produce a periodicity
track over 15,000 GRanges of 150 bp (with default parameters) using `PROCS = 12`.
We highly recommand the user to run this command in a new `screen` session. 

```r
generatePeriodicityTrack(
    ce_seq,
    granges = ce_proms, 
    MOTIF = 'TT',
    FREQ = 1/10,
    PROCS = 12, 
    GENOME.WINDOW.SIZE = 100, 
    GENOME.WINDOW.SLIDING = 2, # can be 1 for single-base resolution
    BIN.WINDOW.SIZE = 60, # Set BIN.WINDOW.SIZE == GENOME.WINDOW.SIZE for no sliding window
    BIN.WINDOW.SLIDING = 5, 
    bw.file = 'TT-10-bp-periodicity_over-proms_gwin100_bwin60_bslide5.bw'
)
```

## Other functions of the package

Please read the [Introduction](vignettes/Introduction.md) vignette 
for a full presentation of the package functions.

