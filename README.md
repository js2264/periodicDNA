# perioDNA

**perioDNA is still in pre-Alpha and has not been thorouhly tested yet.  
Documentation will come soon.**

![perioDNA](examples/png/ubiquitous-promoters_TT-periodicity.png)

## Introduction

This package helps the user to identify very short sequences (e.g. di- or 
tri-nucleotides) present periodically in a set of genomic loci (typically 
regulatory elements). It is not aimed at identifying motifs separated by a 
conserved distance - for this analysis, please visit [MEME](http://meme-suite.org)
website.

## Installation

Most up-to-date version of perioDNA can be ran installed from Github as follow:

```r
install.packages("devtools")
devtools::install_github("js2264/perioDNA")
library(perioDNA)
```

## Overview

To begin, a genome sequence and a set of genomic loci must be defined. Let's 
focus on the C. elegans genome for now, and more specifically around its TSSs. 

```r
require(magrittr)
require(GenomicRanges)
require(ggplot2)
ce_seq <- Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11)
ce_proms <- readRDS(url('http://ahringerlab.com/perioDNA/ce11_annotated_REs.rds')) %>% 
    '['(.$is.prom) %>% 
    deconvolveBidirectionalPromoters() %>% 
    alignToTSS(50, 300)
```

## Dinucleotide periodicity over a set of Genomic Ranges

Let's look at the TT 10-bp periodicity strength around ubiquitous promoters:

```r
MOTIF <- 'TT'
ubiq_TT <- getPeriodicity(
    ce_proms[ce_proms$which.tissues == 'Ubiq.'], 
    genome = ce_seq, 
    motif = MOTIF, 
    cores = 4
)
list_plots <- plotPeriodicityResults(ubiq_TT)
``` 

![TT-periodicity](examples/png/ubiquitous-promoters_TT-periodicity.png)

## Make track of periodicity over a set of Genomic Ranges

Another major use of this package is to generate specific tracks 
over a set of loci, e.g. the strength of WW 10-pb periodicity over promoters.  
**Important note:** We recommand to run this command across multiple processors
(specific by the `PROCS` argument). This command will take several hours and
possibly up to a day to run. It would typically take one day to produce a periodicity
track over 15,000 GRanges of 150 bp (with default parameters) using `PROCS = 12`.
We highly recommand the user to run this command in a new `screen` session. 

```r
generatePeriodicityTrack(
    ce_seq,
    granges = ce_proms, 
    MOTIF = 'TT',
    FREQ = 1/10,
    PROCS = 100, 
    GENOME.WINDOW.SIZE = 100, 
    BIN.WINDOW.SIZE = 60, 
    BIN.WINDOW.SLIDING = 5, 
    bw.file = 'TT-10-bp-periodicity_over-proms_gwin100_bwin60_bslide5.bw'
)
```

## Other functions of the package

Please read the [Introduction](vignettes/Introduction.md) vignette 
for a full presentation of the package functions.

