library(testthat)
library(periodicDNA)
library(BiocParallel)
register(setUpBPPARAM(1), default = TRUE)

test_check("periodicDNA")
