library(testthat)
library(periodicDNA)
library(BiocParallel)
register(SnowParam(workers = 1), default = TRUE)

test_check("periodicDNA")
