context("test-FPI")

test_that("FPI works", {
    expect_equal({
        data(ce11_proms)
        rand_regions <- sampleGRanges(
            ce11_proms, n = 20, w = 300, exclude = FALSE
        )
        fpi <- getFPI(
            rand_regions[1:10],
            genome = 'ce11', 
            motif = 'TT', 
            n_shuffling = 3, 
            cores_shuffling = 1,
            order = 1
        )
        p <- plotFPI(fpi)
        TRUE
    }, TRUE)
})

test_that("FPI order 2 works", {
    skip('skip2')
    expect_equal({
        data(ce11_proms)
        rand_regions <- sampleGRanges(
            ce11_proms, n = 20, w = 300, exclude = FALSE
        )
        fpi <- getFPI(
            rand_regions[1:10],
            genome = 'ce11', 
            motif = 'TT', 
            n_shuffling = 3, 
            cores_shuffling = 1,
            order = 2
        )
        p <- plotFPI(fpi)
        TRUE
    }, TRUE)
})

test_that("FPI test", {
    skip('skip2')
    expect_equal({
        data(ce11_TSSs)
        fpi <- getFPI(
            ce11_TSSs[['Ubiq.']][1:500],
            genome = 'ce11', 
            motif = 'TT', 
            n_shuffling = 1000,
            cores_shuffling = 100
        )
        p <- plotFPI(fpi)
        TRUE
    }, TRUE)
})

test_that("FPI several organisms", {
    skip('skip2')
    expect_equal({
        genomes <- c('sacCer3', 'ce11', 'dm6', 'danRer10', 'mm10', 'hg38')
        rand_regions_1 <- mclapply(mc.cores = 6, genomes, function(genome) {
            set.seed(1)
            sampleGenome(
                genome, n = 2000, w = 1000, exclude = FALSE
            )
        }) %>% setNames(genomes)
        fpis_genomes <- mclapply(mc.cores = 6, genomes, function(genome) {
            getFPI(
                rand_regions_1[[genome]]$seq,
                motif = 'TT', 
                BPPARAM = setUpBPPARAM(20), 
                n_shuffling = 20
            )
        })
        p <- cowplot::plot_grid(plotlist = lapply(fpis_genomes, plotFPI), nrow = 2)
        ggsave('tmp_TT_1to200_seed1.pdf', width = 15, height = 10)
    }, TRUE)
})

test_that("FPI TSSs of several organisms", {
    skip('skip3')
    expect_equal({
        genomes <- c('sacCer3', 'ce11', 'dm6', 'danRer10', 'mm10', 'hg38')
        genome_seqs <- mclapply(mc.cores = 6, genomes, function(genome) {
            seqs <- Biostrings::getSeq(char2BSgenome(genome))
        }) %>% setNames(genomes)
        load('~/20191025_GSEA_organisms/.list_TSSs.RData')
        names(list_TSSs) <- genomes[2:6]
        TSSs <- lapply(genomes[2:6], function(genome) {
            list(
                'lowCV' = list_TSSs[[genome]][list_TSSs[[genome]]$cv_bin %in% c(1,2)],
                'highCV' = list_TSSs[[genome]][list_TSSs[[genome]]$cv_bin %in% c(9, 10)]
            )
        }) %>% setNames(genomes[2:6])
        #
        fpis_genomes <- mclapply(mc.cores = 6, genomes[2:6], function(genome) {
            set.seed(2)
            r1 <- getFPI(
                genome_seqs[[genome]][TSSs[[genome]][['lowCV']]],
                motif = 'WW', 
                BPPARAM = setUpBPPARAM(20), 
                n_shuffling = 20, 
                range_spectrum = 1:75
            )
            set.seed(2)
            r2 <- getFPI(
                genome_seqs[[genome]][TSSs[[genome]][['highCV']]],
                motif = 'WW', 
                BPPARAM = setUpBPPARAM(20), 
                n_shuffling = 20, 
                range_spectrum = 1:75
            )
            return(list(r1, r2))
        })
        # fpis_genomes <- purrr::flatten(fpis_genomes)
        p <- cowplot::plot_grid(plotlist = lapply(fpis_genomes, plotFPI, s = 0.1), nrow = 3)
        ggsave('TSSs_TT-periodicity_1to75.pdf', width = 15, height = 10)
    }, TRUE)
})
