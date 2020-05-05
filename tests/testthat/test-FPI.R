context("test-FPI")

test_that("FPI works", {
    expect_equal({
        data(ce_proms_seqs)
        g <- DNAStringSet2GRanges(ce_proms_seqs)
        rand_regions <- sampleGRanges(
            g, n = 20, w = 200, exclude = FALSE
        )
        rand_regions_seq <- ce_proms_seqs[rand_regions]
        fpi <- getFPI(
            rand_regions_seq,
            motif = 'TT', 
            parallel_shuffling = 1
        )
        p <- plotFPI(fpi)
        methods::is(p, 'gg')
    }, TRUE)
})
