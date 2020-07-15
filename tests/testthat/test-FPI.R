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
