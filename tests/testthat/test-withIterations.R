context("test-iterations")

test_that("getPeriodicity works with iterations", {
    expect_equal({
        data(ce11_proms)
        res <- getPeriodicity(
            ce11_proms[1:10],
            genome = 'BSgenome.Celegans.UCSC.ce11', 
            motif = 'TT', 
            n_shuffling = 3, 
            cores_shuffling = 1,
            order = 1, 
            range_spectrum = seq_len(100)
        )
        p <- plotPeriodicityResults(res)
        TRUE
    }, TRUE)
})
