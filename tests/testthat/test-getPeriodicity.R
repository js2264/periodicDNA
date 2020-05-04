context("test-getPeriodicity")

test_that("getPeriodicity works", {
    expect_equal({
        data(ce_proms)
        periodicity_result <- getPeriodicity(
            ce_proms[seq_len(2)],
            genome = 'ce11',
            motif = 'TT', 
            cores = 1
        )
        list_plots <- plotPeriodicityResults(periodicity_result)
        methods::is(list_plots, "list")
    }, TRUE)
})

test_that("getPeriodicity works with shuffling", {
    expect_equal({
        data(ce_proms)
        periodicity_result <- getPeriodicity(
            ce_proms[seq_len(2)],
            genome = 'ce11',
            motif = 'TT', 
            cores = 1, 
            skip_shuffling = FALSE, 
            doZscore = TRUE
        )
        list_plots <- plotPeriodicityResults(periodicity_result)
        methods::is(list_plots, "list")
    }, TRUE)
})
