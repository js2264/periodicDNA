context("test-getPeriodicity")

test_that("getPeriodicity works", {
    expect_equal({
        data(ce_proms)
        periodicity_result <- getPeriodicity(
            ce_proms[1:2],
            genome = 'ce11',
            motif = 'TT', 
            cores = 1
        )
        list_plots <- plotPeriodicityResults(periodicity_result)
        class(list_plots) == "list"
    }, TRUE)
})
