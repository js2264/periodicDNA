context("test-getPeriodicity")

test_that("getPeriodicity and plotPeriodicityResults works", {
    expect_equal({
        data(ce11_proms_seqs)
        periodicity_result <- getPeriodicity(
            ce11_proms_seqs[[1]],
            motif = 'TT', 
            verbose = TRUE
        )
        #
        data(ce11_proms)
        periodicity_result <- getPeriodicity(
            ce11_proms[seq_len(2)],
            range_spectrum = 1:100,
            genome = 'ce11',
            motif = 'TT'
        )
        #
        list_plots <- plotPeriodicityResults(
            periodicity_result, filter_periods = FALSE, 
            grid = FALSE, axis = TRUE, ticks = TRUE, border = FALSE
        )
        list_plots <- plotPeriodicityResults(
            periodicity_result, filter_periods = FALSE
        )
        list_plots <- plotPeriodicityResults(
            periodicity_result, filter_periods = TRUE
        )
        list_plots <- plotPeriodicityResults(
            periodicity_result, filter_periods = TRUE
        )
        methods::is(list_plots, "gg")
    }, TRUE)
})

test_that("getPeriodicity works with shuffling", {
    expect_equal({
        data(ce11_TSSs)
        data(ce11_proms)
        periodicity_result <- getPeriodicity(
            ce11_proms[seq_len(10)],
            range_spectrum = 1:100,
            genome = 'ce11',
            motif = 'TT',
            n_shuffling = 3
        )
        periodicity_result_2 <- getPeriodicity(
            ce11_TSSs[['Ubiq.']][seq_len(10)] %>% 
                resize(width = 250, fix = 'center'),
            genome = 'ce11', 
            motif = 'TT', 
            n_shuffling = 1, 
            cores_shuffling = 1
        )
        periodicity_result_2 <- getPeriodicity(
            ce11_TSSs[['Ubiq.']][seq_len(10)] %>% 
                resize(width = 250, fix = 'center'),
            genome = 'ce11', 
            motif = 'TT', 
            n_shuffling = 4, 
            cores_shuffling = 1
        )
        list_plots <- plotPeriodicityResults(periodicity_result_2)
        # ggsave('tmp2.pdf', width = 12, height = 4)
        # getPeriodsMetrics(periodicity_result_2)$periodicityMetrics %>% dplyr::mutate(ObservedPSD = formatC(ObservedPSD, format = "e", digits = 2)) %>% dplyr::mutate(pval = formatC(pval, , format = "e", digits = 2)) %>% head(30) %>% knitr::kable()
        methods::is(list_plots, "gg")
    }, TRUE)
})
