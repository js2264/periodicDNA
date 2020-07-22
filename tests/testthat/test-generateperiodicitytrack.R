context("test-getPeriodicityTrack")

test_that("getPeriodicityTrack and plotAggregateCoverage works", {
    expect_equal({
        data(ce11_proms)
        track <- getPeriodicityTrack(
            genome = 'BSgenome.Celegans.UCSC.ce11',
            granges = ce11_proms[2:3], 
            extension = 100, 
            motif = 'TT',
            period = 10,
            window_size = 60, 
            step_size = 30,
            smooth_track = 1,
            bw_file = 'TT-10-bp-periodicity_over-proms.bw', 
            setUpBPPARAM(1)
        )
        scaled_track <- scaleBigWigs(list('test' = track))
        scaled_track2 <- scaleBigWigs(track)
        subtrack <- RleList('chrI' = scaled_track2[[1]][1:12000])
        subtrack <- smoothBigWig(subtrack, k = 10, setUpBPPARAM(1))
        vec <- na.replace(scaled_track, 0)
        vec <- na.remove(scaled_track)
        p <- plotAggregateCoverage(
            track, 
            ce11_proms[1:4]
        )
        q <- plotAggregateCoverage(
            list('test' = track, 'test2' = track), 
            list('g1' = ce11_proms[1:4], 'g2' = ce11_proms[1:4]), 
            split_by_granges = TRUE, split_by_track = TRUE
        )
        r <- plotAggregateCoverage(
            list('test' = track, 'test2' = track), 
            list(ce11_proms[1:4], ce11_proms[1:4]), 
            split_by_granges = FALSE, split_by_track = TRUE, 
            free_scales = FALSE
        )
        s <- plotAggregateCoverage(
            list('test' = track, 'test2' = track), 
            list('g1' = ce11_proms[1:4], 'g2' = ce11_proms[1:4]), 
            split_by_granges = TRUE, split_by_track = FALSE, 
            free_scales = TRUE
        )
        unlink('TT-10-bp-periodicity_over-proms.bw')
        methods::is(s, "gg")
    }, TRUE)
})
