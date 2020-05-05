context("test-generateperiodicitytrack")

test_that("generateperiodicitytrack works", {
    expect_equal({
        data(ce_proms)
        generatePeriodicityTrack(
            Biostrings::getSeq(
                BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
            ),
            granges = ce_proms[1], 
            motif = 'TT',
            freq = 1/10,
            cores = 1, 
            bw_file = 'TT-10-bp-periodicity_over-proms.bw'
        )
        track <- rtracklayer::import(
            'TT-10-bp-periodicity_over-proms.bw', as = 'Rle'
        )
        scaled_track <- scaleBigWigs(list('test' = track))
        scaled_track <- scaleBigWigs(track)
        vec <- na.replace(scaled_track, 0)
        vec <- na.remove(scaled_track)
        p <- plotAggregateCoverage(
            track, 
            ce_proms
        )
        p <- plotAggregateCoverage(
            list('test' = track), 
            ce_proms
        )
        unlink('TT-10-bp-periodicity_over-proms.bw')
        methods::is(p, "gg")
    }, TRUE)
})
