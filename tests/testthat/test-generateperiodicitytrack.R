context("test-generateperiodicitytrack")

test_that("generateperiodicitytrack works", {
    expect_equal({
        data(ce_proms)
        generatePeriodicityTrack(
            Biostrings::getSeq(
                BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
            ),
            granges = ce_proms[1], 
            MOTIF = 'TT',
            FREQ = 1/10,
            PROCS = 1, 
            GENOME.WINDOW.SIZE = 100, 
            GENOME.WINDOW.SLIDING = 20, # can be 1 for single-base resolution
            BIN.WINDOW.SIZE = 60, # Set BIN.WINDOW.SIZE == GENOME.WINDOW.SIZE for no sliding window
            BIN.WINDOW.SLIDING = 5, 
            bw.file = 'TT-10-bp-periodicity_over-proms_gwin100_bwin60_bslide5.bw'
        )
        track <- rtracklayer::import('TT-10-bp-periodicity_over-proms_gwin100_bwin60_bslide5.bw', as = 'Rle')
        scaled_track <- scaleBigWigs(track)
        vec <- na.replace(scaled_track, 0)
        vec <- na.remove(scaled_track)
        p <- plotAggregateCoverage(
            track, 
            ce_proms
        )
        unlink('TT-10-bp-periodicity_over-proms_gwin100_bwin60_bslide5.bw')
        any(class(p) == "gg")
    }, TRUE)
})
