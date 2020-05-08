context("test-getPeriodicityTrack")

test_that("getPeriodicityTrack and plotAggregateCoverage works", {
    expect_equal({
        data(ce11_proms)
        track <- getPeriodicityTrack(
            Biostrings::getSeq(
                BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
            ),
            granges = ce11_proms[1], 
            motif = 'TT',
            period = 10,
            extension = 400, 
            genome_sliding_sliding = 5, 
            cores = 1, 
            bw_file = 'TT-10-bp-periodicity_over-proms.bw'
        )
        scaled_track <- scaleBigWigs(list('test' = track))
        scaled_track2 <- scaleBigWigs(track)
        vec <- na.replace(scaled_track, 0)
        vec <- na.remove(scaled_track)
        p <- plotAggregateCoverage(
            track, 
            ce11_proms
        )
        q <- plotAggregateCoverage(
            list('test' = track, 'test2' = track), 
            list('g1' = ce11_proms, 'g2' = ce11_proms), 
            split_by_granges = TRUE, split_by_track = TRUE
        )
        r <- plotAggregateCoverage(
            list('test' = track, 'test2' = track), 
            list(ce11_proms, ce11_proms), 
            split_by_granges = FALSE, split_by_track = TRUE, 
            free_scales = FALSE
        )
        s <- plotAggregateCoverage(
            list('test' = track, 'test2' = track), 
            list('g1' = ce11_proms, 'g2' = ce11_proms), 
            split_by_granges = TRUE, split_by_track = FALSE, 
            free_scales = TRUE
        )
        unlink('TT-10-bp-periodicity_over-proms.bw')
        methods::is(s, "gg")
    }, TRUE)
})

test_that("getPeriodicityTrack works", {
    skip('skip')
    expect_equal({
        load('~/shared/data/classification_tissue-spe-genes-REs_REs-GRanges.RData')
        ce_seq <- getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11)
        list_params <- mclapply(
            c('Ubiq.', 'Germline', 'Neurons', 'Muscle', 'Hypod.', 'Intest.'), 
            mc.cores = 6, 
            function(TISSUE) {
                proms_center <- granges.all[
                    granges.all$is.prom & granges.all$which.tissues == TISSUE
                ] %>% 
                    deconvolveBidirectionalPromoters() %>% 
                    resize(140, fix = 'center')
                proms_flanking <- c(
                    proms_center %>% resize(1, fix = 'center') %>% 
                        shift(140) %>% resize(140, fix = 'center'),
                    proms_center %>% resize(1, fix = 'center') %>% 
                        shift(-140) %>% resize(140, fix = 'center')
                )
                proms_extended <- c(
                    proms_center %>% resize(1, fix = 'center') %>% 
                        shift(280) %>% resize(140, fix = 'center'),
                    proms_center %>% resize(1, fix = 'center') %>% 
                        shift(-280) %>% resize(140, fix = 'center')
                )
                enhs_center <- granges.all[
                    granges.all$regulatory_class == 'putative_enhancer' & 
                    granges.all$which.tissues == TISSUE
                ] %>% 
                    deconvolveBidirectionalPromoters() %>% 
                    resize(100, fix = 'center')
                enhs_flanking <- c(
                    enhs_center %>% resize(1, fix = 'center') %>% 
                        shift(140) %>% resize(140, fix = 'center'),
                    enhs_center %>% resize(1, fix = 'center') %>% 
                        shift(-140) %>% resize(140, fix = 'center')
                )
                enhs_extended <- c(
                    enhs_center %>% resize(1, fix = 'center') %>% 
                        shift(280) %>% resize(140, fix = 'center'),
                    enhs_center %>% resize(1, fix = 'center') %>% 
                        shift(-280) %>% resize(140, fix = 'center')
                )
                l <- list(
                    proms_center,
                    proms_flanking,
                    proms_extended,
                    enhs_center,
                    enhs_flanking,
                    enhs_extended
                )
            }
        ) %>% setNames(c('Ubiq.', 'Germline', 'Neurons', 'Muscle', 'Hypod.', 'Intest.'))
        TTs <- rtracklayer::import(
            'TT-10-bp-periodicity_over-proms.bw', as = 'Rle'
        )
        #
        vec <- na.replace(scaled_track, 0)
        vec <- na.remove(scaled_track)
        p <- plotAggregateCoverage(
            track, 
            ce11_proms
        )
        p <- plotAggregateCoverage(
            list('test' = track), 
            ce11_proms
        )
        unlink('TT-10-bp-periodicity_over-proms.bw')
        methods::is(p, "gg")
    }, TRUE)
})
