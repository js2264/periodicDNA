context("test-utils")

test_that("shuffleSeq works", {
    expect_true({
        data(ce_proms_seqs)
        shuffleSeq(ce_proms_seqs)
        shuffleSeq(ce_proms_seqs[[1]])
        TRUE
    })
})

test_that("namedListToLongFormat", {
    expect_true({
        data(ce_proms)
        ce_proms <- alignToTSS(deconvolveBidirectionalPromoters(
            ce_proms
        ), 0, 200)
        l <- lapply(
            c('Ubiq.', 'Germline'), 
            function(TISSUE) {
                data.frame(
                    n = names(ce_proms[ce_proms$which.tissues == TISSUE]), 
                    t = 1
                )
            }
        )
        names(l) <- c('Ubiq.', 'Germline')
        d <- namedListToLongFormat(l)
        TRUE
    })
})

test_that("sampleGRanges works", {
    expect_true({
        data(ce_proms)
        g <- sampleGRanges(
            ce_proms, 
            100
        )
        g <- sampleGRanges(
            ce_proms, 
            width = 10, 
            n = 10
        )
        #
        g <- sampleGRanges(
            (BSgenome.Scerevisiae.UCSC.sacCer3::
            BSgenome.Scerevisiae.UCSC.sacCer3), 
            100, 
            exclude = FALSE
        )
        #
        data(ce_proms_seqs)
        g <- sampleGRanges(
            ce_proms_seqs,
            n = 10, 
            width = 10, 
            exclude = FALSE
        )
        TRUE
    })
})
