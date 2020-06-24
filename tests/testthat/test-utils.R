context("test-utils")

test_that("shuffleSeq works", {
    skip('skip')
    expect_true({
        data(ce11_proms_seqs)
        shuffleSeq(ce11_proms_seqs)
        shuffleSeq(ce11_proms_seqs[[1]])
        TRUE
    })
})

test_that("namedListToLongFormat", {
    expect_true({
        data(ce11_proms)
        ce11_proms <- alignToTSS(deconvolveBidirectionalPromoters(
            ce11_proms
        ), 0, 200)
        l <- lapply(
            c('Ubiq.', 'Germline'), 
            function(TISSUE) {
                data.frame(
                    n = names(ce11_proms[ce11_proms$which.tissues == TISSUE]), 
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
    skip('skip')
    expect_true({
        g <- sampleGRanges(
            (BSgenome.Scerevisiae.UCSC.sacCer3::
            BSgenome.Scerevisiae.UCSC.sacCer3), 
            width = 100, 
            n = 1000
        )
        #
        g2 <- sampleGRanges(
            g, 
            w = 100
        )
        g2 <- sampleGRanges(
            g, 
            width = 10, 
            n = 10, 
            exclude = TRUE, 
            avoid_overlap = TRUE
        )
        #
        TRUE
    })
})
