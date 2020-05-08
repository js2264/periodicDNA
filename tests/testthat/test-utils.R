context("test-utils")

test_that("shuffleSeq works", {
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
    expect_true({
        data(ce11_proms)
        g <- sampleGRanges(
            ce11_proms, 
            100
        )
        g <- sampleGRanges(
            ce11_proms, 
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
        data(ce11_proms_seqs)
        g <- sampleGRanges(
            ce11_proms_seqs,
            n = 10, 
            width = 10, 
            exclude = FALSE
        )
        TRUE
    })
})
