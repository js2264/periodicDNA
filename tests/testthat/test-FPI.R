context("test-FPI")

test_that("FPI works", {
    expect_equal({
        data(ce_proms_seqs)
        g <- DNAStringSet2GRanges(ce_proms_seqs)
        rand_regions <- sampleGRanges(g, 100) %>% 
            resize(fix = 'center', width = 200) %>% 
            trim() %>% 
            '['(width(.) == 200)
        rand_regions_seq <- ce_proms_seqs[rand_regions]
        fpi <- FPI(
            rand_regions_seq,
            motif = 'TT'
        )
        p <- plotFPI(fpi)
        class(p)[[1]] == "gg"
    }, TRUE)
})
