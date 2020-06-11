context("test-getPeriodicity")

test_that("getPeriodicity and plotPeriodicityResults works", {
    expect_equal({
        data(ce11_proms_seqs)
        periodicity_result <- getPeriodicity(
            ce11_proms_seqs[[1]],
            motif = 'TT', 
            cores = 1
        )
        #
        data(ce11_proms)
        periodicity_result <- getPeriodicity(
            ce11_proms[seq_len(2)],
            genome = 'ce11',
            motif = 'TT', 
            cores = 1
        )
        #
        list_plots <- plotPeriodicityResults(
            periodicity_result, skip_shuffling = FALSE, filter_periods = FALSE, 
            grid = FALSE, axis = TRUE, ticks = TRUE, border = FALSE
        )
        list_plots <- plotPeriodicityResults(
            periodicity_result, skip_shuffling = TRUE, filter_periods = FALSE
        )
        list_plots <- plotPeriodicityResults(
            periodicity_result, skip_shuffling = FALSE, filter_periods = TRUE
        )
        list_plots <- plotPeriodicityResults(
            periodicity_result, skip_shuffling = TRUE, filter_periods = TRUE
        )
        methods::is(list_plots, "gg")
    }, TRUE)
})

test_that("getPeriodicity works with shuffling", {
    expect_equal({
        data(ce11_proms)
        periodicity_result <- getPeriodicity(
            ce11_proms[seq_len(50)],
            genome = 'ce11',
            motif = 'TT', 
            cores = 1, 
            skip_shuffling = TRUE, 
            doZscore = TRUE
        )
        list_plots <- plotPeriodicityResults(periodicity_result)
        methods::is(list_plots, "gg")
    }, TRUE)
})

test_that("getPeriodicity on lists of proms and mutlitple dinucleotides", {
    skip('skip')
    expect_equal({
        data(ce11_proms)
        dinucs <- c('WW', 'SS', 'RR', 'YY', 'KK')
        periodicity_results <- mclapply(mc.cores = length(dinucs), dinucs, function(MOTIF) {
            mclapply(mc.cores = 6, c('Ubiq.', 'Germline', 'Muscle'), function(TISSUE) {
                message(MOTIF, ' -- ', TISSUE)
                granges <- alignToTSS(ce11_proms[ce11_proms$which.tissues == TISSUE], 30, 170)
                r <- getPeriodicity(
                    granges, 
                    genome = 'ce11', 
                    motif = MOTIF, 
                    range_spectrum = 1:100,
                    verbose = FALSE
                )
                # res <- r$PSD %>% 
                #     dplyr::slice(which.min(abs(period - 10))) %$% 
                #     PSD
                d <- data.frame(
                    PSD = res, 
                    tissue = TISSUE, 
                    dinuc = MOTIF
                )
                return(d)
            }) %>% do.call(rbind, .) 
        }) %>% do.call(rbind, .)
        p <- ggplot(periodicity_results, aes(x = dinuc, y = PSD, fill = dinuc)) + 
            labs(x = '', fill = 'Dinucleotide') +
            geom_col(position= "dodge") + 
            facet_wrap(~tissue, nrow = 2) + 
            theme_bw() + 
            theme(legend.position = 'bottom')
        TRUE
    }, TRUE)
})

test_that("getPeriodicity for sacCer3 random loci", {
    skip('figure_paper_20200507')
    expect_equal({
        sacCer3_random_regions <- sampleGenome(
            'sacCer3', n = 10000, w = 1000
        )
        # 
        PSDs_yeast <- getPeriodicity(
            sacCer3_random_regions, 
            motif = 'WW', 
            period = 10, 
            cores = 100, 
            skip_shuffling = FALSE
        )
        p <- plotPeriodicityResults(PSDs_yeast, xlim = 150)
        ggsave('PSDs_yeast_WW_1000longT.pdf', width = 32, height = 8, unit = 'cm')
        TRUE
    }, TRUE)
})

test_that("getPeriodicity for sacCer3 MNase", {
    skip('figure_paper_20200507')
    expect_equal({
        bam_MNase_sacCer3 <- readRDS(
            url('http://ahringerlab.com/VplotR/MNase_sacCer3_Henikoff2011.rds')
        )
        bam_MNase_sacCer3 <- unlist(GRangesList(bam_MNase_sacCer3))
        ws <- width(bam_MNase_sacCer3)
        frags <- bam_MNase_sacCer3[ws >= 147 & ws <= 152]
        # 
        PSDs_yeast <- getPeriodicity(
            frags[1:1000000], 
            genome = 'sacCer3',
            motif = 'TT', 
            period = 10, 
            cores = 20, 
            skip_shuffling = FALSE
        )
        p <- plotPeriodicityResults(PSDs_yeast, xlim = 150)
        ggsave('PSDs_yeast_MNase_TTs.pdf', width = 32, height = 8, unit = 'cm')
        TRUE
    }, TRUE)
})

test_that("getPeriodicity for ce11 proms/enhancers", {
    skip('figure_paper_20200507')
    expect_equal({
        load('~/shared/data/classification_tissue-spe-genes-REs_REs-GRanges.RData')
        ce_seq <- getSeq(
            BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
        )
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
                proms_distal <- c(
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
                    resize(140, fix = 'center')
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
                    ce_seq[proms_center],
                    ce_seq[proms_flanking],
                    ce_seq[proms_distal],
                    ce_seq[enhs_center],
                    ce_seq[enhs_flanking],
                    ce_seq[enhs_extended]
                )
            }
        ) %>% setNames(c('Ubiq.', 'Germline', 'Neurons', 'Muscle', 'Hypod.', 'Intest.'))
        #
        PSDs <- mclapply(names(list_params), mc.cores = 6, function(TISSUE) {
            list_res <- mclapply(mc.cores = 6, list_params[[TISSUE]], function(seqs) {
                res <- getPeriodicity(
                    seqs, 
                    motif = 'TT', 
                    period = 10,
                    cores = 12
                )
                res <- res$PSD$PSD[which.min(abs(res$PSD$period - 10))]
            })
            df <- data.frame(
                psd = unlist(list_res), 
                Loc = c('proms_center', 'proms_flanking', 'proms_distal', 'enhs_center', 'enhs_flanking', 'enhs_extended'), 
                Class = c('Center', 'Flanking', 'Distal', 'Center', 'Flanking', 'Distal'), 
                Type = c('Prom', 'Prom', 'Prom', 'Enh', 'Enh', 'Enh'),
                tissue = TISSUE
            )
            return(df)
        }) %>% do.call(rbind, .) %>% 
            mutate(Class = factor(Class, levels = c('Center', 'Flanking', 'Distal'))) %>% 
            mutate(Type = factor(Type, levels = c('Prom', 'Enh'))) %>% 
            mutate(tissue = factor(tissue, levels = names(list_params))) %>% 
            mutate(Loc = factor(Loc, levels = c('proms_center', 'proms_flanking', 'proms_distal', 'enhs_center', 'enhs_flanking', 'enhs_extended')))
        p <- ggplot(PSDs, aes(x = Loc, y = psd)) + 
            geom_col(
                position = "dodge", aes(col = Type, alpha = Class), size = 0.25
            ) + 
            theme_ggplot2() + 
            facet_wrap(~tissue, nrow = 1) + 
            labs(
                x = '', 
                y = 'TT PSD @ 10-bp'
            ) + 
            scale_x_discrete(labels = c('', 'Proms', '', '', 'Enhs', '')) + 
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
            ) + 
            scale_color_manual(values = c(NA, 'black'), guide = "none") + 
            scale_alpha_manual(values = c(1, 0.6, 0.3), name = '')
        ggsave(
            'PSDs_ce11_TT_proms-enhs.pdf', width = 21, height = 8, unit = 'cm'
        )
        #
        FPIs <- mclapply(names(list_params), mc.cores = 6, function(TISSUE) {
            list_res <- mclapply(mc.cores = 6, list_params[[TISSUE]], function(seqs) {
                res <- getFPI(
                    seqs, 
                    genome = 'ce11', 
                    motif = 'TT', 
                    period = 10, 
                    cores_shuffling = 4, 
                    n_shuffling = 16
                )$FPI
            })
            df <- data.frame(
                fpi = unlist(list_res), 
                Loc = c('proms_center', 'proms_flanking', 'proms_distal', 'enhs_center', 'enhs_flanking', 'enhs_extended'), 
                Class = c('Center', 'Flanking', 'Distal', 'Center', 'Flanking', 'Distal'), 
                Type = c('Prom', 'Prom', 'Prom', 'Enh', 'Enh', 'Enh'),
                tissue = TISSUE
            )
            return(df)
        }) %>% do.call(rbind, .) %>% 
            mutate(Class = factor(Class, levels = c('Center', 'Flanking', 'Distal'))) %>% 
            mutate(Type = factor(Type, levels = c('Prom', 'Enh'))) %>% 
            mutate(tissue = factor(tissue, levels = names(list_params))) %>% 
            mutate(Loc = factor(Loc, levels = c('proms_center', 'proms_flanking', 'proms_distal', 'enhs_center', 'enhs_flanking', 'enhs_extended')))
        p <- ggplot(FPIs, aes(x = Loc, y = fpi)) + 
            geom_col(position = "dodge", aes(col = Type, alpha = Class), size = 0.25) + 
            theme_ggplot2() + 
            facet_wrap(~tissue, nrow = 1) + 
            labs(
                x = '', 
                y = 'TT FPI @ 10-bp'
            ) + 
            scale_x_discrete(labels = c('', 'Proms', '', '', 'Enhs', '')) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
            scale_color_manual(values = c(NA, 'black'), guide = "none") + 
            scale_alpha_manual(values = c(1, 0.6, 0.3), name = '')
        ggsave(
            'FPIs_ce11_TT_proms-enhs.pdf', width = 21, height = 8, unit = 'cm'
        )
        #
        TRUE
    }, TRUE)
})
