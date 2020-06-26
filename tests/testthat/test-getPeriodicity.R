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
        data(ce11_proms)
        periodicity_result <- getPeriodicity(
            ce11_proms[seq_len(10)],
            range_spectrum = 1:100,
            genome = 'ce11',
            motif = 'TT',
            n_shuffling = 3, 
            doZscore = TRUE
        )
        list_plots <- plotPeriodicityResults(periodicity_result)
        methods::is(list_plots, "gg")
    }, TRUE)
})

test_that("getPeriodicity + FPI Yeast", {
    skip('skip2')
    expect_equal({
        set.seed(52) 
        ####
        ####
        #### getPeriodicity with shuffling
        set.seed(52) 
        yeast_seqs <- sampleGenome('sacCer3', n = 1000, w = 800)
        yeast_PSDs <- getPeriodicity(
            yeast_seqs, 
            motif = 'WW', 
            n_shuffling = 30, 
            cores_shuffling = 30
        )
        plots <- plotPeriodicityResults(yeast_PSDs, xlim = 150)
        ggsave('tmp.pdf', width = 12, height = 4)
        ####
        ####
        #### Order 1
        set.seed(52) 
        yeast_seqs <- sampleGenome('sacCer3', n = 1000, w = 800)
        FPI <- getFPI(
            yeast_seqs, 
            motif = 'WW', 
            n_shuffling = 30, 
            cores_shuffling = 30
        )
        p <- plotFPI(FPI)
        ggsave('tmp2.pdf')
        ####
        ####
        #### Order 2
        set.seed(52) 
        yeast_seqs <- sampleGenome('sacCer3', n = 1000, w = 800)
        FPI <- getFPI(
            yeast_seqs, 
            motif = 'WW', 
            n_shuffling = 10, 
            cores_shuffling = 1, 
            order = 2
        )
        p <- plotFPI(FPI)
        ggsave('tmp3.pdf')
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
        set.seed(51)
        sacCer3_random_regions <- sampleGenome(
            'sacCer3', n = 1000, w = 800
        )
        # 
        PSDs_yeast <- getPeriodicity(
            sacCer3_random_regions, 
            motif = 'WW', 
            n_shuffling = 10, 
            cores_shuffling = 10
        )
        #
        p <- plotPeriodicityResults(PSDs_yeast, xlim = 150)
        ggsave('PSDs_yeast_WW_800long.pdf', width = 28, height = 9, unit = 'cm')
        #
        # Value for varying sequences lengths
        g <- sampleGenome('sacCer3', n = 1000, w = 1500)
        PSDs_yeast_1500 <- getPeriodicity(
            g, 
            motif = 'WW', 
            period = 10, 
            BPPARAM = setUpBPPARAM(12)
        )
        p <- plotPeriodicityResults(PSDs_yeast_1500, xlim = 150)
        ggsave('PSDs_yeast_WW_1500long.pdf', width = 32, height = 8, unit = 'cm')
        g_500 <- c(
            resize(g, fix = 'start', width = 500), 
            resize(g, fix = 'center', width = 500), 
            resize(g, fix = 'end', width = 500)
        )
        PSDs_yeast_500 <- getPeriodicity(
            g_500, 
            motif = 'WW', 
            period = 10, 
            BPPARAM = setUpBPPARAM(12)
        )
        p <- plotPeriodicityResults(PSDs_yeast_500, xlim = 150)
        ggsave('PSDs_yeast_WW_500long.pdf', width = 32, height = 8, unit = 'cm')
        #
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
            BPPARAM = setUpBPPARAM(20), 
            range_spectrum = 1:100
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
                proms_full <- resize(proms_center, 700, fix = 'center')
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
                enhs_full <- resize(enhs_center, 700, fix = 'center')
                l <- list(
                    ce_seq[proms_center],
                    ce_seq[proms_flanking],
                    ce_seq[proms_distal],
                    ce_seq[proms_full],
                    ce_seq[enhs_center],
                    ce_seq[enhs_flanking],
                    ce_seq[enhs_extended], 
                    ce_seq[enhs_full]
                )
            }
        ) %>% setNames(c('Ubiq.', 'Germline', 'Neurons', 'Muscle', 'Hypod.', 'Intest.'))
        library(BiocParallel)
        #
        PSDs <- bplapply(BPPARAM = setUpBPPARAM(1), names(list_params), function(TISSUE) {
            list_res <- bplapply(BPPARAM = setUpBPPARAM(8), list_params[[TISSUE]], function(seqs) {
                res <- getPeriodicity(
                    seqs, 
                    motif = 'TT', 
                    period = 10,
                    range_spectrum = 1:100, 
                    n_shuffling = 10, 
                    cores_shuffling = 10
                )
                list(
                    psd = res$PSD$PSD[which.min(abs(res$PSD$period - 10))], 
                    pvalue = getSignificantPeriods(res$FPI) %>% filter(Period == 10) %>% select(pval) %>% '[['(1),
                    fpi = res$FPI$FPI
                )
            })
            df <- data.frame(
                PSD = unlist(lapply(list_res, '[[', 'psd')), 
                pval = unlist(lapply(list_res, '[[', 'pvalue')), 
                FPI = unlist(lapply(list_res, '[[', 'fpi')), 
                Loc = c('proms_center', 'proms_flanking', 'proms_distal', 'proms_full', 'enhs_center', 'enhs_flanking', 'enhs_extended', 'enhs_full'), 
                Class = c('Core', 'Flanking', 'Distal', 'Full', 'Core', 'Flanking', 'Distal', 'Full'), 
                Type = c('Prom', 'Prom', 'Prom', 'Prom', 'Enh', 'Enh', 'Enh', 'Prom'),
                tissue = TISSUE
            )
            return(df)
        })
        PSDs_2 <- PSDs %>% 
            do.call(rbind, .) %>% 
            mutate(Class = factor(Class, levels = c('Core', 'Flanking', 'Distal', 'Full'))) %>% 
            dplyr::filter(Class %in% c('Core', 'Flanking', 'Distal')) %>% 
            mutate(Type = factor(Type, levels = c('Prom', 'Enh'))) %>% 
            mutate(tissue = factor(tissue, levels = names(list_params))) %>% 
            mutate(Loc = factor(Loc, levels = c('proms_center', 'proms_flanking', 'proms_distal', 'proms_full', 'enhs_center', 'enhs_flanking', 'enhs_extended', 'enhs_full')))
        #
        p <- ggplot(PSDs_2, aes(x = Loc, y = PSD)) + 
            geom_col(
                position = "dodge", aes(col = Type, alpha = Class), size = 0.25
            ) + 
            theme_ggplot2() + 
            facet_wrap(~tissue, nrow = 1) + 
            labs(
                x = '', 
                y = 'TT PSD @ 10-bp'
            ) + 
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
            ) + 
            scale_color_manual(values = c(NA, 'black'), guide = "none") + 
            # scale_x_discrete(labels = c('', 'Proms', '', '', '', 'Enhs', '', '')) + 
            # scale_alpha_manual(values = c(1, 0.6, 0.4, 0.2), name = '')
            scale_x_discrete(labels = c('', 'Proms', '', '', 'Enhs', '')) + 
            scale_alpha_manual(values = c(1, 0.6, 0.3), name = '')
        ggsave(
            'PSDs_ce11_TT_proms-enhs.pdf', width = 21, height = 12, unit = 'cm'
        )
        #
        p <- ggplot(PSDs_2, aes(x = Loc, y = FPI)) + 
            geom_col(
                position = "dodge", aes(col = Type, fill = Class), size = 0.25
            ) + 
            theme_ggplot2(panel_spacing = unit(0.5, units = 'cm')) + 
            facet_wrap(~tissue, nrow = 1) + 
            labs(
                x = '', 
                y = 'TT FPI @ 10-bp'
            ) + 
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
            ) + 
            scale_color_manual(values = c(NA, 'black'), guide = "none") + 
            scale_x_discrete(labels = c('', 'Proms', '', '', 'Enhs', '')) + 
            scale_fill_manual(values = c('grey10', 'grey50', 'grey80'), name = '')
        ggsave(
            'FPIs_ce11_TT_proms-enhs.pdf', width = 21, height = 8, unit = 'cm'
        )
        #
        TRUE
    }, TRUE)
})
