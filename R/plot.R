plotPeriodicityResults <- function(results, HIST_YLIM = c(-1.5, 1.5)) {
    require(ggplot2)
    p1 <- ggplot(results$normalized_hist, aes(x = distance, y = norm_counts)) + 
        geom_line() +
        ylim(HIST_YLIM) + 
        theme_classic() + 
        labs(x = paste0('Distance between pairs of ', results$motif), y = 'Normalized frequency', title = paste0('Normalized frequency of\ndistances between pairs of ', results$motif))
    p2 <- ggplot(results$PSD, aes(x = freq, y = PSD)) + 
        geom_point() +
        geom_line() +
        theme_classic() + 
        labs(x = paste0(results$motif, ' frequency'), y = 'Power Spectrum Density', title = paste0('Power Spectrum Density of\n', results$motif, ' frequencies'))
    p3 <- ggplot(results$signal_to_noise_ratio, aes(x = period, y = SNR)) + 
        geom_point() +
        geom_line() +
        theme_classic() + 
        labs(x = paste0(results$motif, ' period'), y = 'Signal to noise ratio', title = paste0('Strenght of\n', results$motif, ' periodicity'))
    return(list(p1, p2, p3))
}
