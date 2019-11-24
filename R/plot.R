#' Plotting function
#'
#' @param results The output of getPeriodicity function.
#' @param HIST_YLIM Vector a numerical vector of length 2, to specify the y-axis 
#' limits of the first plot (normalized distribution of pairwise distances).
#' 
#' @return list A list containing three ggplots
#' 
#' @import ggplot2
#' @export

plotPeriodicityResults <- function(results, HIST_YLIM = c(-1.5, 1.5)) {
    p1 <- ggplot2::ggplot(results$normalized_hist, ggplot2::aes(x = distance, y = norm_counts)) + 
        ggplot2::geom_line() +
        ggplot2::ylim(HIST_YLIM) + 
        ggplot2::theme_classic() + 
        ggplot2::labs(x = paste0('Distance between pairs of ', results$motif), y = 'Normalized distribution', title = paste0('Normalized distribution of\ndistances between pairs of ', results$motif))
    p2 <- ggplot2::ggplot(tail(results$PSD, -2), ggplot2::aes(x = freq, y = PSD)) + 
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_classic() + 
        ggplot2::labs(x = paste0(results$motif, ' frequency'), y = 'Power Spectrum Density', title = paste0('Power Spectrum Density of\n', results$motif, ' frequencies'))
    p3 <- ggplot2::ggplot(results$signal_to_noise_ratio, ggplot2::aes(x = period, y = SNR)) + 
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_classic() + 
        ggplot2::labs(x = paste0(results$motif, ' period'), y = 'Signal to noise ratio', title = paste0('Strenght of\n', results$motif, ' periodicity'))
    return(list(p1, p2, p3))
}
