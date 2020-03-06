#' Plotting function
#'
#' @param results The output of getPeriodicity function.
#' @param norm_hist_ylim Vector a numerical vector of length 2, to specify the y-axis 
#' limits of the first plot (normalized distribution of pairwise distances).
#' 
#' @return list A list containing four ggplots
#' 
#' @import ggplot2
#' @export

plotPeriodicityResults <- function(results, norm_hist_ylim = c(-1.5, 1.5)) {
    p0 <- ggplot2::ggplot(results$hist[-c(1:5),], ggplot2::aes(x = distance, y = counts)) + 
        ggplot2::geom_line() +
        ggplot2::theme_classic() + 
        ggplot2::labs(x = paste0('Distance between pairs of ', results$motif), y = 'Distribution', title = paste0('Distribution of\ndistances between pairs of ', results$motif))
    p1 <- ggplot2::ggplot(results$normalized_hist, ggplot2::aes(x = distance, y = norm_counts)) + 
        ggplot2::geom_line() +
        ggplot2::ylim(norm_hist_ylim) + 
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
    return(list(p0, p1, p2, p3))
}
