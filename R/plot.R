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

plotPeriodicityResults <- function(results, periods = c(2, 20), filter_periods = FALSE, skip_shuffling = FALSE) {
    if (skip_shuffling == FALSE & 'dists_shuffled' %in% names(results)) {
        d <- rbind(
            cbind(results$hist[-c(1:5),], type = 'observed'),
            cbind(results$hist_shuffled[-c(1:5),], type = 'shuffled')
        )
        d$type = factor(d$type)
        p0 <- ggplot2::ggplot(d, ggplot2::aes(x = distance, y = counts, col = type)) + 
            ggplot2::geom_line(alpha = c(1, 0.4)[d$type]) +
            ggplot2::theme_classic() + 
            ggplot2::labs(x = paste0('Distance between pairs of ', results$motif), y = 'Distribution', title = paste0('Distribution of\ndistances between pairs of ', results$motif))
        #
        d <- rbind(
            cbind(results$normalized_hist %>% '['(-c(1:5, (nrow(.)-1):nrow(.)),), type = 'observed'),
            cbind(results$normalized_hist_shuffled %>% '['(-c(1:5, (nrow(.)-1):nrow(.)),), type = 'shuffled')
        )
        d$type = factor(d$type)
        p1 <- ggplot2::ggplot(d) + 
            ggplot2::aes(x = distance, y = norm_counts, col = type) + 
            ggplot2::geom_line(alpha = c(1, 0.4)[d$type]) +
            ggplot2::theme_classic() + 
            ggplot2::labs(x = paste0('Distance between pairs of ', results$motif), y = 'Normalized distribution', title = paste0('Normalized distribution of\ndistances between pairs of ', results$motif))
        #
        if (filter_periods == FALSE) {
            df <- data.frame(
                x = c(rev(1/results$PSD$freq), rev(1/results$PSD_shuffled$freq)),
                y = c(rev(results$PSD$PSD), rev(results$PSD_shuffled$PSD)), 
                type = c(
                    rep('observed', length(results$PSD$PSD)),
                    rep('shuffled', length(results$PSD_shuffled$PSD))
                )
            )
            df <- df[df$x >= periods[1] & df$x <= periods[2],]
        }
        else {
            df <- data.frame(
                x = c(periods[1]:periods[2], periods[1]:periods[2]),
                y = c(
                    results$PSD$PSD[unlist(lapply(1/(periods[1]:periods[2]), function (freq) {
                        idx <- which.min(abs(freq - results$PSD$freq))
                    }))], 
                    results$PSD_shuffled$PSD[unlist(lapply(1/(periods[1]:periods[2]), function (freq) {
                        idx <- which.min(abs(freq - results$PSD_shuffled$freq))
                    }))]
                ), 
                type = c(
                    rep('observed', length(periods[1]:periods[2])),
                    rep('shuffled', length(periods[1]:periods[2]))
                )
            )
            df <- df[df$x >= periods[1] & df$x <= periods[2],]
        }
        p2 <- ggplot2::ggplot(df) + 
            ggplot2::aes(x = x, y = y, col = type, fill = type) + 
            ggplot2::geom_point(alpha = c(1, 0.4)[df$type]) +
            ggplot2::geom_segment(aes(x=x, xend=x, y=0, yend=y), alpha = c(1, 0.4)[df$type]) +
            ggplot2::xlim(periods) +
            ggplot2::theme_classic() + 
            ggplot2::labs(
                x = paste0(results$motif, ' frequency'), 
                y = 'Power Spectrum Density', 
                title = paste0('Power Spectrum Density of\n', results$motif, ' frequencies')
            )
        #
        return(list(p0, p1, p2))
    }
    else {
        p0 <- ggplot2::ggplot(results$hist[-c(1:5),], ggplot2::aes(x = distance, y = counts)) + 
            ggplot2::geom_line() +
            ggplot2::theme_classic() + 
            ggplot2::labs(x = paste0('Distance between pairs of ', results$motif), y = 'Distribution', title = paste0('Distribution of\ndistances between pairs of ', results$motif))
        #
        df <- results$normalized_hist %>% '['(-c(1:5, (nrow(.)-15):nrow(.)),)
        p1 <- ggplot2::ggplot(df, ggplot2::aes(x = distance, y = norm_counts)) + 
            ggplot2::geom_line() +
            ggplot2::theme_classic() + 
            ggplot2::labs(x = paste0('Distance between pairs of ', results$motif), y = 'Normalized distribution', title = paste0('Normalized distribution of\ndistances between pairs of ', results$motif))
        #
        if (filter_periods == FALSE) {
            df <- data.frame(
                x = rev(1/results$PSD$freq),
                y = rev(results$PSD$PSD)
            )
        }
        else {
            df <- data.frame(
                x = periods[1]:periods[2],
                y = results$PSD$PSD[unlist(lapply(1/(periods[1]:periods[2]), function (freq) {
                    idx <- which.min(abs(freq - results$PSD$freq))
                }))]
            )
        }
        p2 <- ggplot2::ggplot(df[df$x >= periods[1] & df$x <= periods[2],], ggplot2::aes(x = x, y = y)) + 
            ggplot2::geom_point() +
            ggplot2::geom_segment(aes(x=x, xend=x, y=0, yend=y)) +
            ggplot2::xlim(periods) +
            ggplot2::theme_classic() + 
            ggplot2::labs(x = paste0(results$motif, ' frequency'), y = 'Power Spectrum Density', title = paste0('Power Spectrum Density of\n', results$motif, ' frequencies'))
        #
        return(list(p0, p1, p2))
    }
}

#' Plotting function
#'
#' @param results The output of FPI function.
#' @param periods Vector a numerical vector of length 2, to specify the
#' x-axis limits 
#' 
#' @return ggplot A ggplot
#' 
#' @import ggplot2
#' @export

plotFPI <- function(fpi, periods = c(2, 20)) {
    n_shuffled <- length(fpi$shuffled_spectra)
    x = c(
        1/fpi$observed_spectra$PSD$freq,
        lapply(1:n_shuffled, function(k) {
            1/fpi$shuffled_spectra[[k]]$PSD$freq
        }) %>% unlist()
    )
    psds = c(
        fpi$observed_spectra$PSD$PSD,
        lapply(1:n_shuffled, function(k) {
            fpi$shuffled_spectra[[k]]$PSD$PSD
        }) %>% unlist()
    )
    groups = c(
        rep(1, length(fpi$observed_spectra$PSD$PSD)),
        lapply(1:n_shuffled, function(k) {
            rep(k+1, length(fpi$shuffled_spectra[[k]]$PSD$PSD))
        }) %>% unlist()
    )
    df <- data.frame(
        x = x, 
        y = psds, 
        type = c(
            rep('observed', length(fpi$observed_spectra$PSD$freq)), 
            rep('shuffled', length(x)-length(fpi$observed_spectra$PSD$freq))
        ), 
        group = groups
    )
    df <- df[x >= periods[1] & x <= periods[2],]
    p <- ggplot2::ggplot(df) + 
        ggplot2::aes(x = x, y = y) + 
        ggplot2::geom_point(
            alpha = c(1, 0.3)[df$type], 
            col = c('red', 'grey50')[df$type]
        ) +
        ggplot2::geom_line(
            data = df[df$type != 'observed',],
            aes(x = x, y = y, group = group), 
            stat = "smooth", 
            method = "loess", 
            span = 0.05, 
            col = 'grey50',
            alpha = 0.5, 
        ) +
        ggplot2::geom_smooth(
            data = df[df$type == 'observed',],
            aes(x = x, y = y, group = group), 
            col = 'red',
            alpha = 0, 
            method = "loess",
            span = 0.05
        ) +
        ggplot2::xlim(periods) +
        ggplot2::theme_classic() + 
        ggplot2::labs(
            x = paste0(fpi$motif, ' periods'), 
            y = 'Power Spectrum Density', 
            title = paste0('FPI @ ', fpi$period, 'bp: ', round(fpi$FPI, 1))
        )
    return(p)
}
