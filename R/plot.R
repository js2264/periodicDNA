#' Plot the output of getPeriodicity()
#' 
#' This function plots some results from the result of getPeriodicity(). 
#' It plots the raw distogram, the distance-decay normalized 
#' distogram and the resulting PSD values. If a shuffled control has
#' been performed by getPeriodicity(), it also displays it. 
#'
#' @param results The output of getPeriodicity function.
#' @param periods Vector a numerical vector of length 2, to specify the
#' x-axis limits 
#' @param filter_periods Boolean Should the x-axis be 
#' constrained to the periods?
#' @param skip_shuffling Boolean should the shuffling sequences be done?
#' @param facet_control Boolean should the shuffling plots be faceted?
#' @param xlim Integer x axis upper limit in raw and norm. distograms
#' @param ... Additional theme arguments passed to theme_ggplot2()
#' @return list A list containing four ggplots
#' 
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' periodicity_result <- getPeriodicity(
#'     ce11_proms[1:100],
#'     genome = 'ce11',
#'     motif = 'TT', 
#'     cores = 1, 
#'     skip_shuffling = FALSE
#' )
#' head(periodicity_result$PSD)
#' plotPeriodicityResults(periodicity_result)
#' plotPeriodicityResults(periodicity_result, xlim = 150)
#' plotPeriodicityResults(
#'     periodicity_result, xlim = 150, filter_periods = FALSE
#' )
#' plotPeriodicityResults(
#'     periodicity_result, xlim = 150, facet_control = FALSE
#' )

plotPeriodicityResults <- function(
    results, 
    periods = c(2, 20), 
    filter_periods = TRUE, 
    skip_shuffling = FALSE, 
    facet_control = TRUE,
    xlim = NULL, 
    ...
) 
{
    if (is.null(xlim)) {
        xlim <- nrow(results$hist) - 6
    }
    
    if (skip_shuffling == FALSE & 'dists_shuffled' %in% names(results)) {
        d <- rbind(
            head(
                cbind(results$hist[-c(seq_len(5)),], type = 'observed'), 
                xlim
            ),
            head(
                cbind(
                    results$hist_shuffled[-c(seq_len(5)),], type = 'shuffled'
                ), 
                xlim
            )
        )
        d$type = factor(d$type)
        p0 <- ggplot2::ggplot(
            d, ggplot2::aes(x = distance, y = counts, col = type)
        ) + 
            ggplot2::geom_line(alpha = c(1, 0.4)[d$type]) +
            theme_ggplot2(...) + 
            ggplot2::labs(
                x = paste0('Distance between pairs of ', results$motif), 
                y = 'Distribution', 
                title = paste0(
                    'Distribution of\ndistances between pairs of ',
                    results$motif
                )
            ) + 
            scale_color_manual(values = c('black', 'grey30')) +
            theme(legend.position = 'none')
        #
        d <- rbind(
            head(
                cbind(
                    results$normalized_hist[-c(seq_len(5)),],
                    type = 'observed'
                ), 
                xlim
            ),
            head(
                cbind(
                    results$normalized_hist_shuffled[-c(seq_len(5)),], 
                    type = 'shuffled'
                ), 
                xlim
            )
        )
        d$type = factor(d$type)
        p1 <- ggplot2::ggplot(d, ggplot2::aes(
            x = distance, y = norm_counts, col = type
        )) + 
            ggplot2::geom_line(alpha = c(1, 0.4)[d$type]) +
            theme_ggplot2(...) + 
            ggplot2::labs(
                x = paste0('Distance between pairs of ', results$motif), 
                y = 'Normalized distribution', 
                title = paste0(
                    'Normalized distribution of\ndistances between pairs of ',
                    results$motif
                )
            ) + 
            scale_color_manual(values = c('black', 'grey30')) +
            theme(legend.position = 'none')
        if (facet_control) p1 <- p1 + ggplot2::facet_wrap(~type, nrow = 2) 
        #
        if (filter_periods == FALSE) {
            df <- data.frame(
                x = c(
                    rev(1/results$PSD$freq), rev(1/results$PSD_shuffled$freq)
                ),
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
                    results$PSD$PSD[
                        unlist(lapply(1/(periods[1]:periods[2]), 
                        function (freq) {
                            idx <- which.min(abs(freq - results$PSD$freq))
                        }
                    ))
                    ], 
                    results$PSD_shuffled$PSD[
                    unlist(lapply(
                        1/(periods[1]:periods[2]), 
                        function (freq) {
                            idx <- which.min(
                                abs(freq - results$PSD_shuffled$freq)
                            )
                        }
                    ))
                    ]
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
            ggplot2::geom_segment(
                aes(x=x, xend=x, y=0, yend=y), alpha = c(1, 0.4)[df$type]
            ) +
            ggplot2::xlim(periods) +
            theme_ggplot2(...) + 
            ggplot2::labs(
                x = paste0(results$motif, ' periods'), 
                y = 'Power Spectral Density', 
                title = paste0(
                    results$motif, ' power spectral density'
                )
            ) + 
            scale_color_manual(values = c('black', 'grey30')) +
            theme(legend.position = 'none')
        if (facet_control) p2 <- p2 + ggplot2::facet_wrap(~type) 
        #
        p <- cowplot::plot_grid(
            plotlist = list(p0, p1, p2), nrow = 1, 
            rel_widths = c(0.5, 0.5, 1)
        )
        return(p)
    }
    else {
        p0 <- ggplot2::ggplot(
            head(results$hist[-c(seq_len(5)),], xlim), 
            ggplot2::aes(x = distance, y = counts)
        ) + 
            ggplot2::geom_line() +
            theme_ggplot2(...) + 
            ggplot2::labs(
                x = paste0('Distance between pairs of ', results$motif),
                y = 'Distribution', 
                title = paste0(
                    'Distribution of\ndistances between pairs of ', 
                    results$motif
                )
            )
        #
        df <- head(results$normalized_hist[-c(seq_len(5)),], xlim)
        p1 <- ggplot2::ggplot(
            df, ggplot2::aes(x = distance, y = norm_counts)
        ) +
            ggplot2::geom_line() +
            theme_ggplot2(...) + 
            ggplot2::labs(
                x = paste0('Distance between pairs of ', results$motif), 
                y = 'Normalized distribution', 
                title = paste0(
                    'Normalized distribution of\ndistances between pairs of ',
                    results$motif
                )
            )
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
                y = results$PSD$PSD[
                unlist(lapply(1/(periods[1]:periods[2]), function (freq) {
                    idx <- which.min(abs(freq - results$PSD$freq))
                }))
                ]
            )
        }
        p2 <- ggplot2::ggplot(
            df[df$x >= periods[1] & df$x <= periods[2],], 
            ggplot2::aes(x = x, y = y)
        ) + 
            ggplot2::geom_point() +
            ggplot2::geom_segment(aes(x=x, xend=x, y=0, yend=y)) +
            ggplot2::xlim(periods) +
            theme_ggplot2(...) + 
            ggplot2::labs(
                x = paste0(results$motif, ' periods'), 
                y = 'Power Spectral Density', 
                title = paste0(
                    'Power Spectral Density of\n', 
                    results$motif, 
                    ' at different periods'
                )
            )
        #
        p <- cowplot::plot_grid(plotlist = list(p0, p1, p2), nrow = 1)
        return(p)
    }
}

#' Plot the output of getFPI()
#'
#' This function plots some results from the result of getFPI(). 
#' It plots the observed PSD values in red and all the shuffled 
#' runs in grey.
#'
#' @param fpi The output of getFPI function.
#' @param periods Vector a numerical vector of length 2, to specify the
#' x-axis limits 
#' @param s float span parameter for loess smooth
#' @return ggplot A ggplot
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' fpi <- getFPI(
#'     ce11_proms_seqs[1:10], 
#'     genome = 'ce11', 
#'     motif = 'TT', 
#'     cores_shuffling = 1
#' )
#' plotFPI(fpi)

plotFPI <- function(fpi, periods = c(2, 20), s = 0.05) {
    
    n_shuffled <- length(fpi$shuffled_spectra)
    x = c(
        1/fpi$observed_spectra$PSD$freq,
        lapply(seq_len(n_shuffled), function(k) {
            1/fpi$shuffled_spectra[[k]]$PSD$freq
        }) %>% unlist()
    )
    psds = c(
        fpi$observed_spectra$PSD$PSD,
        lapply(seq_len(n_shuffled), function(k) {
            fpi$shuffled_spectra[[k]]$PSD$PSD
        }) %>% unlist()
    )
    groups = c(
        rep(1, length(fpi$observed_spectra$PSD$PSD)),
        lapply(seq_len(n_shuffled), function(k) {
            rep(k+1, length(fpi$shuffled_spectra[[k]]$PSD$PSD))
        }) %>% unlist()
    )
    df <- data.frame(
        'x' = x, 
        'y' = psds, 
        'type' = c(
            rep('observed', length(fpi$observed_spectra$PSD$freq)), 
            rep('shuffled', length(x)-length(fpi$observed_spectra$PSD$freq))
        ), 
        'group' = groups
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
        span = s, 
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
        y = 'Power Spectral Density', 
        title = paste0('FPI @ ', fpi$period, '-bp period: ', round(fpi$FPI, 1))
    )
    return(p)
}

#' Personal ggplot2 theming function, adapted from roboto-condensed 
#' at https://github.com/hrbrmstr/hrbrthemes/
#'
#' @param base_family,base_size base font family and size
#' @param plot_title_family,plot_title_face, plot title family, face
#' @param plot_title_size,plot_title_margin, plot title size and margin
#' @param subtitle_face,subtitle_size plot subtitle family, 
#' face and size
#' @param subtitle_margin plot subtitle margin bottom (single numeric value)
#' @param strip_text_family,strip_text_face,strip_text_size facet label font 
#' family, face and size
#' @param caption_face,caption_size,caption_margin plot caption
#' family, face, size and margin
#' @param axis_title_family,axis_title_face,axis_title_size axis title font 
#' family, face and size
#' @param axis_title_just axis title font justificationk one of `[blmcrt]`
#' @param axis_text_size font size of axis text
#' @param plot_margin plot margin (specify with [ggplot2::margin])
#' @param panel_spacing panel spacing (use `unit()`)
#' @param grid_col grid color
#' @param grid panel grid (`TRUE`, `FALSE`, or a combination of 
#' `X`, `x`, `Y`, `y`)
#' @param axis_col axis color
#' @param axis add x or y axes? `TRUE`, `FALSE`, "`xy`"
#' @param ticks ticks if `TRUE` add ticks
#' @param border border if `TRUE` add border
#' @return theme A ggplot theme
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   labs(x="Fuel effiency (mpg)", y="Weight (tons)",
#'        title="Seminal ggplot2 scatterplot example") +
#'   theme_ggplot2()

theme_ggplot2 <- function(
    grid = TRUE,
    border = TRUE, 
    base_family = NULL, base_size = 8,
    plot_title_family = base_family, plot_title_size = 12,
    plot_title_face = "plain", plot_title_margin = 5,
    subtitle_size = 11,
    subtitle_face = "plain", subtitle_margin = 5,
    strip_text_family = base_family, strip_text_size = 10,
    strip_text_face = "bold",
    caption_size = 9,
    caption_face = "plain", caption_margin = 3,
    axis_text_size = base_size,
    axis_title_family = base_family,
    axis_title_size = 9,
    axis_title_face = "plain",
    axis_title_just = "rt",
    panel_spacing = grid::unit(2, "lines"),
    grid_col = "#cccccc", 
    plot_margin = margin(12, 12, 12, 12),
    axis_col = "#cccccc", 
    axis = FALSE, 
    ticks = FALSE
) 
{
    
    ret <- ggplot2::theme_minimal(
        base_family = base_family, base_size = base_size
    )
    
    ret <- ret + theme(legend.background = element_blank())
    ret <- ret + theme(legend.key = element_blank())
    ret <- ret + theme(legend.position = 'bottom')
    ret <- ret + theme(legend.title = element_text(size = 10, face="bold")) 
    ret <- ret + theme(legend.text = element_text(size = 9))
    
    ret <- ret + theme(plot.margin = plot_margin)
    
    ret <- ret + theme(panel.spacing = panel_spacing)
    
    if (inherits(grid, "character") | grid == TRUE) {
        ret <- ret + theme(
            panel.grid = element_line(color = grid_col, size = 0.2), 
            panel.grid.major = element_line(color = grid_col, size = 0.2), 
            panel.grid.minor = element_line(color = grid_col, size = 0.15)
        )
        
        if (inherits(grid, "character")) {
            if (regexpr("X", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.major.x = element_blank())
            if (regexpr("Y", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.major.y = element_blank())
            if (regexpr("x", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.minor.x = element_blank())
            if (regexpr("y", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.minor.y = element_blank())
        }
    } else {
        ret <- ret + theme(panel.grid = element_blank())
        ret <- ret + theme(panel.grid.major  = element_blank())
        ret <- ret + theme(panel.grid.major.x  = element_blank())
        ret <- ret + theme(panel.grid.major.y  = element_blank())
        ret <- ret + theme(panel.grid.minor  = element_blank())
        ret <- ret + theme(panel.grid.minor.x  = element_blank())
        ret <- ret + theme(panel.grid.minor.y  = element_blank())
    }
    
    if (border == TRUE) {
        ret <- ret + theme(
            panel.border = element_rect(
                colour = "black", fill = NA, size = 0.5
            )
        )
    }
    
    if (inherits(axis, "character") | axis == TRUE) {
        ret <- ret + theme(
            axis.line = element_line(color = axis_col, size = 0.15)
        )
        if (inherits(axis, "character")) {
            axis <- tolower(axis)
            if (regexpr("x", axis)[1] < 0) {
                ret <- ret + theme(axis.line.x = element_blank())
            } else {
                ret <- ret + theme(
                    axis.line.x = element_line(color = axis_col, size = 0.15)
                )
            }
            if (regexpr("y", axis)[1] < 0) {
                ret <- ret + theme(axis.line.y = element_blank())
            } else {
                ret <- ret + theme(
                    axis.line.y = element_line(color = axis_col, size = 0.15)
                )
            }
        } else {
            ret <- ret + theme(
                axis.line.x = element_line(color = axis_col, size = 0.15), 
                axis.line.y = element_line(color = axis_col, size = 0.15)
            )
        }
    } else {
        ret <- ret + theme(axis.line = element_blank())
    }
    
    if (!ticks) {
        ret <- ret + theme(axis.ticks = element_blank())
        ret <- ret + theme(axis.ticks.x = element_blank())
        ret <- ret + theme(axis.ticks.y = element_blank())
    } else {
        ret <- ret + theme(axis.ticks = element_line(size = 0.15))
        ret <- ret + theme(axis.ticks.x = element_line(size = 0.15))
        ret <- ret + theme(axis.ticks.y = element_line(size = 0.15))
        ret <- ret + theme(axis.ticks.length = grid::unit(5, "pt"))
    }
    
    xj <- switch(
        tolower(substr(axis_title_just, 1, 1)), 
        b = 0, l = 0, m = 0.5, c = 0.5, r = 1, t = 1
    )
    yj <- switch(
        tolower(substr(axis_title_just, 2, 2)), 
        b = 0, l = 0, m = 0.5, c = 0.5, r = 1, t = 1
    )
    
    ret <- ret + theme(axis.text = element_text(
        size = axis_text_size, margin = margin(t = 0, r = 0)
    ))
    ret <- ret + theme(axis.text.x = element_text(
        size = axis_text_size, margin = margin(t = 0)
    ))
    ret <- ret + theme(axis.text.y = element_text(
        size = axis_text_size, margin = margin(r = 0)
    ))
    
    ret <- ret + theme(axis.title = element_text(
        size = axis_title_size, family = axis_title_family)
    )
    ret <- ret + theme(axis.title.x = element_text(
        # hjust = xj, 
        size = axis_title_size,
        family = axis_title_family, face = axis_title_face
    ))
    ret <- ret + theme(axis.title.y = element_text(
        # hjust = yj, 
        size = axis_title_size,
        family = axis_title_family, face = axis_title_face
    ))
    ret <- ret + theme(axis.title.y.right = element_text(
        # hjust = yj, 
        size = axis_title_size, angle = 90,
        family = axis_title_family, face = axis_title_face
    ))
    
    ret <- ret + theme(strip.text = element_text(
        hjust = 0, size = strip_text_size,
        face = strip_text_face, family = strip_text_family
    ))
    
    ret <- ret + theme(plot.title = element_text(
        hjust = 0, size = plot_title_size,
        margin = margin(b = plot_title_margin),
        family = plot_title_family, face = plot_title_face
    ))
    ret <- ret + theme(plot.subtitle = element_text(
        hjust = 0, size = subtitle_size,
        margin = margin(b = subtitle_margin),
        face = subtitle_face
    ))
    ret <- ret + theme(plot.caption = element_text(
        hjust = 1, size = caption_size,
        margin = margin(t = caption_margin),
        face = caption_face
    ))
    
    ret
    
}