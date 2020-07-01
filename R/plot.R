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
#' @param facet_control Boolean should the shuffling plots be faceted?
#' @param xlim Integer x axis upper limit in raw and norm. distograms
#' @param pval_threshold Float, significance threshold
#' @param ... Additional theme arguments passed to theme_ggplot2()
#' @return list A list containing four ggplots
#' 
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @export
#' 
#' @examples
#' data(ce11_TSSs)
#' periodicity_result <- getPeriodicity(
#'     ce11_TSSs[['Ubiq.']][1:100],
#'     genome = 'ce11',
#'     motif = 'TT', 
#'     BPPARAM = setUpBPPARAM(1)
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
    facet_control = TRUE,
    xlim = NULL, 
    pval_threshold = 0.01,
    ...
) 
{
    if (is.null(xlim)) {
        xlim <- nrow(results$hist) - 6
    }
    p0 <- ggplot2::ggplot(
        head(results$hist[-c(seq_len(5)),], xlim), 
        ggplot2::aes(x = distance, y = counts)
    ) + 
        ggplot2::geom_line() +
        theme_ggplot2() + 
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
        theme_ggplot2() + 
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
    if (!("FPI" %in% names(results))) {
        p2 <- ggplot2::ggplot(
            df[df$x >= periods[1] & df$x <= periods[2],], 
            ggplot2::aes(x = x, y = y)
        ) + 
            ggplot2::geom_point() +
            ggplot2::geom_segment(aes(x=x, xend=x, y=0, yend=y)) +
            ggplot2::xlim(periods) +
            theme_ggplot2() + 
            ggplot2::labs(
                x = paste0(results$motif, ' periods'), 
                y = 'Power Spectral Density', 
                title = paste0(
                    'Power Spectral Density of\n', 
                    results$motif, 
                    ' at different periods'
                )
            )
    }
    else {
        p2 <- plotFPI(results$PSD_withShuffling, 
            periods = periods, 
            pval_threshold = pval_threshold
        )
        p2 <- p2 + ggplot2::labs(
            x = paste0(results$motif, ' periods'), 
            y = 'Power Spectral Density', 
            title = paste0(
                'Observed PSD (red) vs. PSD in\nshuffled sequences (grey, n=',
                length(results$PSD_withShuffling$shuffled_PSD), ')'
            )
        )
    }
    #
    p <- cowplot::plot_grid(plotlist = list(p0, p1, p2), nrow = 1)
    return(p)
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
#' @param pval_threshold Float, significance threshold
#' @return ggplot A ggplot
#' 
#' @import ggplot2
#' @importFrom stats t.test
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' fpi <- getFPI(
#'     ce11_proms_seqs[1:10], 
#'     motif = 'TT', 
#'     cores_computing = 1
#' )
#' plotFPI(fpi)

plotFPI <- function(
    fpi, 
    periods = c(2, 20), 
    pval_threshold = 0.01
) {
    periods_bounds <- periods 
    n_shuffled <- length(fpi$shuffled_spectra)
    periods <- c(
        1/fpi$observed_spectra$PSD$freq,
        lapply(seq_len(n_shuffled), function(k) {
            1/fpi$shuffled_spectra[[k]]$PSD$freq
        }) %>% unlist()
    )
    psds <- c(
        fpi$observed_spectra$PSD$PSD,
        lapply(seq_len(n_shuffled), function(k) {
            fpi$shuffled_spectra[[k]]$PSD$PSD
        }) %>% unlist()
    )
    groups <- c(
        rep(1, length(fpi$observed_spectra$PSD$PSD)),
        lapply(seq_len(n_shuffled), function(k) {
            rep(k+1, length(fpi$shuffled_spectra[[k]]$PSD$PSD))
        }) %>% unlist()
    )
    # Calculate CI
    mat <- matrix(psds[groups != 1], ncol = n_shuffled, byrow = FALSE)
    conint <- apply(mat, 1, function (n) {
        Q <- stats::qnorm(0.975, sum(!is.na(n))-1)
        S <- sd(n, na.rm=TRUE)
        sq <- sqrt(sum(!is.na(n)))
        Q * S / sq
    })
    ribbon_coords <- data.frame(
        periods = unique(periods), 
        means = rowMeans(mat), 
        meansUp = rowMeans(mat) + conint, 
        meansDown = rowMeans(mat) - conint, 
        type = 'shuffled'
    )
    ribbon_coords$meansDown[ribbon_coords$meansDown < 0] <- 0
    # Making data frame
    df <- data.frame(
        'periods' = periods, 
        'y' = psds, 
        'type' = c(
            rep('observed', length(fpi$observed_spectra$PSD$freq)), 
            rep('shuffled', length(periods)-length(fpi$observed_spectra$PSD$freq))
        ), 
        'group' = groups
    )
    df <- df[periods >= periods_bounds[1] & periods <= periods_bounds[2],]
    ribbon_coords <- ribbon_coords[
        ribbon_coords$periods >= periods_bounds[1] & 
        ribbon_coords$periods <= periods_bounds[2],
    ]
    # Adding significance
    pvals <- fpi$significantPeriods
    if (!("pval" %in% colnames(pvals))) {
        pvals$pval <- 1
    }
    pvals <- pvals[
        pvals$Period >= periods_bounds[1] & pvals$Period <= periods_bounds[2],
    ]
    df$pval <- 1
    df$pval[seq_len(nrow(pvals))] <- pvals$pval
    df$isSign <- factor(FALSE, levels = c(FALSE, TRUE))
    df$isSign[df$pval <= pval_threshold] <- TRUE
    df$isSign[df$type != 'observed'] <- NA
    # Plotting
    p <- ggplot2::ggplot(df) + 
        ggplot2::aes(x = periods, y = y) + 
        ggplot2::geom_line(
            data = df[df$type != 'observed',],
            aes(x = periods, y = y, group = group), 
            col = 'grey50',
            alpha = 0.5, 
            size = 0.5
        ) +
        ggplot2::geom_ribbon(
            data = ribbon_coords,
            ggplot2::aes(x = periods, y = means, ymin = meansDown, ymax = meansUp), 
            alpha = 0.2, col = NA
        ) + 
        ggplot2::geom_line(
            data = df[df$type == 'observed',],
            aes(x = periods, y = y, group = group), 
            col = 'red',
            alpha = 1, 
            size = 1.2
        ) +
        ggplot2::geom_point(
            data = df[df$type == 'observed' & df$isSign == FALSE,], 
            ggplot2::aes(
                x = periods, 
                y = y,
            ), 
            alpha = 1,
            col = 'black', 
            fill = 'white', 
            shape = 21, 
            size = 1.5
        ) +
        ggplot2::geom_point(
            data = df[df$type == 'observed' & df$isSign == TRUE,], 
            ggplot2::aes(
                x = periods, 
                y = y,
            ), 
            alpha = 1,
            col = 'black', 
            fill = 'red', 
            shape = 21, 
            size = 2
        ) +
        ggplot2::guides(alpha = FALSE, col = FALSE) + 
        ggplot2::xlim(periods_bounds) +
        theme_ggplot2() + 
        ggplot2::labs(
            x = paste0(fpi$motif, ' periods'), 
            y = 'Power Spectral Density', 
            title = paste0(
                'FPI @ ', fpi$period, '-bp period: ', round(fpi$FPI, 1), 
                '\n(p = ', df$pval[which.min(abs(df$periods - fpi$period))], ')'
            )
        ) + 
        ggplot2::theme(legend.position = "null")
    return(p)
}

#' Personal ggplot2 theming function, adapted from roboto-condensed 
#' at https://github.com/hrbrmstr/hrbrthemes/
#'
#' @param base_size base font size
#' @param plot_title_face, plot title face
#' @param plot_title_size,plot_title_margin, plot title size and margin
#' @param subtitle_face,subtitle_size plot subtitle face and size
#' @param subtitle_margin plot subtitle margin bottom (single numeric value)
#' @param strip_text_face,strip_text_size facet label font face and size
#' @param caption_face,caption_size,caption_margin plot caption face, size 
#' and margin
#' @param axis_title_face,axis_title_size axis title font face and size
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
    base_size = 8,
    plot_title_size = 12,
    plot_title_face = "plain", plot_title_margin = 5,
    subtitle_size = 11,
    subtitle_face = "plain", subtitle_margin = 5,
    strip_text_size = 10,
    strip_text_face = "bold",
    caption_size = 9,
    caption_face = "plain", caption_margin = 3,
    axis_text_size = base_size,
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
        base_size = base_size
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
        size = axis_title_size)
    )
    ret <- ret + theme(axis.title.x = element_text(
        # hjust = xj, 
        size = axis_title_size,
        face = axis_title_face
    ))
    ret <- ret + theme(axis.title.y = element_text(
        # hjust = yj, 
        size = axis_title_size,
        face = axis_title_face
    ))
    ret <- ret + theme(axis.title.y.right = element_text(
        # hjust = yj, 
        size = axis_title_size, angle = 90,
        face = axis_title_face
    ))
    
    ret <- ret + theme(strip.text = element_text(
        hjust = 0, size = strip_text_size,
        face = strip_text_face
    ))
    
    ret <- ret + theme(plot.title = element_text(
        hjust = 0, size = plot_title_size,
        margin = margin(b = plot_title_margin),
        face = plot_title_face
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