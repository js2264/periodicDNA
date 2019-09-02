extendBothSides <- function(granges, upstream, downstream) {
    granges.bak <- granges
    granges <- GenomicRanges::resize(granges, width = GenomicRanges::width(granges.bak) + downstream, fix = 'start')
    granges <- GenomicRanges::resize(granges, width = upstream + GenomicRanges::width(granges.bak) + downstream, fix = 'end')
    return(granges)
}
AlignToTSS <- function(granges, upstream, downstream) {
    if (any(GenomicRanges::strand(granges) == '*')) {
        granges <- deconvolveBidirectionalPromoters(granges)
    }
    if (is.null(granges$TSS.fwd)) {
        granges$TSS.fwd <- IRanges::start(granges)
        granges$TSS.rev <- IRanges::end(granges)
    }
    GenomicRanges::ranges(granges) <- IRanges::IRanges(
        start = ifelse(as.vector(GenomicRanges::strand(granges)) == '+', (granges$TSS.fwd - upstream), (granges$TSS.rev - downstream + 1)),
        width = downstream + upstream,
        names = names(IRanges::ranges(granges))
    )
    return(granges)
}
deconvolveBidirectionalPromoters <- function(granges) {
    unid <- granges[GenomicRanges::strand(granges) == '+' | GenomicRanges::strand(granges) == '-']
    bid <- granges[GenomicRanges::strand(granges) == '*']
    bid.fwd <- bid
    GenomicRanges::strand(bid.fwd) <- '+'
    bid.rev <- bid
    GenomicRanges::strand(bid.rev) <- '-'
    granges_shifted <- sort(c(unid, bid.fwd, bid.rev), ignore.strand = T)
    return(granges_shifted)
}
withSeq <- function(granges, genome) {
    granges$seq <- genome[granges]
    return(granges)
}

# ----------------

getCovMatrix <- function(granges, bw.file, norm = 'none', center = FALSE, flank = NULL, stranded = TRUE, split.bid.proms = TRUE, verbose = TRUE) {
    scores <- bw.file
    # Deconvolute and resize GRanges 
    if (any(GenomicRanges::strand(granges) == '*')) {
        if (stranded == T & split.bid.proms == T) {
            if (verbose) message('.. Bi-directional segments found and split...')
            granges <- deconvolveBidirectionalPromoters(granges)
        } else {
            if (verbose) message('.. Bi-directional promoters found but not split...')
        }
    } else {
        if (verbose) message('.. No bi-directional segments found. Continuing...')
    }
    if (center == T) {
        if (verbose) message('.. Centering GRanges...')
        granges <- IRanges::resize(granges, width = 1, fix = 'center')
    }
    if (length(flank) > 0) {
        if (verbose) message('.. Resizing GRanges...')
        granges <- IRanges::resize(granges, fix = 'end', width = flank[1])
        granges <- IRanges::resize(granges, fix = 'start', width = flank[1] + flank[2])
    }
    # Subset scores over GRanges
    scores.subset <- scores[granges]
    # Turn it into a rectangular matrix and flip reverse strand scores
    if (verbose) message('.. Converting scores into matrix...')
    scores.subset <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(scores.subset)), nrow = length(granges), byrow = T)))
    scores.subset.flipped <- t(apply(scores.subset, 1, rev))
    if (stranded) {
        scores.subset <- matrix(sapply(1:nrow(scores.subset), function(K) {if((as.vector(GenomicRanges::strand(granges)) == '-')[K]) {scores.subset.flipped[K,]} else {scores.subset[K,]}}), nrow = length(granges), byrow = T)
    }
    
    # Normalize matrix
    if (norm == 'zscore') {
        scores.subset <- apply(scores.subset, 1, scale) %>% t() %>% na.replace(0)
    } else if (norm == 'log2') {
        scores.subset <- log2(scores.subset)
    }
    # Return matrix
    return(scores.subset)
    
}

getMotifMatrix <- function(granges, motif, seqs = Biostrings::readDNAStringSet("~/shared/sequences/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"), center = FALSE, flank = NULL, stranded = TRUE, split.bid.proms = TRUE, verbose = FALSE) {
    # Deconvolute and resize GRanges 
    if (any(GenomicRanges::strand(granges) == '*')) {
        if (stranded == T & split.bid.proms == T) {
            if (verbose) message('.. Bi-directional segments found and split...')
            granges <- deconvolveBidirectionalPromoters(granges)
        } else {
            if (verbose) message('.. Bi-directional promoters found but not split...')
        }
    } else {
        if (verbose) message('.. No bi-directional segments found. Continuing...')
    }
    if (center == T) {
        if (verbose) message('.. Centering GRanges...')
        granges <- IRanges::resize(granges, width = 1, fix = 'center')
    }
    if (length(flank) > 0) {
        if (verbose) message('.. Resizing GRanges...')
        granges <- IRanges::resize(granges, fix = 'end', width = flank[1])
        granges <- IRanges::resize(granges, fix = 'start', width = flank[1] + flank[2])
    }
    # Get matrix
    pattern <- Biostrings::DNAString(motif)
    seqs_g <- seqs[granges]
    if (!all(sapply(unlist(strsplit(as.character(pattern), '')), function(CHAR) {CHAR %in% c('A', 'T', 'C', 'G')}))) {
        ALGO <- "shift-or"
        FIXED <- FALSE
    }
    else {
        ALGO <- 'naive-exact'
        FIXED <- TRUE
    }
    hits <- Biostrings::vmatchPattern(pattern, seqs_g, algorithm = ALGO, fixed = FIXED) %>% 
        as("CompressedIRangesList") %>%
        IRanges::resize(1, fix = 'center')
    rhits <- Biostrings::vmatchPattern(pattern, seqs_g, algorithm = ALGO, fixed = FIXED) %>% 
        as("CompressedIRangesList") %>%
        resize(1, fix = 'center')
    out <- lapply(seq_along(granges), function(K) {
        if (as.vector(GenomicRanges::strand(granges))[K] == '+') {
            hits[[K]]
        } else {
            rhits[[K]]
        }
    }) %>% 
        as("CompressedIRangesList") %>%
        GenomicRanges::coverage(width = GenomicRanges::width(seqs_g)) %>% 
        sapply(., function(x) {
            as.numeric(x)
        }) %>% 
        t()
    return(out)
}

plotCoverageLine <- function(score, use.mean = TRUE, BIN = 10, plotEE = TRUE, XLIM = NULL, YLIM = NULL, XLAB = NULL, YLAB = NULL, COL = NULL, verbose = FALSE, return.extrema = FALSE, p = NULL, ...) {
    # Compute specific metrics to plot
    medians <- apply(score, 2, median)
    means <- apply(score, 2, mean)
    smoothed_line <- if (use.mean) {zoo::rollmean(means, BIN)} else {zoo::rollmean(medians, BIN)}
    stderror <- zoo::rollmean(apply(score, 2, function(n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n)) ) }), BIN)
    conint  <- zoo::rollmean(apply(score, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), BIN)
    topEE <- smoothed_line + conint
    bottomEE <- smoothed_line - conint
    MIN <- min(bottomEE) * 0.9
    MAX <- max(topEE) * 1.1
    df <- data.frame(
        coords = seq(round(-unique(length(smoothed_line) / 2)), round(unique(length(smoothed_line) / 2)))[1:length(smoothed_line)]+round(BIN/2),
        smoothed_line = smoothed_line,
        topEE = topEE,
        bottomEE = bottomEE
    )
    if (verbose) message('.. Stats computed')
    if (!return.extrema) {
        # Get other parameters
        YLIM <- if (!is.null(YLIM)) {YLIM} else if (plotEE) {c(min(bottomEE), max(topEE))} else {c(min(smoothed_line), max(smoothed_line))}
        XLIM <- if(!is.null(XLIM)) {XLIM} else {c(df$coords[1], df$coords[length(df$coords)])}
        XLAB <- ifelse(!is.null(XLAB), XLAB, 'Distance from beginning of GRanges')
        YLAB <- ifelse(!is.null(YLAB), YLAB, 'Score')
        COL <- ifelse(!is.null(COL), COL, '#000000')
        # Proceed to plotting
        if (is.null(p)) {
            p <- ggplot2::ggplot(data = df, ggplot2::aes(x = df$coords, y = smoothed_line)) + 
                ggplot2::geom_line(col = COL) + 
                ggplot2::geom_ribbon(aes(x = df$coords, ymin = bottomEE, ymax = topEE), alpha = 0.2, fill = COL, col = NA) +
                ggplot2::xlim(XLIM) + 
                ggplot2::ylim(YLIM) +
                ggplot2::theme_classic() + 
                ggplot2::xlab(XLAB) + 
                ggplot2::ylab(YLAB)
            return(p)
        } else {
            p <- p + 
                ggplot2::geom_line(data = df, ggplot2::aes(x = df$coords, y = smoothed_line), col = COL) + 
                ggplot2::geom_ribbon(data = df, ggplot2::aes(x = df$coords, ymin = bottomEE, ymax = topEE), alpha = 0.2, fill = COL, col = NA)
        }
    } else {
        return(c(MIN, MAX))
    }
} 

plotAggregateCoverage <- function(x, ...) {
    UseMethod("plotAggregateCoverage")
}
plotAggregateCoverage.default <- function(x, ...) {
    plotAggregateCoverage.GRanges(x, ...)
}
plotAggregateCoverage.GRanges <- function(
    granges, # The regions of interest
    list.bw.files, # Can be either a bw.file or an alreday imported Rle
    list.COL = NULL,
    split.bid.proms = TRUE,
    center = FALSE,
    flank = NULL,
    stranded = TRUE,
    use.mean = TRUE,
    BIN = 10,
    plotEE = TRUE, 
    auto.scale = c(0.05, 0.95), 
    scale.max = TRUE, 
    YLIM = NULL, 
    XLIM = NULL,
    XLAB = NULL,
    YLAB = NULL,
    verbose = TRUE,
    by.granges = TRUE,
    plot.central = TRUE,
    plot.legend = FALSE,
    ...) {
    
    .opts.plotAggregateCoverage <- function() {
        list.COL = NULL;
        split.bid.proms = TRUE;
        center = FALSE;
        flank = NULL;
        stranded = TRUE;
        use.mean = TRUE;
        BIN = 10;
        plotEE = TRUE;
        auto.scale = c(0.05, 0.95);
        scale.max = TRUE;
        YLIM = NULL; 
        XLIM = NULL;
        XLAB = NULL;
        YLAB = NULL;
        verbose = TRUE;
        by.granges = TRUE;
        plot.central = TRUE;
        plot.legend = FALSE;
    }
    .opts.plotAggregateCoverage.extra <- list(...)
    
    if (class(list.bw.files) == "SimpleRleList") { # Case where there is only 1 bigWig file and 1 set of granges
        # Define starting parameters
        bw.file <- list.bw.files
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        if (verbose) message('>>> Plotting granges')
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- rep(RColorBrewer::brewer.pal('Set2', n = 8), 5)
        }
        # Determining YLIM 
        if (verbose) message('.. Determining Y boundaries...')
        if (!is.null(YLIM)) {
            YLIM <- YLIM
        } else if (!is.null(auto.scale)) {
            max <- quantile(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE), auto.scale[2])
            min <- quantile(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE), auto.scale[1])
            YLIM <- c(min, max)
        } else if (scale.max) {
            max <- plotCoverageLine(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE), verbose = FALSE, return.extrema = TRUE)[2]
            min <- plotCoverageLine(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE),  verbose = FALSE, return.extrema = TRUE)[1]
            YLIM <- c(min, max)
        } else {
            stop('Please provide some Y limits. Aborting.')
        }
        if (verbose) message('.. YLIM: [', YLIM[1], '; ', YLIM[2], ']...')
        # Initiate plot
        if (verbose) message('.. Initiating graph...')
        if (is.null(XLIM)) XLIM <- c(1, flank[1] + flank[2])
        plot(
            NA, 
            bty = 'l',
            xlim = XLIM,
            ylim = YLIM,
            xlab = ifelse(length(XLAB) == 1, XLAB, XLAB[[k.g]]), 
            ylab = YLAB, 
            axes = F, 
            ...
        )
        box(bty='l')
        axis(2)
        axis(1, at = seq(XLIM[1], XLIM[2], length.out = 3), labels = c(-unique(GenomicRanges::width(granges) / 2), '0', unique(GenomicRanges::width(granges) / 2)))
        # Add central line
        if (plot.central) abline(v = seq(XLIM[1], XLIM[2], length.out = 3)[2], lty = 3)
        # Plot track
        if (verbose) message('.. Computing values...')
        mat <- getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE)
        par(new = T)
        if (verbose) message('.. Plot track...')
        p <- plotCoverageLine(mat, use.mean = use.mean, BIN = BIN, plotEE = plotEE, YLIM = YLIM, XLAB = XLAB, YLAB = YLAB, COL = list.COL[1], p = p)
        # Add central line
        if (plot.central) p <- p + geom_vline(xintercept = 0, colour="black", linetype = "longdash", size = 0.1)
        # Plot legend
        if (plot.legend == T) {
            par(new = T)
            legend('topleft', legend = names(list.bw.files), fill = list.COL, bty = 'n')
        }
        # Return plot
        return(p)
    } 
    else if (class(list.bw.files) == "list") { # Case where there are several bigWig files in a list
        # Define starting parameters
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        k.b <- 0
        p <- NULL
        # Detemine limits
        if (verbose) message('.. Determining Y boundaries...')
        if (!is.null(YLIM)) {
            YLIM <- YLIM
        } else if (!is.null(auto.scale)) {
            tmp <- lapply(list.bw.files, function(bw.file) {getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F)})
            max <- max(sapply(tmp, function(vec) {quantile(vec, auto.scale[2])}))
            min <- min(sapply(tmp, function(vec) {quantile(vec, auto.scale[1])}))
            YLIM <- c(min, max)
        } else if (scale.max) {
            max <- max(sapply(list.bw.files, function(bw.file) {plotCoverageLine(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[2]}))
            min <- min(sapply(list.bw.files, function(bw.file) {plotCoverageLine(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[1]}))
            YLIM <- c(min, max)
        } else if (!is.null(YLIM)) {
            YLIM <- YLIM
        } else {
            stop('Please provide some Y limits. Aborting.')
        }
        if (verbose) message('.. YLIM: [', YLIM[1], '; ', YLIM[2], ']...')
        if (is.null(XLIM)) XLIM <- c(1, (unique(GenomicRanges::width(granges)) - BIN))
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- rep(RColorBrewer::brewer.pal('Set2', n = 8), 5)
        }
        # Plotting
        for (k.b in 1:length(list.bw.files)) {
            if (verbose) message('>>> Plotting bw.file ', k.b, '/', length(list.bw.files))
            bw.file <- list.bw.files[[k.b]]
            # Determining matrix 
            mat <- getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE)
            # Plot
            p <- plotCoverageLine(
                mat, 
                use.mean = use.mean, 
                BIN = BIN, 
                plotEE = plotEE, 
                YLIM = YLIM, 
                XLAB = ifelse(length(XLAB) == 1, XLAB, XLAB[[k.b]]), 
                YLAB = YLAB, 
                COL = list.COL[k.b], 
                p = p
            )
        }
        # Add central line
        if (plot.central) p <- p + geom_vline(xintercept = 0, colour="black", linetype = "longdash", size = 0.1)
        # Plot legend
        if (plot.legend == T) {
            p <- p + 
                scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
                scale_colour_manual(
                    name = 'the colour', 
                    values =c('black'='black','red'='red'), 
                    labels = c('c2','c1')
                )
        }
        # Return plot
        return(p)
    }
    # End of loops
    message('All graphs successfully plotted')
}
plotAggregateCoverage.list <- function(list.granges, list.bw.files, by.granges = TRUE, ...) {
    if (by.granges == FALSE & length(list.bw.files) == 1) { # Case where there is only 1 bigWig file, that I want to plot over multiple sets of granges in the same graph
        # Fetch arguments from ellipsis
        .opts.ellipsis <- list(...)
        if (!is.null(.opts.ellipsis$list.COL)) {list.COL <- .opts.ellipsis$list.COL} else {list.COL <- TRUE}
        if (!is.null(.opts.ellipsis$split.bid.proms)) {split.bid.proms <- .opts.ellipsis$split.bid.proms} else {split.bid.proms <- TRUE}
        if (!is.null(.opts.ellipsis$center)) {center <- .opts.ellipsis$center} else {center <- FALSE}
        if (!is.null(.opts.ellipsis$flank)) {flank <- .opts.ellipsis$flank} else {flank <- NULL}
        if (!is.null(.opts.ellipsis$stranded)) {stranded <- .opts.ellipsis$stranded} else {stranded <- TRUE}
        if (!is.null(.opts.ellipsis$use.mean)) {use.mean <- .opts.ellipsis$use.mean} else {use.mean <- TRUE}
        if (!is.null(.opts.ellipsis$BIN)) {BIN <- .opts.ellipsis$BIN} else {BIN <- 10}
        if (!is.null(.opts.ellipsis$plotEE)) {plotEE <- .opts.ellipsis$plotEE} else {plotEE <- TRUE}
        if (!is.null(.opts.ellipsis$auto.scale)) {auto.scale <- .opts.ellipsis$auto.scale} else {auto.scale <- c(0.05, 0.95)}
        if (!is.null(.opts.ellipsis$scale.max)) {scale.max <- .opts.ellipsis$scale.max} else {scale.max <- TRUE}
        if (!is.null(.opts.ellipsis$YLIM)) {YLIM <- .opts.ellipsis$YLIM} else {YLIM <- NULL}
        if (!is.null(.opts.ellipsis$XLIM)) {XLIM <- .opts.ellipsis$XLIM} else {XLIM <- NULL}
        if (!is.null(.opts.ellipsis$YLAB)) {YLAB <- .opts.ellipsis$YLAB} else {YLAB <- NULL}
        if (!is.null(.opts.ellipsis$XLAB)) {XLAB <- .opts.ellipsis$XLAB} else {XLAB <- NULL}
        if (!is.null(.opts.ellipsis$verbose)) {verbose <- .opts.ellipsis$verbose} else {verbose <- TRUE}
        if (!is.null(.opts.ellipsis$plot.central)) {plot.central <- .opts.ellipsis$plot.central} else {plot.central <- TRUE}
        if (!is.null(.opts.ellipsis$plot.legend)) {plot.legend <- .opts.ellipsis$plot.legend} else {plot.legend <- FALSE}
        # Define starting parameters
        bw.file <- list.bw.files[[1]]
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        k.g <- 0
        p <- NULL
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- rep(RColorBrewer::brewer.pal('Set2', n = 8), 5)
        }
        # Determining YLIM 
        if (verbose) message('.. Determining Y boundaries...')
        if (!is.null(auto.scale)) {
            max <- max(sapply(list.granges, function(granges) quantile(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE), auto.scale[2])))
            min <- min(sapply(list.granges, function(granges) quantile(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE), auto.scale[1])))
            YLIM <- c(min, max)
        } 
        else if (scale.max) {
            max <- max(sapply(list.granges, function(granges) plotCoverageLine(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE), verbose = FALSE, return.extrema = TRUE)[2]))
            min <- min(sapply(list.granges, function(granges) plotCoverageLine(getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE),  verbose = FALSE, return.extrema = TRUE)[1]))
            YLIM <- c(min, max)
        } 
        else if (!is.null(YLIM)) {
            YLIM <- YLIM
        } 
        else {
            stop('Please provide some Y limits. Aborting.')
        }
        print(YLIM)
        # Plotting
        for (k.g in 1:length(list.granges)) {
            if (verbose) message('>>> Plotting granges ', k.g, '/', length(list.granges))
            granges <- list.granges[[k.g]]
            # Determining matrix 
            mat <- getCovMatrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = FALSE)
            # Plot
            p <- plotCoverageLine(mat, use.mean = use.mean, BIN = BIN, plotEE = plotEE, YLIM = YLIM, XLAB = XLAB, YLAB = YLAB, COL = list.COL[k.g], p = p)
        }
        # Add central line
        if (plot.central) p <- p + geom_vline(xintercept = 0, colour="black", linetype = "longdash", size = 0.1)
        # Plot legend
        if (plot.legend == T) {
            p <- p + 
                scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
                scale_colour_manual(
                    name = 'the colour', 
                    values =c('black'='black','red'='red'), 
                    labels = c('c2','c1')
                )
        }
        # Return plot
        return(p)
    } 
    else { # Case where there is 1 or more bigWig files, that I want to plot over multiple sets of granges separately
        lapply(list.granges, function(granges) {
            plotAggregateCoverage(
                granges,
                list.bw.files, 
                ...
            )
        })
    }
}

plotHeatmaps <- function(
    MAT, 
    colorspace = NULL, 
    autoscale = TRUE, 
    s = 0.05, 
    zmin = 0, 
    zmax = 10, 
    titles=rep('', length(MAT)),
    xlab='',
    ylab="", 
    main='', 
    cex.lab=12.0, 
    cex.axis=12.0, 
    cex.legend=12.0, 
    pointsize=12,
    cex.title=20,
    ...
) 
{
    if (class(MAT) == 'matrix') MAT <- list(MAT)
    if (length(colorspace)) {
        gcol <- colorRampPalette(colorspace)
    } else {
        gcol <- colorRampPalette(c(
            "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
            "#FF7F00", "red", "#7F0000"
        ))     
    }
    NP <- length(MAT)
    datapoints <- unlist(MAT)
    min <- min(datapoints, na.rm = TRUE)
    max <- max(datapoints, na.rm=TRUE) 
    plots <- list()
    for (i in seq(NP)) {
        data <- MAT[[i]]
        bins <- 1:ncol(data) 
        colnames(data) <- bins
        # Rescale outliers
        if (autoscale) {
            zlim <- quantile(data, c(s,1-s), na.rm=TRUE)
            zmin <- zlim[1]
            zmax <- zlim[2]
        } 
        data[data < zmin] <- zmin
        data[data > zmax] <- zmax
        # Define color scale limits
        keycolor_lim <- range(data, na.rm = TRUE)
        keycolor_lim[1] <- zmin 
        keycolor_lim[2] <- zmax 
        # Plot graph
        p <- ggplot(reshape2::melt(data), aes(Var2, Var1, fill = value))
        p <- p + geom_raster()
        p <- p + scale_fill_gradientn(
            colours = gcol(50), limits = keycolor_lim, 
            breaks = keycolor_lim, labels = format(keycolor_lim, digits = 2)
        )
        p <- p + scale_x_continuous(
            breaks = c(min(bins), max(bins)),
            labels = paste0(c('-', ''), (length(bins) / 2)),
            expand = c(0.05, 0.05)
        ) 
        p <- p + geom_vline(
            xintercept = length(bins)/2, size = .5, colour = 'black', linetype = 'dashed'
        ) 
        p <- p + scale_y_reverse(
            breaks=c(1, nrow(data)),
            labels=c(1, nrow(data)),
            expand = c(0, 0)
        ) 
        p <- p + ggtitle(titles[i]) 
        p <- p + xlab(xlab) 
        p <- p + ylab(ylab) 
        p <- p + theme(
            legend.position = "bottom", 
            axis.text=element_text(size=cex.axis), 
            axis.title=element_text(size=cex.lab),
            title=element_text(size=cex.lab),
            legend.text=element_text(size=cex.legend) 
        ) 
        p <- p + guides(
            fill = guide_colorbar(frame.colour = 'black', frame.size = 0.25, title = "", raster = TRUE)
        )
        p <- p + theme(
            panel.background = element_blank(), 
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_blank(),
        ) 
        p <- p + geom_rect(
            aes(
                xmin = min(Var2), 
                xmax = max(Var2), 
                ymin = min(Var1), 
                ymax = max(Var1)
            ), 
            fill = NA, color = 'black', size = 0.5
        )
        plots[[i]] <- p
    }
    if (length(plots) == 1) plots <- plots[[1]]
    return(plots)
}

num2bp <- function( n ) {
    if(n == 0)  return(paste0(0, 'bp'))
    
    number <- abs(n)
    if(number == n) sign <- '' else sign <- '-'
    
    lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06, 
             0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21, 
             1e+24)
    pre <- c("yb", "zb", "ab", "fb", "pb", "nb", "ub", "mb", "", "kb", 
             "Mb", "Gb", "Tb", "Pb", "Eb", "Zb", "Yb")
    ix <- findInterval(number, lut)
    if (lut[ix]!=1) {
        sistring <- paste0(
            sign, format(round(number/lut[ix], 3), digits = 3), pre[ix]
        )

    } else {
        sistring <- paste0(format(n, digits=3), 'bp')
    }
    return(sistring)
}
