#### ---- GRanges utils ---- ####

#' Core function
#'
#' @param granges A GRanges object with a TSS column or TSS.rev and TSS.fwd columns
#' @param upstream How many bases upstream of the TSS?
#' @param downstream How many bases downstream of the TSS?
#' 
#' @export
#' @return GRanges aligned to the TSS column or to TSS.rev and TSS.fwd columns

alignToTSS <- function(granges, upstream, downstream) {
    if (any(GenomicRanges::strand(granges) == '*')) {
        granges <- deconvolveBidirectionalPromoters(granges)
    }
    if (!is.null(granges$TSS)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(as.vector(GenomicRanges::strand(granges)) == '+', (granges$TSS - upstream), (granges$TSS - downstream + 1)),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else if (!is.null(granges$TSS.fwd) & !is.null(granges$TSS.rev)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(as.vector(GenomicRanges::strand(granges)) == '+', (granges$TSS.fwd - upstream), (granges$TSS.rev - downstream + 1)),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else {
        error("No TSS column found. Aborting.")
    }
    return(granges)
}

#' Core function
#'
#' @param granges A GRanges object 
#' 
#' @export
#' @return GRanges with only '+' and '-' strands. GRanges with '*' strand 
#' have been duplicated and split into forward and reverse strands.

deconvolveBidirectionalPromoters <- function(granges) {
    unid <- granges[GenomicRanges::strand(granges) == '+' | GenomicRanges::strand(granges) == '-']
    bid <- granges[GenomicRanges::strand(granges) == '*']
    bid.fwd <- bid
    GenomicRanges::strand(bid.fwd) <- '+'
    bid.rev <- bid
    GenomicRanges::strand(bid.rev) <- '-'
    granges_shifted <- sort(c(unid, bid.fwd, bid.rev), ignore.strand = TRUE)
    return(granges_shifted)
}

#' Core function
#'
#' @param granges A GRanges object 
#' @param genome DNAStringSet object. Ideally, the sequence of an entire genome, 
#' obtained for instance by running 
#' \code{Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11)}.
#' 
#' @export
#' @return GRanges with a seq column containing the sequence of the GRanges.

withSeq <- function(granges, genome) {
    granges$seq <- genome[granges]
    return(granges)
}

#### ---- Line coverage plot ---- ####

#' Internal function
#'
#' @param granges A GRanges to map a track onto
#' @param bw_file a track in RleList format
#' @param norm Should the signal be normalized ('none', 'zscore' or 'log2')?
#' @param verbose Boolean
#' 
#' @export
#' @return A square numerical matrix of signal values over the GRanges

getCovMatrix <- function(granges, bw_file, norm = 'none', verbose = TRUE) {
    scores <- bw_file
    # Subset scores over GRanges
    granges <- subsetByOverlaps(granges, GRanges(seqlevels(bw_file), IRanges(1, width = lengths(bw_file))))
    scores.subset <- scores[granges]
    # Turn it into a rectangular matrix and flip reverse strand scores
    if (verbose) message('.. Converting scores into matrix...')
    scores.subset <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(scores.subset)), nrow = length(granges), byrow = T)))
    scores.subset.flipped <- t(apply(scores.subset, 1, rev))
    scores.subset <- matrix(sapply(1:nrow(scores.subset), function(K) {if((as.vector(GenomicRanges::strand(granges)) == '-')[K]) {scores.subset.flipped[K,]} else {scores.subset[K,]}}), nrow = length(granges), byrow = T)
    # Normalize matrix
    if (norm == 'zscore') {
        scores.subset <- apply(scores.subset, 1, scale) %>% t() %>% na.replace(0)
    } 
    else if (norm == 'log2') {
        scores.subset <- log2(scores.subset)
    }
    # Return matrix
    return(scores.subset)
}

#' Core function
#'
#' @param x a single signal track (in SimpleRleList or CompressedRleList class) 
#' or a set of signal tracks (in list class)
#' 
#' @export
#' @return A plot

plotAggregateCoverage <- function(x, ...) {
    UseMethod("plotAggregateCoverage")
}

#' Core function
#'
#' @param bw a single signal track CompressedRleList class
#' 
#' @return A plot

plotAggregateCoverage.CompressedRleList <- function(bw, granges, ...) {
    bw_file <- as(bw_file, 'SimpleRleList')
    plotAggregateCoverage(bw_file, granges, ...)
}

#' Core function
#'
#' @param bw a single signal track SimpleRleList class
#' @param granges a GRanges object or a list of GRanges
#' @param colors a vector of colors
#' @param xlab x axis label
#' @param ylab y axis label
#' @param xlim y axis limits
#' @param ylim y axis limits
#' @param quartiles Which quantiles to use to determine y scale automatically?
#' @param verbose Boolean
#' @param bin Integer Width of the window to use to smooth values
#' @param plot_central Boolean Draw a vertical line at 0
#' 
#' @return A plot

plotAggregateCoverage.SimpleRleList <- function(
    bw, 
    granges, 
    colors = c('#991919', '#1232D9', '#3B9B46', '#D99B12', '#7e7e7e', '#D912D4', '#9E7FCC', '#B0E797', '#D1B3B3', '#23A4A3', '#000000'), 
    xlab = 'Center of elements',
    ylab = 'Score', 
    xlim = NULL, 
    ylim  = NULL, 
    quartiles = c(0.025, 0.975),
    verbose = FALSE,
    bin = 1,
    plot_central = TRUE,
    ...
) 
{
    if (class(granges) == 'GRanges') granges <- list('granges' = granges)
    # Compute scores
    lists <- lapply(seq_along(granges), function(K) {
        g <- granges[[K]]
        if (length(unique(width(g))) > 1) {
            error('Please provide GRanges that are all the same width. Aborting.')
        }
        mat <- getCovMatrix(g, bw, verbose = FALSE)
        colnames(mat) <- c(
            -unique(width(g))/2, 
            rep('', (ncol(mat) - 3)/2), 
            '1', 
            rep('', (ncol(mat) - 2)/2), 
            paste0('+', unique(width(g))/2)
        )
        row.names(mat) <- paste0('locus_', as.character(1:nrow(mat)))
        means <- c(rep(NA, bin/2), zoo::rollmean(apply(mat, 2, mean), bin), rep(NA, bin/2))[1:ncol(mat)]
        conint <- c(rep(NA, bin/2), zoo::rollmean(apply(mat, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), bin), rep(NA, bin/2))[1:ncol(mat)]
        topEE <- means + conint
        bottomEE <- means - conint
        coords <- round(seq(-unique(width(g))/2, unique(width(g))/2, length.out = 1000))
        EE <- data.frame(
            means = means, 
            meansUp = topEE, 
            meansDown = bottomEE, 
            x = coords,
            grange = rep(ifelse(!is.null(names(granges)), names(granges)[K], as.character(K)), length(means)), 
            stringsAsFactors = FALSE
        )
        return(EE)
    }) 
    EEs <- do.call(rbind, lists)
    if (!is.null(names(granges))) {
        EEs$grange <- factor(EEs$grange, levels = names(granges))
    } 
    else {
        EEs$grange <- factor(EEs$grange, levels = as.character(1:length(granges)))
    }
    # Determining YLIM 
    if (is.null(ylim)) {
        max <- max(EEs$meansUp) * 1.05
        min <- min(EEs$meansDown) * 0.95
        ylim <- c(min, max)
    }
    # Start plotting
    require(ggplot2)
    p <- ggplot(data = EEs, aes(
        x = x, 
        y = means, 
        ymin = meansDown, 
        ymax = meansUp, 
        group = grange, 
        col = grange, 
        fill = grange
    ))
    p <- p + geom_line()
    p <- p + geom_ribbon(alpha = 0.2, col = NA)
    p <- p + theme_bw()
    p <- p + scale_x_continuous(
        limits = xlim, 
        expand = c(0, 0)
    )
    p <- p + scale_y_continuous(
        limits = ylim, 
        expand = c(0, 0)
    )
    p <- p + labs(x = xlab, y = ylab, col = 'Sets:')
    p <- p + guides(fill = FALSE)
    if (length(granges) == 1) p <- p + guides(color = FALSE)
    if (length(granges) <= length(colors)) {
        p <- p + scale_fill_manual(values = colors)
        p <- p + scale_color_manual(values = colors)
    }
    if (plot_central) p <- p + geom_vline(xintercept = 0, colour="black", linetype = "longdash", size = 0.1)
    # Return plot
    return(p)
}

#' Core function
#'
#' @param bw several signal tracks (SimpleRleList or CompressedRleList class) into
#' a class object
#' @param granges a GRanges object or a list of GRanges
#' @param colors a vector of colors
#' @param xlab x axis label
#' @param ylab y axis label
#' @param xlim y axis limits
#' @param ylim y axis limits
#' @param quartiles Which quantiles to use to determine y scale automatically?
#' @param verbose Boolean
#' @param bin Integer Width of the window to use to smooth values
#' @param plot_central Boolean Draw a vertical line at 0
#' @param split_by_granges Boolean Split plots over the sets of GRanges
#' @param split_by_granges Boolean Split plots over the sets of GRanges
#' @param split_by_track Boolean Split plots over the sets of signal tracks
#' 
#' @return A plot

plotAggregateCoverage.list <- function(
    bw_list, 
    granges, 
    colors = c('#991919', '#1232D9', '#3B9B46', '#D99B12', '#7e7e7e', '#D912D4', '#9E7FCC', '#B0E797', '#D1B3B3', '#23A4A3', '#000000'), 
    xlab = 'Center of elements',
    ylab = 'Score', 
    xlim = NULL, 
    ylim  = NULL, 
    quartiles = c(0.025, 0.975),
    verbose = FALSE,
    bin = 1,
    plot_central = TRUE,
    split_by_granges = TRUE,
    split_by_track = FALSE,
    ...
) 
{
    if (!all(unlist(lapply(bw_list, function(L) class(L) %in% c('SimpleRleList', 'CompressedRleList'))))) 
        error('Some objects in the bw_list are not bigwig tracks (in SimpleRleList or CompressedRleList format). Aborting.')
    llists <- lapply(seq_along(bw_list), function(B) {
        bw <- bw_list[[B]]
        # Compute scores
        lists <- lapply(seq_along(granges), function(K) {
            g <- granges[[K]]
            if (length(unique(width(g))) > 1) {
                error('Please provide GRanges that are all the same width. Aborting.')
            }
            mat <- getCovMatrix(g, bw, verbose = FALSE)
            colnames(mat) <- c(
                -unique(width(g))/2, 
                rep('', (ncol(mat) - 3)/2), 
                '1', 
                rep('', (ncol(mat) - 2)/2), 
                paste0('+', unique(width(g))/2)
            )
            row.names(mat) <- paste0('locus_', as.character(1:nrow(mat)))
            means <- c(rep(NA, bin/2), zoo::rollmean(apply(mat, 2, mean), bin), rep(NA, bin/2))[1:ncol(mat)]
            conint <- c(rep(NA, bin/2), zoo::rollmean(apply(mat, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), bin), rep(NA, bin/2))[1:ncol(mat)]
            topEE <- means + conint
            bottomEE <- means - conint
            coords <- round(seq(-unique(width(g))/2, unique(width(g))/2, length.out = 1000))
            EE <- data.frame(
                means = means, 
                meansUp = topEE, 
                meansDown = bottomEE, 
                x = coords,
                grange = rep(ifelse(!is.null(names(granges)), names(granges)[K], as.character(K)), length(means)), 
                bw = rep(ifelse(!is.null(names(bw_list)), names(bw_list)[B], as.character(K)), length(means)), 
                stringsAsFactors = FALSE
            )
            return(EE)
        }) 
    })
    EEs <- do.call(rbind, do.call(rbind, llists))
    if (!is.null(names(granges))) {
        EEs$grange <- factor(EEs$grange, levels = names(granges))
    } 
    else {
        EEs$grange <- factor(EEs$grange, levels = as.character(1:length(granges)))
    }
    if (!is.null(names(bw_list))) {
        EEs$bw <- factor(EEs$bw, levels = names(bw_list))
    } 
    else {
        EEs$bw <- factor(EEs$bw, levels = as.character(1:length(bw_list)))
    }
    # Determining YLIM 
    if (is.null(ylim)) {
        max <- max(EEs$meansUp) * 1.05
        min <- min(EEs$meansDown) * 0.95
        ylim <- c(min, max)
    }
    # Start plotting
    require(ggplot2)
    if (split_by_track) split_by_granges <- FALSE
    if (split_by_granges) split_by_track <- FALSE
    if (split_by_granges & split_by_track) {
        p <- ggplot(data = EEs, aes(
            x = x, 
            y = means, 
            ymin = meansDown, 
            ymax = meansUp, 
            group = grange, 
            col = grange, 
            fill = grange
        ))
        p <- p + labs(x = xlab, y = ylab, col = 'Sets:')
        p <- p + facet_grid(bw~grange)
    } 
    else if (!split_by_granges & split_by_track) {
        p <- ggplot(data = EEs, aes(
            x = x, 
            y = means, 
            ymin = meansDown, 
            ymax = meansUp, 
            group = grange, 
            col = grange, 
            fill = grange
        ))
        p <- p + labs(x = xlab, y = ylab, col = 'Sets:')
        p <- p + facet_grid(~bw)
    } 
    else if (split_by_granges & !split_by_track) {
        p <- ggplot(data = EEs, aes(
            x = x, 
            y = means, 
            ymin = meansDown, 
            ymax = meansUp, 
            group = bw, 
            col = bw, 
            fill = bw
        ))
        p <- p + labs(x = xlab, y = ylab, col = 'Signal:')
        p <- p + facet_grid(~grange)
    }
    p <- p + geom_line()
    p <- p + geom_ribbon(alpha = 0.2, col = NA)
    p <- p + theme_bw()
    p <- p + scale_x_continuous(
        limits = xlim, 
        expand = c(0, 0)
    )
    p <- p + scale_y_continuous(
        limits = ylim, 
        expand = c(0, 0)
    )
    p <- p + guides(fill = FALSE)
    if (length(granges) == 1) p <- p + guides(color = FALSE)
    if (length(granges) <= length(colors)) {
        p <- p + scale_fill_manual(values = colors)
        p <- p + scale_color_manual(values = colors)
    }
    if (plot_central) p <- p + geom_vline(xintercept = 0, colour="black", linetype = "longdash", size = 0.1)
    p <- p + theme(panel.spacing = unit(2, "lines"))
    # Return plot
    return(p)
}

#### ---- Heatmap coverage plot ---- ####

#' Core function - NOT WORKING YET
#'
#' @export

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
    BIN = 1,
    seriate = FALSE, 
    clusters = 0,
    ...
) 
{
    require(ggplot2)
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
    plots <- list()
    for (i in seq(NP)) {
        data <- MAT[[i]]
        bins <- 1:ncol(data)
        colnames(data) <- 1:ncol(data) 
        # Bin data
        if (BIN > 1) {
            data <- apply(data, 1, function(ROW) {
                zoo::rollmean(ROW, BIN)
            }) %>% t()
            bins <- 1:ncol(data) + BIN/2 - 1
            colnames(data) <- bins 
        }
        # Rescale outliers
        if (autoscale) {
            zlim <- quantile(data, c(s,1-s), na.rm=TRUE)
            zmin <- zlim[1]
            zmax <- zlim[2]
        } 
        data[data < zmin] <- zmin
        data[data > zmax] <- zmax
        # Seriate data
        if (seriate & clusters == 0) {
            data <- data[seriation::get_order(seriation::seriate(data)),]
        }
        # Cluster data
        if (clusters > 0 & !seriate) {
            set.seed(222)
            data <- data[order(kmeans(data, centers = clusters)$cluster),]
        }
        # Define color scale limits
        keycolor_lim <- range(data, na.rm = TRUE)
        keycolor_lim[1] <- zmin 
        keycolor_lim[2] <- zmax 
        # Plot graph
        df <- reshape2::melt(data)
        rect_df <- data.frame(
            xmin = min(df$Var2), 
            xmax = max(df$Var2), 
            ymin = min(df$Var1), 
            ymax = max(df$Var1)
        )
        p <- ggplot(df, aes(Var2, Var1, fill = value))
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
            data = rect_df,
            mapping = aes(
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax
            ), 
            inherit.aes = FALSE,
            fill = NA, 
            color = 'black', 
            size = 0.5
        )
        plots[[i]] <- p
    }
    if (length(plots) == 1) plots <- plots[[1]]
    return(plots)
}

#### ---- Motif mapping ---- ####

#' Core function - NOT WORKING YET
#' 
#' @export

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

