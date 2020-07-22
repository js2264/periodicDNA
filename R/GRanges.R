#' A function to re-align a GRanges object to TSSs
#' 
#' This function re-aligns ranges (typically regulatory elements)
#' to a set of coordinates, either the TSS column or the
#' TSS.fwd and TSS.rev columns. If none are found, the function
#' assumes the ranges are promoters and that the end or the ranges
#' are the TSSs. 
#'
#' @param granges A stranded GRanges object with a TSS column 
#' or TSS.rev and TSS.fwd columns
#' @param upstream How many bases upstream of the TSS should the GRanges
#' object by extended by? [Default: 0]
#' @param downstream How many bases downstream of the TSS should the GRanges
#' object by extended by? [Default: 1]
#' @return GRanges aligned to the TSS column or to TSS.rev 
#' and TSS.fwd columns, and extended by upstream/downstream bp. 
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import BSgenome
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' ce11_proms
#' alignToTSS(ce11_proms)

alignToTSS <- function(granges, upstream = 0, downstream = 1) {
    if (any(GenomicRanges::strand(granges) == '*')) {
        granges <- deconvolveBidirectionalPromoters(granges)
    }
    if (!is.null(granges$TSS)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(
                as.vector(GenomicRanges::strand(granges)) == '+',
                (granges$TSS - upstream), 
                (granges$TSS - downstream + 1)
            ),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else if (!is.null(granges$TSS.fwd) & !is.null(granges$TSS.rev)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(
                as.vector(GenomicRanges::strand(granges)) == '+', 
                (granges$TSS.fwd - upstream), 
                (granges$TSS.rev - downstream + 1)
            ),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else {
        TSSs <- GenomicRanges::start(
            GenomicRanges::resize(granges, fix = 'end', width = 1)
        )
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(
                as.vector(GenomicRanges::strand(granges)) == '+', 
                (TSSs - upstream), 
                (TSSs - downstream + 1)
            ),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    return(granges)
}

#' A function to duplicate bi-directional GRanges
#' 
#' This function splits bi-directional ranges into + and - 
#' stranded ranges. It duplicates the ranges which are '*'.
#'
#' @param granges A stranded GRanges object 
#' @return GRanges with only '+' and '-' strands. GRanges with '*' strand 
#' have been duplicated and split into forward and reverse strands.
#' 
#' @import GenomicRanges
#' @export
#' 
#' @examples
#' data(ce11_all_REs)
#' library(GenomicRanges)
#' proms <- ce11_all_REs[grepl('prom', ce11_all_REs$regulatory_class)]
#' proms
#' table(strand(proms))
#' proms <- deconvolveBidirectionalPromoters(proms)
#' proms
#' table(strand(proms))

deconvolveBidirectionalPromoters <- function(granges) {
    filt <- GenomicRanges::strand(granges) == '+' | 
        GenomicRanges::strand(granges) == '-'
    unid <- granges[filt]
    bid <- granges[GenomicRanges::strand(granges) == '*']
    bid.fwd <- bid
    GenomicRanges::strand(bid.fwd) <- '+'
    bid.rev <- bid
    GenomicRanges::strand(bid.rev) <- '-'
    granges_shifted <- sort(c(unid, bid.fwd, bid.rev), ignore.strand = TRUE)
    return(granges_shifted)
}

#' A function to add sequence to an existing GRanges
#' 
#' This function is 'tidyverse-friendly', i.e. it adds a .$seq 
#' columns to a GRanges of interest, according the provided BSgenome
#' or DNAStringSet object. 
#'
#' @param granges A GRanges object 
#' @param genome UCSC ID, BSgenome of DNAStringSet object. 
#' @return GRanges with a seq column containing the sequence of the GRanges.
#' 
#' @importFrom methods is
#' @import GenomicRanges
#' @import Biostrings
#' @export
#' 
#' @examples
#' data(ce11_all_REs)
#' ce11_all_REs
#' ce11_all_REs <- withSeq(ce11_all_REs, 'BSgenome.Celegans.UCSC.ce11')
#' ce11_all_REs
#' ce11_all_REs$seq

withSeq <- function(granges, genome) {
    if (methods::is(genome, 'character')) {
        if (genome %in% c(
            'sacCer3', 'ce11', 'dm6', 'mm10', 'hg38', 'danRer10'
        ) | genome %in% BSgenome::installed.genomes()) {
            genome <- BSgenome::getBSgenome(genome)
        }
        else {
            return(stop(
                'Only sacCer3, ce11, dm6, mm10, hg38 
                and danRer10 are supported'
            ))
        }
    }
    if (methods::is(genome, 'BSgenome')) {
        genome <- Biostrings::getSeq(genome)
    }
    granges$seq <- genome[granges]
    return(granges)
}

getCovMatrix <- function(g, bw, norm = 'none', verbose = FALSE) {
    scores <- bw
    # Subset scores over GRanges
    g <- IRanges::subsetByOverlaps(
        g, 
        GenomicRanges::GRanges(
            GenomeInfoDb::seqlevels(bw), 
            IRanges::IRanges(1, width = lengths(bw))
        ), 
        type = 'within'
    )
    scores.subset <- scores[g]
    # Turn it into a rectangular matrix and flip reverse strand scores
    if (verbose) message('.. Converting scores into matrix...')
    scores.subset <- suppressWarnings(suppressMessages(
        matrix(
            as.vector(unlist(scores.subset)), nrow = length(g),
            byrow = TRUE
        )
    ))
    scores.subset.flipped <- t(apply(scores.subset, 1, rev))
    scores.subset <- matrix(
        lapply(
            seq_len(nrow(scores.subset)), 
            function(K) {
                if((as.vector(GenomicRanges::strand(g)) == '-')[K]) {
                    scores.subset.flipped[K,]
                } else {
                    scores.subset[K,]
                }
            }
        ) %>% unlist(), 
        nrow = length(g), 
        byrow = TRUE
    )
    # Normalize matrix
    if (norm == 'zscore') {
        scores.subset <- apply(
            scores.subset, 1, scale
        ) %>% t() %>% na.replace(0)
    } 
    else if (norm == 'log2') {
        scores.subset <- log2(scores.subset + 1)
    }
    # Return matrix
    return(scores.subset)
}

#' A function to plot aggregated signals over sets of GRanges
#'
#' This function takes one or several RleList genomic tracks (imported
#' by rtraklayer::import(..., as = 'Rle')) and one or several GRanges 
#' objects. It computes coverage of the GRanges by the genomic tracks
#' and returns an aggregate coverage plot. 
#'
#' @param x CompressedRleList or SimpleRleList (in lists or single)
#' @param ... additional parameters
#' @return An aggregate coverage plot. 
#' 
#' @export
#' 
#' @examples
#' data(ce11_ATACseq)
#' data(ce11_proms)
#' p <- plotAggregateCoverage(
#'     ce11_ATACseq, 
#'     resize(ce11_proms[1:100], fix = 'center', width = 1000)
#' )
#' p

plotAggregateCoverage <- function(x, ...) {
    UseMethod("plotAggregateCoverage")
}

#' A function to plot aggregated signals over sets of GRanges
#'
#' This function takes one or several RleList genomic tracks (imported
#' by rtraklayer::import(..., as = 'Rle')) and one or several GRanges 
#' objects. It computes coverage of the GRanges by the genomic tracks
#' and returns an aggregate coverage plot. 
#'
#' @param x a single signal track (CompressedRleList class)
#' @param granges a GRanges object or a list of GRanges
#' @param ... additional parameters
#' @return An aggregate coverage plot. 
#' 
#' @importFrom methods as
#' @export
#' 
#' @examples
#' data(ce11_ATACseq)
#' data(ce11_proms)
#' class(ce11_ATACseq)
#' p <- plotAggregateCoverage(
#'     ce11_ATACseq, 
#'     resize(ce11_proms[1:100], fix = 'center', width = 1000)
#' )
#' p

plotAggregateCoverage.CompressedRleList <- function(
    x,
    granges, 
    ...
)
{
    bw <- methods::as(x, 'SimpleRleList')
    plotAggregateCoverage(bw, granges, ...)
}

#' A function to plot aggregated signals over sets of GRanges
#'
#' This function takes one or several RleList genomic tracks (imported
#' by rtraklayer::import(..., as = 'Rle')) and one or several GRanges 
#' objects. It computes coverage of the GRanges by the genomic tracks
#' and returns an aggregate coverage plot. 
#'
#' @param x a single signal track (SimpleRleList class)
#' @param granges a GRanges object or a list of GRanges
#' @param colors a vector of colors
#' @param xlab x axis label
#' @param ylab y axis label
#' @param xlim y axis limits
#' @param ylim y axis limits
#' @param quartiles Which quantiles to use to determine y scale automatically?
#' @param verbose Boolean
#' @param bin Integer Width of the window to use to smooth values 
#' by zoo::rollMean
#' @param plot_central Boolean Draw a vertical line at 0
#' @param run_in_parallel Boolean Should the plots be computed in parallel
#' using mclapply?
#' @param split_by_granges Boolean Facet plots over the sets of GRanges
#' @param norm character Should the signal be normalized 
#' ('none', 'zscore' or 'log2')?
#' @param ... additional parameters
#' @return A plot of aggregated signals
#' 
#' @import GenomicRanges
#' @import ggplot2
#' @importFrom zoo rollmean
#' @importFrom parallel mclapply
#' @importFrom stats qt
#' @export
#' 
#' @examples
#' data(ce11_ATACseq)
#' data(ce11_proms)
#' class(ce11_ATACseq)
#' p1 <- plotAggregateCoverage(
#'     ce11_ATACseq, 
#'     resize(ce11_proms[1:100], fix = 'center', width = 1000)
#' )
#' p1
#' proms <- resize(ce11_proms[1:100], fix = 'center', width = 500)
#' p2 <- plotAggregateCoverage(
#'     ce11_ATACseq, 
#'     list(
#'         'Ubiq & Germline promoters' = 
#'             proms[proms$which.tissues %in% c('Ubiq.', 'Germline')],
#'         'Other promoters' = 
#'             proms[!(proms$which.tissues %in% c('Ubiq.', 'Germline'))]
#'     )
#' )
#' p2

plotAggregateCoverage.SimpleRleList <- function(
    x, 
    granges, 
    colors = NULL, 
    xlab = 'Center of elements',
    ylab = 'Score', 
    xlim = NULL, 
    ylim  = NULL, 
    quartiles = c(0.025, 0.975),
    verbose = FALSE,
    bin = 1,
    plot_central = TRUE,
    run_in_parallel = FALSE,
    split_by_granges = FALSE,
    norm = 'none',
    ...
) 
{
    bw <- x
    if (is.null(colors)) {
        colors <- rep(c(
            '#991919', 
            '#1232D9', 
            '#3B9B46', 
            '#D99B12', 
            '#7e7e7e', 
            '#D912D4', 
            '#9E7FCC', 
            '#B0E797', 
            '#D1B3B3', 
            '#23A4A3', 
            '#000000'
        ), 5)
    }
    
    if (methods::is(granges, 'GRanges')) granges <- list('granges' = granges)
    # Compute scores
    lists <- parallel::mclapply(seq_along(granges), function(K) {
        g <- granges[[K]]
        if (length(unique(GenomicRanges::width(g))) > 1) {
            stop('Please provide only GRanges that 
            are all the same width. Aborting.')
        }
        if (unique(GenomicRanges::width(g)) < 2) {
            stop('Please provide only GRanges with widths >=2. Aborting.')
        }
        mat <- getCovMatrix(g, bw, verbose = FALSE, norm = norm)
        colnames(mat) <- c(
            -unique(GenomicRanges::width(g))/2, 
            rep('', (ncol(mat) - 3)/2), 
            '1', 
            rep('', (ncol(mat) - 2)/2), 
            paste0('+', unique(GenomicRanges::width(g))/2)
        )
        row.names(mat) <- paste0('locus_', as.character(seq_len(nrow(mat))))
        means <- c(
            rep(NA, bin/2), 
            zoo::rollmean(apply(mat, 2, mean), bin), 
            rep(NA, bin/2)
        )[seq_len(ncol(mat))]
        conint <- c(
            rep(NA, bin/2), 
            zoo::rollmean(
                apply(mat, 2, function (n) {
                    Q <- stats::qt(0.975,sum(!is.na(n)))
                    S <- sd(n,na.rm=TRUE)
                    sq <- sqrt(sum(!is.na(n)))
                    Q * S / sq
                }), 
                bin
            ), 
            rep(NA, bin/2)
        )[seq_len(ncol(mat))]
        topEE <- means + conint
        bottomEE <- means - conint
        coords <- round(
            seq(
                -unique(GenomicRanges::width(g))/2, 
                unique(GenomicRanges::width(g))/2, 
                length.out = unique(width(g))
            )
        )
        EE <- data.frame(
            'means' = means, 
            'meansUp' = topEE, 
            'meansDown' = bottomEE, 
            'x' = coords,
            'grange' = rep(
                ifelse(
                    !is.null(names(granges)), 
                    names(granges)[K], 
                    as.character(K)
                ), 
                length(means)
            ), 
            stringsAsFactors = FALSE
        )
        return(EE)
    }, mc.cores = ifelse(run_in_parallel, length(granges), 1)) 
    EEs <- do.call(rbind, lists)
    if (!is.null(names(granges))) {
        EEs$grange <- factor(EEs$grange, levels = names(granges))
    } 
    else {
        EEs$grange <- factor(
            EEs$grange, levels = as.character(seq_along(granges))
        )
    }
    # Determining YLIM 
    if (is.null(ylim)) {
        max <- max(na.remove(EEs$meansUp)) * 1.05
        min <- min(na.remove(EEs$meansDown)) * 0.95
        ylim <- c(min, max)
    }
    # Start plotting
    p <- ggplot2::ggplot(data = EEs, ggplot2::aes(
        x = x, 
        y = means, 
        ymin = meansDown, 
        ymax = meansUp, 
        group = grange, 
        col = grange, 
        fill = grange
    ))
    p <- p + ggplot2::geom_line()
    p <- p + ggplot2::geom_ribbon(alpha = 0.2, col = NA)
    p <- p + theme_ggplot2()
    p <- p + ggplot2::scale_x_continuous(
        limits = xlim, 
        expand = c(0, 0)
    )
    p <- p + ggplot2::scale_y_continuous(
        limits = ylim, 
        expand = c(0, 0)
    )
    p <- p + ggplot2::labs(x = xlab, y = ylab, col = 'Sets:')
    p <- p + ggplot2::guides(fill = FALSE)
    if (length(granges) == 1) p <- p + ggplot2::guides(color = FALSE)
    if (length(granges) <= length(colors)) {
        p <- p + ggplot2::scale_fill_manual(values = colors)
        p <- p + ggplot2::scale_color_manual(values = colors)
    }
    if (plot_central) p <- p + ggplot2::geom_vline(
        xintercept = 0, colour="black", linetype = "longdash", size = 0.1
    )
    if (split_by_granges) p <- p + ggplot2::facet_wrap(~ grange)
    # Return plot
    return(p)
}

#' A function to plot aggregated signals over sets of GRanges
#'
#' This function takes one or several RleList genomic tracks (imported
#' by rtraklayer::import(..., as = 'Rle')) and one or several GRanges 
#' objects. It computes coverage of the GRanges by the genomic tracks
#' and returns an aggregate coverage plot. 
#'
#' @param x several signal tracks (SimpleRleList or CompressedRleList class) 
#' grouped in a named list
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
#' @param split_by_granges Boolean Facet plots by the sets of GRanges
#' @param split_by_track Boolean Facet plots by the sets of signal tracks
#' @param free_scales Boolean Should each facet have independent y-axis scales?
#' @param run_in_parallel Boolean Should the plots be computed in parallel 
#' using mclapply?
#' @param norm character Should the signals be normalized 
#' ('none', 'zscore' or 'log2')?
#' @param ... additional parameters
#' @return A plot of aggregated signals
#' 
#' @import GenomicRanges
#' @import ggplot2
#' @importFrom zoo rollmean
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom stats qt
#' @export
#' 
#' @examples
#' library(GenomicRanges)
#' data(ce11_ATACseq)
#' data(ce11_WW_10bp)
#' data(ce11_proms)
#' proms <- resize(ce11_proms[1:100], fix = 'center', width = 400)
#' p1 <- plotAggregateCoverage(
#'     list(
#'         'atac' = ce11_ATACseq, 
#'         'WW_10bp' = ce11_WW_10bp
#'     ), 
#'     proms,
#'     norm = 'zscore'
#' )
#' p1
#' p2 <- plotAggregateCoverage(
#'     list(
#'         'ATAC-seq' = ce11_ATACseq, 
#'         'WW 10-bp periodicity' = ce11_WW_10bp
#'     ), 
#'     list(
#'         'Ubiq & Germline promoters' = 
#'             proms[proms$which.tissues %in% c('Ubiq.', 'Germline')],
#'         'Other promoters' = 
#'             proms[!(proms$which.tissues %in% c('Ubiq.', 'Germline'))]
#'     ), 
#'     norm = 'zscore'
#' )
#' p2
#' p3 <- plotAggregateCoverage(
#'     list(
#'         'ATAC-seq' = ce11_ATACseq, 
#'         'WW 10-bp periodicity' = ce11_WW_10bp
#'     ), 
#'     list(
#'         'Ubiq & Germline promoters' = 
#'             proms[proms$which.tissues %in% c('Ubiq.', 'Germline')],
#'         'Other promoters' = 
#'             proms[!(proms$which.tissues %in% c('Ubiq.', 'Germline'))]
#' ), 
#' split_by_granges = FALSE,
#' split_by_track = TRUE,
#' norm = 'zscore'
#' )
#' p3

plotAggregateCoverage.list <- function(
    x, 
    granges, 
    colors = NULL, 
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
    free_scales = FALSE,
    run_in_parallel = FALSE,
    norm = 'none',
    ...
) 
{
    bw_list <- x
    
    if (is.null(colors)) {
        colors <- rep(c(
            '#991919', 
            '#1232D9', 
            '#3B9B46', 
            '#D99B12', 
            '#7e7e7e', 
            '#D912D4', 
            '#9E7FCC', 
            '#B0E797', 
            '#D1B3B3', 
            '#23A4A3', 
            '#000000'
        ), 5)
    }
    if (!all(unlist(
        lapply(
            bw_list, 
            function(L) {
                class(L) %in% c('SimpleRleList', 'CompressedRleList')
            }
        )
    ))) {
        stop('Some objects in the bw_list are not bigwig tracks 
        (in SimpleRleList or CompressedRleList format). Aborting.')
    }
    if (methods::is(granges, 'GRanges')) granges <- list('granges' = granges)
    llists <- parallel::mclapply(seq_along(bw_list), function(B) {
        bw <- bw_list[[B]]
        # Compute scores
        lists <- parallel::mclapply(seq_along(granges), function(K) {
            g <- granges[[K]]
            if (length(unique(GenomicRanges::width(g))) > 1) {
                return(stop('Please provide GRanges that are 
                all the same width. Aborting.'))
            }
            mat <- getCovMatrix(g, bw, norm = norm, verbose = FALSE)
            colnames(mat) <- c(
                -unique(GenomicRanges::width(g))/2, 
                rep('', (ncol(mat) - 3)/2), 
                '1', 
                rep('', (ncol(mat) - 2)/2), 
                paste0('+', unique(GenomicRanges::width(g))/2)
            )
            row.names(mat) <- paste0(
                'locus_', as.character(seq_len(nrow(mat)))
            )
            means <- c(
                rep(NA, bin/2), 
                zoo::rollmean(apply(mat, 2, mean), bin), rep(NA, bin/2)
            )[seq_len(ncol(mat))]
            conint <- c(
                rep(NA, bin/2), 
                zoo::rollmean(
                    apply(
                        mat, 
                        2, 
                        function (n) { 
                            qt <- stats::qt(0.975,sum(!is.na(n)))
                            sd <- sd(n,na.rm=TRUE)
                            sqrt <- sqrt(sum(!is.na(n)))
                            return(qt*sd/sqrt)
                        }
                    ), 
                    bin
                ), 
                rep(NA, bin/2)
            )[seq_len(ncol(mat))]
            topEE <- means + conint
            bottomEE <- means - conint
            coords <- round(
                seq(
                    -unique(GenomicRanges::width(g))/2, 
                    unique(GenomicRanges::width(g))/2, 
                    length.out = unique(width(g))
                )
            )
            EE <- data.frame(
                'means' = means, 
                'meansUp' = topEE, 
                'meansDown' = bottomEE, 
                'x' = coords,
                'grange' = rep(
                    ifelse(
                        !is.null(names(granges)), 
                        names(granges)[K], 
                        as.character(K)
                    ), length(means)
                ), 
                'bw' = rep(
                    ifelse(
                        !is.null(names(bw_list)), 
                        names(bw_list)[B], 
                        as.character(K)
                    ), 
                    length(means)
                ), 
                stringsAsFactors = FALSE
            )
            return(EE)
        }) 
    }, mc.cores = ifelse(run_in_parallel, length(granges), 1)) 
    EEs <- do.call(rbind, do.call(rbind, llists))
    if (!is.null(names(granges))) {
        EEs$grange <- factor(EEs$grange, levels = names(granges))
    } 
    else {
        EEs$grange <- factor(
            EEs$grange, levels = as.character(seq_along(granges))
        )
    }
    if (!is.null(names(bw_list))) {
        EEs$bw <- factor(EEs$bw, levels = names(bw_list))
    } 
    else {
        EEs$bw <- factor(EEs$bw, levels = as.character(seq_along(bw_list)))
    }
    # Determining YLIM 
    if (is.null(ylim)) {
        max <- max(na.remove(EEs$meansUp)) * 1.05
        min <- min(na.remove(EEs$meansDown)) * 0.95
        ylim <- c(min, max)
    }
    # Start plotting
    if (split_by_granges & split_by_track) {
        p <- ggplot2::ggplot(data = EEs, ggplot2::aes(
            x = x, 
            y = means, 
            ymin = meansDown, 
            ymax = meansUp, 
            group = grange, 
            col = grange, 
            fill = grange
        ))
        p <- p + ggplot2::labs(x = xlab, y = ylab, col = 'Sets:')
        if (length(granges) > 1) {
            if (free_scales) {
                p <- p + ggplot2::facet_grid(bw~grange, scales = 'free')
            } 
            else {
                p <- p + ggplot2::facet_grid(bw~grange)
            }
        }
    } 
    else if (!split_by_granges & split_by_track) {
        p <- ggplot2::ggplot(data = EEs, ggplot2::aes(
            x = x, 
            y = means, 
            ymin = meansDown, 
            ymax = meansUp, 
            group = grange, 
            col = grange, 
            fill = grange
        ))
        p <- p + ggplot2::labs(x = xlab, y = ylab, col = 'Sets:')
        if (length(granges) > 1) {
            if (free_scales) {
                        p <- p + ggplot2::facet_grid(~bw, scales = 'free')
                    } 
            else {
                p <- p + ggplot2::facet_grid(~bw)
            }
        }
    } 
    else if (split_by_granges & !split_by_track) {
        p <- ggplot2::ggplot(data = EEs, ggplot2::aes(
            x = x, 
            y = means, 
            ymin = meansDown, 
            ymax = meansUp, 
            group = bw, 
            col = bw, 
            fill = bw
        ))
        p <- p + ggplot2::labs(x = xlab, y = ylab, col = 'Signal:')
        if (length(granges) > 1) {
            if (free_scales) {
                p <- p + ggplot2::facet_grid(~grange, scales = 'free')
            } 
            else {
                p <- p + ggplot2::facet_grid(~grange)
            }
        }
    }
    p <- p + ggplot2::geom_line()
    p <- p + ggplot2::geom_ribbon(alpha = 0.2, col = NA)
    p <- p + theme_ggplot2()
    p <- p + ggplot2::scale_x_continuous(
        limits = xlim, 
        expand = c(0, 0)
    )
    if (!free_scales) {
        p <- p + ggplot2::scale_y_continuous(
            limits = ylim, 
            expand = c(0, 0)
        )
    }
    p <- p + ggplot2::guides(fill = FALSE)
    if (length(granges) == 1) p <- p + ggplot2::guides(color = FALSE)
    if (length(granges) <= length(colors)) {
        p <- p + ggplot2::scale_fill_manual(values = colors)
        p <- p + ggplot2::scale_color_manual(values = colors)
    }
    if (plot_central) p <- p + ggplot2::geom_vline(
        xintercept = 0, colour="black", linetype = "longdash", size = 0.1
    )
    p <- p + ggplot2::theme(panel.spacing = unit(2, "lines"))
    # Return plot
    return(p)
}

###### ------- EXPERIMENTAL - DO NOT RUN ------- ######
# 
# 
# #' Core function - NOT WORKING YET
# #'
# #' @import ggplot2
# #' @import magrittr
# #' @importFrom zoo rollmean
# #' @importFrom reshape2 melt
# #' @export
# 
# plotHeatmaps <- function(
#     MAT, 
#     colorspace = NULL, 
#     autoscale = TRUE, 
#     s = 0.05, 
#     zmin = 0, 
#     zmax = 10, 
#     titles=rep('', length(MAT)),
#     xlab='',
#     ylab='', 
#     main='', 
#     reorder_decreasing = TRUE,
#     bin = 1,
#     clusters = 0,
#     #
#     cex.lab=12.0, 
#     cex.axis=12.0, 
#     cex.legend=12.0, 
#     cex.title=20,
#     pointsize=12,
#     ...
# )
# {
#     if (is(MAT, 'matrix') MAT <- list(MAT)
#     # Define color space
#     if (length(colorspace)) {
#         gcol <- colorRampPalette(colorspace)
#     } else {
#         gcol <- colorRampPalette(c(
#             "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
#             "#FF7F00", "red", "#7F0000"
#         ))     
#     }
#     NP <- length(MAT)
#     datapoints <- unlist(MAT)
#     plots <- list()
#     # Process each matrix 
#     for (i in seq(NP)) {
#         data <- MAT[[i]]
#         bins <- seq_len(ncol(data))
#         colnames(data) <- seq_len(ncol(data) )
#         # Bin data
#         if (bin > 1) {
#             data <- apply(data, 1, function(ROW) {
#                 zoo::rollmean(ROW, bin)
#             }) %>% t()
#             bins <- seq_len(ncol(data)) + bin/2 - 1
#             colnames(data) <- bins 
#         }
#         # Rescale outliers
#         if (autoscale) {
#             zlim <- quantile(data, c(s,1-s), na.rm=TRUE)
#             zmin <- zlim[1]
#             zmax <- zlim[2]
#         } 
#         data[data < zmin] <- zmin
#         data[data > zmax] <- zmax
#         # Cluster data
#         if (clusters > 0) {
#             data <- data[order(kmeans(data, centers = clusters)$cluster),]
#         } 
#         # Reorder data 
#         else if (reorder_decreasing) {
#             data <- data[order(rowSums(data), decreasing = TRUE),]
#         }
#         # Define color scale limits
#         keycolor_lim <- range(data, na.rm = TRUE)
#         keycolor_lim[1] <- zmin 
#         keycolor_lim[2] <- zmax 
#         # Plot graph
#         df <- reshape2::melt(data)
#         rect_df <- data.frame(
#             xmin = min(df$Var2), 
#             xmax = max(df$Var2), 
#             ymin = min(df$Var1), 
#             ymax = max(df$Var1)
#         )
#         p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = value))
#         p <- p + ggplot2::geom_raster()
#         p <- p + ggplot2::scale_fill_gradientn(
#             colours = gcol(50), limits = keycolor_lim, 
#             breaks = keycolor_lim, labels = format(keycolor_lim, digits = 2)
#         )
#         p <- p + ggplot2::scale_x_continuous(
#             breaks = c(min(bins), max(bins)),
#             labels = paste0(c('-', ''), (length(bins) / 2)),
#             expand = c(0.05, 0.05)
#         ) 
#         p <- p + ggplot2::geom_vline(
#             xintercept = length(bins)/2, 
#             size = .5, colour = 'black', linetype = 'dashed'
#         ) 
#         p <- p + ggplot2::scale_y_reverse(
#             breaks=c(1, nrow(data)),
#             labels=c(1, nrow(data)),
#             expand = c(0, 0)
#         ) 
#         p <- p + ggplot2::ggtitle(titles[i]) 
#         p <- p + ggplot2::xlab(xlab) 
#         p <- p + ggplot2::ylab(ylab) 
#         p <- p + ggplot2::theme(
#             legend.position = "bottom", 
#             axis.text=ggplot2::element_text(size=cex.axis), 
#             axis.title=ggplot2::element_text(size=cex.lab),
#             title=ggplot2::element_text(size=cex.lab),
#             legend.text=ggplot2::element_text(size=cex.legend) 
#         ) 
#         p <- p + ggplot2::guides(
#             fill = ggplot2::guide_colorbar(
#                 frame.colour = 'black', 
#                 frame.size = 0.25, title = "", raster = TRUE
#             )
#         )
#         p <- p + ggplot2::theme(
#             panel.background = ggplot2::element_blank(), 
#             panel.border = ggplot2::element_blank(), 
#             panel.grid.major = ggplot2::element_blank(),
#             panel.grid.minor = ggplot2::element_blank(), 
#             axis.line = ggplot2::element_blank(),
#         ) 
#         p <- p + ggplot2::geom_rect(
#             data = rect_df,
#             mapping = ggplot2::aes(
#                 xmin = xmin,
#                 xmax = xmax,
#                 ymin = ymin,
#                 ymax = ymax
#             ), 
#             inherit.aes = FALSE,
#             fill = NA, 
#             color = 'black', 
#             size = 0.5
#         )
#         plots[[i]] <- p
#     }
#     if (length(plots) == 1) plots <- plots[[1]]
#     return(plots)
# }
# 
# #' Core function - NOT WORKING YET
# #' 
# #' @import magrittr
# #' @import Biostrings
# #' @import GenomicRanges
# #' @export
# 
# getMotifMatrix <- function(
#     granges, 
#     motif, 
#     seqs = Biostrings::readDNAStringSet(
#         "~/shared/sequences/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
#     ), 
#     N_MISMATCHES = 0, 
#     ...
# ) 
# {
#     pattern <- Biostrings::DNAString(motif)
#     seqs_g <- seqs[granges]
#    if (!all(sapply(
#        unlist(strsplit(as.character(pattern), '')), 
#        function(CHAR) {CHAR %in% c('A', 'T', 'C', 'G')}
#    ))) {
#         ALGO <- "shift-or"
#         FIXED <- FALSE
#     } else {
#         if (N_MISMATCHES == 0) {
#             ALGO <- 'naive-exact'
#             FIXED <- TRUE
#         } else {
#             ALGO <- 'auto'
#             FIXED <- TRUE
#         }
#     }
#     hits <- Biostrings::vmatchPattern(
#         pattern, seqs_g, algorithm = ALGO, fixed = FIXED, 
#         max.mismatch = N_MISMATCHES
#     ) %>% 
#         as("CompressedIRangesList") %>%
#         IRanges::resize(1, fix = 'center')
#     rhits <- Biostrings::vmatchPattern(
#         pattern, seqs_g, algorithm = ALGO, fixed = FIXED,
#         max.mismatch = N_MISMATCHES
#     ) %>% 
#         as("CompressedIRangesList") %>%
#         resize(1, fix = 'center')
#     out <- lapply(seq_along(granges), function(K) {
#         if (as.vector(GenomicRanges::strand(granges))[K] == '+') {
#             hits[[K]]
#         } else {
#             rhits[[K]]
#         }
#     }) %>% 
#         as("CompressedIRangesList") %>%
#         GenomicRanges::coverage(width = GenomicRanges::width(seqs_g)) %>% 
#         sapply(., function(x) {
#             as.numeric(x)
#         }) %>% 
#         t()
#     return(out)
# }
# 
