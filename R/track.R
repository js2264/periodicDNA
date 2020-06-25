#' Function to generate a k-mer periodicity track
#' 
#' This function takes a set of GRanges in a genome, recover 
#' the corresponding sequences and divides them using a sliding window. 
#' For each sub-sequence, it then computes the PSD value of a k-mer 
#' of interest at a chosen period, and generates a linear .bigWig
#' track from these values. 
#'
#' @param genome DNAStringSet, BSgenome or genome ID
#' @param granges GRanges object
#' @param motif character, k-mer of interest. 
#' @param extension Integer, the width the GRanges are going to 
#' be extended to (default 1000).
#' @param window_size Integer, the width of the bins to split 
#' the GRanges objects in (default 100).
#' @param step_size Integer, the increment between bins 
#' over GRanges (default 2).
#' @param range_spectrum Numeric vector, the distances between nucleotides
#' to take into consideration when performing Fast Fourier Transform 
#' (default seq_len(50)).
#' @param period Integer, the period of the dinucleotide to study 
#' (default=10).
#' @param smooth_track Integer, smooth the resulting track
#' @param BPPARAM split the workload over several processors using 
#' BiocParallel
#' @param bw_file character, the name of the output bigWig track
#' @return Rlelist and a bigWig track in the working directory. 
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import BiocParallel
#' @importFrom rtracklayer export.bw
#' @importFrom methods is
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' getPeriodicityTrack(
#'     genome = 'ce11', 
#'     ce11_proms[1], 
#'     motif = 'WW',
#'     period = 10,
#'     BPPARAM = setUpBPPARAM(1)
#' )
#' unlink('WW_10-bp-periodicity_g-100^2_smooth-20.bw')

getPeriodicityTrack <- function(
    genome = NULL,
    granges,
    motif = 'WW',
    period = 10,
    BPPARAM = bpparam(),
    extension = 1000,
    window_size = 100,
    step_size = 2,
    range_spectrum = seq(5, 50),
    smooth_track = 20,
    bw_file = NULL
)
{
    cores <- BPPARAM$workers
    freq <- 1/period
    #
    if (is.null(genome)) {
        genome <- 'ce11'
    }
    if (methods::is(genome, 'character')) {
        if (genome %in% c(
            'sacCer3', 'ce11', 'dm6', 'mm10', 'hg38', 'danRer10'
        )) {
            genome <- char2BSgenome(genome)
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
    #
    isCircular(seqinfo(granges)) <- NA
    granges <- reduce(GRanges(granges, strand = '*'))
    granges_extended <- GenomicRanges::resize(
        granges, extension, fix = 'center'
    )
    granges_extended_large <- GenomicRanges::resize(
        granges, 2*extension, fix = 'center'
    )
    genome_partionned <- partitionGenome(
        genome, granges_extended_large, granges_extended, 
        window_size, step_size
    )
    genome_partionned_filtered <- IRanges::subsetByOverlaps(
        genome_partionned, granges_extended
    )
    chunks <- getChunks(genome_partionned_filtered, cores = cores)
    #
    message('
    The genome has been divided in ', window_size, '-bp long windows 
    with a sliding window of ', step_size, '-bp.
    ', length(genome_partionned_filtered), ' windows overlap the ', 
    length(granges), ' input loci (extended to ', extension, '-bp).
    The mapping will be split into ', cores, ' cores. Each core will 
    process ', lengths(chunks)[1], ' windows.
    
    Now starting [', Sys.time(), ']...
    ')
    if (is.null(bw_file)) bw_file <- paste0(
        motif, '_', period, '-bp-periodicity',
        '_g-', window_size, '^', step_size, 
        '_smooth-', smooth_track,
        '.bw'
    )
    
    # Prepare log file
    logfile <- paste0('log.', bw_file, '.txt')
    writeLines('>> Progress:', logfile)
    sink(logfile, append = TRUE)
    for (K in seq_along(chunks)) {
        procID <- K
        cat(sprintf(
            '- Proc. %s:  %s windows processed (%s remaining)\n', 
            procID, 0, length(chunks[[K]])
        ), file = logfile, append = TRUE)
    }
    sink()
    
    # Process each chunk separately
    list_results <- BiocParallel::bplapply(
        chunks, 
        chunkWrapper,
        chunks,
        genome, 
        motif, 
        step_size,
        range_spectrum, 
        freq, 
        logfile,
        BPPARAM = BPPARAM
    )
    
    # Merge all the chunks
    message('    Merging the results...
    ')
    list_results_2 <- unlist(GRangesList(list_results))
    
    # Generate final track
    res <- coverage(list_results_2, weight = list_results_2$score)
    if (smooth_track > 1) {
        message('    Smoothing the track...
        ')
        res <- smoothBigWig(res, k = smooth_track, BPPARAM)
    }
    rtracklayer::export.bw(
        res, 
        bw_file
    )
    message('\n   SUCCESS: ', bw_file, ' has been created!')
    message(
        '   ', length(list_results_2)*step_size, 
        ' bases covered by the generated track.'
    )
    # Remove tmp files
    cleanUpDirectory(logfile)
    
    message('
    Finished without errors [', Sys.time(), ']
    ')
    
    # Return track
    return(res)
}

#' Internal function
#'
#' @param genome DNAStringSet, BSgenome or genome ID
#' @param granges_extended_large A GRanges object.
#' @param granges_extended A GRanges object.
#' @param window_size An integer. The width of the bins to split 
#' the GRanges objects in (default 100).
#' @param step_size An integer. The increment between bins 
#' over GRanges (default 2).
#' @return granges_partionned A Granges object corresponding to the binned 
#' genome. 
#' 
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges

partitionGenome <- function(
    genome, 
    granges_extended_large, 
    granges_extended, 
    window_size, 
    step_size
) 
{
    genome_Seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = names(genome), 
        seqlengths = lengths(genome), 
        isCircular = rep(FALSE, length(genome))
    )
    genome_granges <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(genome_Seqinfo), 
        IRanges::IRanges(start = rep(1, length(genome_Seqinfo)), 
        width = GenomeInfoDb::seqlengths(genome_Seqinfo))
    )
    GenomeInfoDb::seqinfo(genome_granges) <- genome_Seqinfo
    granges_partionned <- GenomicRanges::slidingWindows(
        GenomicRanges::reduce(granges_extended_large), 
        window_size, 
        step_size
    ) %>% 
        GenomicRanges::GRangesList() %>% 
        unlist()
    GenomeInfoDb::seqinfo(granges_partionned) <- genome_Seqinfo
    granges_partionned <- GenomicRanges::trim(granges_partionned)
    granges_partionned <- granges_partionned[
        GenomicRanges::width(granges_partionned) == window_size
    ]
    return(granges_partionned)
}

#' Internal function
#'
#' @param genome_partionned_filtered Output from partitionGenome() filtered 
#' over GRanges of interest. 
#' @param cores An integer Split the workload over several processors 
#' (default 12).
#' @return list Chunks list of primary bins to measure.

getChunks <- function(genome_partionned_filtered, cores) {
    chunks <- as.numeric(
        cut(
            seq_len(length(genome_partionned_filtered)), 
            breaks = seq(
                1, length(genome_partionned_filtered), length.out = cores + 1
            ), 
            include.lowest = TRUE)
        )
    list_granges <- list()
    for (K in seq_len(cores)) {
        list_granges[[K]] <- genome_partionned_filtered[chunks == K]
    }
    return(list_granges)
}

#' Internal function
#'
#' @param chunk A GRanges from the chunks list obtained with getChunks
#' @param chunks GRanges split in a list with getChunks
#' @param genome DNAStringSet, BSgenome or genome ID
#' @param motif String Oligonucleotide of interest. 
#' @param step_size Integer The increment between 
#' bins over GRanges (default 2).
#' @param range_spectrum Numeric vector The distances between nucleotides
#' to take into consideration when performing Fast Fourier Transform 
#' (default seq_len(50)).
#' @param freq Float The frequence of the dinucleotide to study 
#' (default 1/10).
#' @param logfile character, name of log file 
#' @return GRanges a GRanges object with a score column containing 
#' estimated periodicity
#' 
#' @import GenomicRanges
#' @importFrom rtracklayer export.bw

chunkWrapper <- function(
    chunk,
    chunks, 
    genome, 
    motif, 
    step_size, 
    range_spectrum, 
    freq, 
    logfile
)
{
    res <- chunk
    res$score <- 0
    res <- GenomicRanges::resize(
        res, 
        width = step_size, 
        fix = 'center'
    )
    checks <- unique(round(seq(1, length(chunk), length.out = 50)))
    procID <- which(lengths(which(chunks == chunk)) > 0)
    for (K in seq_len(length(chunk))) {
        if (K %in% checks) {
            pct <- which(checks == K)
            progbar <- paste0(
                '|', 
                paste(rep('-', pct), collapse = ''), 
                paste(rep('.', 50-pct), collapse = ''), 
                '|'
            )
            txt <- sprintf(
                'Latest processor (#%s): %s %s%% [%s processed, %s remaining]', 
                procID, progbar, round(K/length(chunk)*100, 0), 
                K, length(chunk) - K
            )
            writeLines(txt, logfile)
        }
        res$score[K] <- binWrapper(
            chunk[K], genome, motif, range_spectrum, freq
        )
    }
    return(res)
}

#' Internal function
#'
#' @param grange A GRanges object of length 1
#' @param genome DNAStringSet, BSgenome or genome ID
#' @param motif String Oligonucleotide of interest. 
#' @param range_spectrum Numeric vector The distances between nucleotides
#' to take into consideration when performing Fast Fourier Transform 
#' (default seq_len(50)).
#' @param freq Float The frequence of the dinucleotide to study
#' (default 1/10).
#' @return Float The estimated periodicity score of the input GRanges
#' 
#' @import magrittr
#' @importFrom stats dist
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @importFrom stats spectrum
#' @importFrom zoo rollmean

binWrapper <- function(
    grange, 
    genome, 
    motif,
    range_spectrum, 
    freq
) 
{
    seqs <- withSeq(grange, genome)$seq
    dists <- unlist(lapply(seq_along(seqs), function(k) {
        Biostrings::vmatchPattern(
            motif, seqs[k], max.mismatch = 0, fixed = FALSE
        ) %>% 
        IRanges::start() %>% 
        '[['(1) %>% 
        dist() %>% 
        c()
    }))
    dists <- dists[dists <= max(range_spectrum)]
    hist <- hist(dists, breaks = seq(0, max(range_spectrum), 1), plot = FALSE)
    hist <- hist$counts[range_spectrum]
    hist <- zoo::rollmean(
        hist, 
        k = 3, fill = 0, na.pad = TRUE, align = 'center'
    )
    s <- stats::spectrum(hist, plot = FALSE)
    spec <- round(s$spec[which.min(abs(s$freq - freq))], 4)
    return(spec)
}

#' Internal function
#'
#' @param bw A RleList object
#' @param k Integer, smoothing window
#' @param BPPARAM BPPARAM
#' @return bw RleList smoothed (for each chromosome) using a window of k.
#' 
#' @import BiocParallel
#' @import S4Vectors
#' @importFrom methods as
#' @importFrom zoo rollmean

smoothBigWig <- function(bw, k = 20, BPPARAM) {
    BPPARAM$workers <- min(length(bw), BPPARAM$workers)
    v <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        names(bw), 
        function(K) {
            S4Vectors::Rle(zoo::rollmean(
                as.vector(bw[[K]]), k = k, fill = 0, align = 'center'
            ))
        }
    )
    names(v) <- names(bw)
    v <- methods::as(v, "RleList")
    return(v)
}

#' Internal function
#' 
#' @param logfile character, name of log file 
#' @return NULL Clean-up temporary files generated when running 
#' getPeriodicityTrack.
#' 
#' @import magrittr

cleanUpDirectory <- function(logfile) {
    list.files(pattern = '.*tmp.*') %>% file.remove()
    file.remove(logfile)
    return(TRUE)
}