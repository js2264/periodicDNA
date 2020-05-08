#' Core function to generate a track of periodicity strenght of a given 
#' dinucleotide at a given frequency, over a set of chosen GRanges.
#'
#' @param genome DNAStringSet, BSgenome or genome ID
#' @param granges GRanges object (with seqnames overlapping 
#'     the names of genome).
#' @param motif String Oligonucleotide of interest. 
#' @param extension Integer The width the GRanges are going to 
#'     be extended to (default 1000).
#' @param genome_sliding_size Integer The width of the bins to split the GRanges
#'     objects in (default 100).
#' @param genome_sliding_sliding Integer The increment between bins over GRanges
#'     (default 2).
#' @param window_sliding_size Integer The width of the bins to split each primary 
#'     bin in (default 60).
#' @param window_sliding_bin Integer The increment between secondary bins
#'     (default 5).
#' @param range_spectrum Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform 
#'     (default seq_len(50)).
#' @param period Integer The period of the dinucleotide to study (default=10).
#' @param cores Integer Split the workload over several processors 
#'     (default 12).
#' @param bw_file String. The name of the output bigWig track
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @importFrom parallel mclapply
#' @importFrom rtracklayer export.bw
#' @importFrom methods is
#' @export
#' @return NULL A bigWig track in the working directory. 

getPeriodicityTrack <- function(
    genome = NULL, 
    granges, 
    motif = 'WW', 
    extension = 1000, 
    genome_sliding_size = 100, 
    genome_sliding_sliding = 2, 
    window_sliding_size = 100, 
    window_sliding_bin = 5, 
    range_spectrum = seq(5, 50), 
    period = 10, 
    cores = 12, 
    bw_file = NULL
)
{
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
    granges.extended <- GenomicRanges::resize(
        granges, extension, fix = 'center'
    )
    granges.extended_large <- GenomicRanges::resize(
        granges, 2*extension, fix = 'center'
    )
    genome.partionned <- partitionGenome(
        genome, granges.extended_large, granges.extended, 
        genome_sliding_size, genome_sliding_sliding
    )
    genome.partionned.filtered <- IRanges::subsetByOverlaps(
        genome.partionned, granges.extended
    )
    chunks <- getChunks(genome.partionned.filtered, cores = cores)
    
    message('
    The genome has been divided in ', genome_sliding_size, '-bp bins
    with a sliding window of ', genome_sliding_sliding, '-bp.
    
    ', length(genome.partionned.filtered), ' windows overlap the input 
    loci (extended to ', extension, '-bp).
    
    The mapping will be split among ', cores, ' jobs.
    
    Each job will analyse ', lengths(chunks)[1], ' bins.
    
    Now starting...
    ')
    
    # Process each chunk separately
    list.results <- parallel::mclapply(
        chunks, 
        chunkWrapper,
        genome, 
        motif, 
        genome_sliding_sliding,
        window_sliding_size, 
        window_sliding_bin,
        range_spectrum, 
        freq,
        mc.cores = length(chunks)
    )
    
    # Merge all the chunks
    list.results.2 <- unlist(
        GRangesList(lapply(list.files(pattern = 'FULL'), rtracklayer::import))
    )
    
    # Generate final track
    if (is.null(bw_file)) bw_file = paste0(
        motif, 
        '-periodicity_g-', genome_sliding_size, '^', genome_sliding_sliding, 
        '_b-', window_sliding_size, '^', window_sliding_bin ,'.bw'
    )
    res <- coverage(list.results.2, weight = list.results.2$score)
    rtracklayer::export.bw(
        res, 
        bw_file
    )
    message('\n   SUCCESS: ', bw_file, ' has been created!')
    message(
        '   ', sum(sum(res != 0)), ' bases covered by the generated track.'
    )
    # Remove tmp files
    cleanUpDirectory()
    
    # Return track
    return(bw_file)
}

#' Internal function
#'
#' @param genome A DNAStringSet object. Ideally, the sequence of an 
#'     entire genome, obtained for instance by running 
#'     \code{Biostrings::getSeq(
#'       BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
#'     )}.
#' @param granges.extended_large A GRanges object.
#' @param granges.extended A GRanges object.
#'     \code{Biostrings::getSeq(
#'       BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
#'     )}.
#' @param genome_sliding_size An integer. The width of the bins to split 
#'     the GRanges objects in (default 100).
#' @param genome_sliding_sliding An integer. The increment between bins 
#'     over GRanges (default 2).
#' 
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @return granges.partionned A Granges object corresponding to the binned 
#'   genome. 

partitionGenome <- function(
    genome, 
    granges.extended_large, 
    granges.extended, 
    genome_sliding_size, 
    genome_sliding_sliding
) 
{
    genome.Seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = names(genome), 
        seqlengths = lengths(genome), 
        isCircular = rep(FALSE, length(genome))
    )
    genome.granges <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(genome.Seqinfo), 
        IRanges::IRanges(start = rep(1, length(genome.Seqinfo)), 
        width = GenomeInfoDb::seqlengths(genome.Seqinfo))
    )
    GenomeInfoDb::seqinfo(genome.granges) <- genome.Seqinfo
    granges.partionned <- GenomicRanges::slidingWindows(
        GenomicRanges::reduce(granges.extended_large), 
        genome_sliding_size, 
        genome_sliding_sliding
    ) %>% 
        GenomicRanges::GRangesList() %>% 
        unlist()
    GenomeInfoDb::seqinfo(granges.partionned) <- genome.Seqinfo
    granges.partionned <- GenomicRanges::trim(granges.partionned)
    granges.partionned <- granges.partionned[
        GenomicRanges::width(granges.partionned) == genome_sliding_size
    ]
    return(granges.partionned)
}

#' Internal function
#'
#' @param genome.partionned.filtered Output from partitionGenome() filtered 
#'     over GRanges of interest. 
#' @param cores An integer Split the workload over several processors 
#'     (default 12).
#' 
#' @return list Chunks list of primary bins to measure.

getChunks <- function(genome.partionned.filtered, cores) {
    chunks <- as.numeric(
        cut(
            seq_len(length(genome.partionned.filtered)), 
            breaks = seq(
                1, length(genome.partionned.filtered), length.out = cores + 1
            ), 
            include.lowest = TRUE)
        )
    list.granges <- list()
    for (K in seq_len(cores)) {
        list.granges[[K]] <- genome.partionned.filtered[chunks == K]
    }
    return(list.granges)
}

#' Internal function
#'
#' @param chunk A GRanges from the output list obtained with getChunks
#' @param genome A DNAStringSet object. Ideally, the sequence of an 
#'     entire genome, obtained for instance by running 
#'     \code{Biostrings::getSeq(
#'       BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
#'     )}.
#' @param motif String Oligonucleotide of interest. 
#' @param genome_sliding_sliding Integer The increment between bins over GRanges
#'     (default 2).
#' @param window_sliding_size Integer The width of the bins to split each primary 
#'     bin in (default 60).
#' @param window_sliding_bin Integer The increment between secondary bins
#'     (default 5).
#' @param range_spectrum Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform 
#'     (default seq_len(50)).
#' @param freq Float The frequence of the dinucleotide to study (default 1/10).
#' 
#' @import GenomicRanges
#' @importFrom rtracklayer export.bw
#' @return GRanges a GRanges object with a score column containing estimated
#'     periodicity

chunkWrapper <- function(
    chunk, 
    genome, 
    motif, 
    genome_sliding_sliding, 
    window_sliding_size, 
    window_sliding_bin, 
    range_spectrum, 
    freq
) 
{
    res <- chunk
    res$score <- 0
    res <- GenomicRanges::resize(
        res, width = genome_sliding_sliding, fix = 'center'
    )
    intervals <- seq(1, length(res), length.out = 100)
    for (K in seq_len(length(chunk))) {
        res$score[K] <- binWrapper(
            chunk[K], genome, motif, window_sliding_size, 
            window_sliding_bin, range_spectrum, freq
        )
        if (K %in% intervals) {
            message('>> ', K, ' bins computed /// TIME: ', Sys.time())
            rtracklayer::export.bw(
                res, paste0("tmp.", K, ".", Sys.getpid(), ".bw")
            )
        }
    }
    rtracklayer::export.bw(res, paste0("tmp.FULL.", Sys.getpid(), ".bw"))
    return(res)
}

#' Internal function
#'
#' @param grange A GRanges object of length 1
#' @param genome A DNAStringSet object. Ideally, the sequence of an entire 
#'     genome
#' @param motif String Oligonucleotide of interest. 
#' @param window_sliding_size Integer The width of the bins to split each primary
#'     bin in (default 60).
#' @param window_sliding_bin Integer The increment between secondary bins
#'     (default 5).
#' @param range_spectrum Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform 
#'     (default seq_len(50)).
#' @param freq Float The frequence of the dinucleotide to study (default 1/10).
#' 
#' @import magrittr
#' @importFrom stats dist
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @importFrom stats spectrum
#' @return Float The estimated periodicity score of the input GRanges

binWrapper <- function(
    grange, 
    genome, 
    motif,
    window_sliding_size,
    window_sliding_bin, 
    range_spectrum, 
    freq
) 
{
    bins <- withSeq(
        partitionBin(grange, genome, window_sliding_size, window_sliding_bin), 
        genome
    )
    seqs <- bins$seq
    dists <- unlist(lapply(seq_along(seqs), function(k) {
        Biostrings::vmatchPattern(
            motif, seqs[k], max.mismatch = 0, fixed = FALSE
        ) %>% 
        IRanges::start() %>% 
        '[['(1) %>% 
        dist() %>% 
        c()
    }))
    hist <- hist(dists, breaks = seq(0, 1000, 1), plot = FALSE)
    s <- stats::spectrum(hist$counts[range_spectrum], plot = FALSE)
    spec <- round(s$spec[which.min(abs(s$freq - freq))], 4)
    return(spec)
}

#' Internal function
#'
#' @param grange A GRanges object of length 1
#' @param genome A DNAStringSet object. Ideally, the sequence of an entire 
#'     genome
#' @param window_sliding_size Integer The width of the bins to split each primary
#'     bin in (default 60).
#' @param window_sliding_bin Integer The increment between secondary bins
#'     (default 5).
#' 
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @return list A GRanges object

partitionBin <- function(grange, genome, window_sliding_size, window_sliding_bin) {
    genome.Seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = names(genome), 
        seqlengths = lengths(genome), 
        isCircular = rep(FALSE, length(genome)), 
        genome = 'custom'
    )
    genome.granges <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(genome.Seqinfo), 
        IRanges::IRanges(
            start = rep(1, length(genome.Seqinfo)), 
            width = GenomeInfoDb::seqlengths(genome.Seqinfo)
        )
    )
    GenomeInfoDb::seqinfo(genome.granges) <- genome.Seqinfo
    starts <- seq(
        IRanges::start(grange), 
        IRanges::start(grange)+GenomicRanges::width(grange)-window_sliding_size, 
        window_sliding_bin
    )
    granges <- GenomicRanges::GRanges(
        seqnames = rep(GenomicRanges::seqnames(grange), length(starts)),
        IRanges::IRanges(start = starts, width = window_sliding_size)
    )
    return(granges)
}

#' Internal function
#'
#' @param list.results List 
#' 
#' @return GRanges

unlistResults <- function(list.results) {
    granges <- list.results[[1]][0]
    for (K in seq_len(length((list.results)))) {
        granges <- c(granges, list.results[[K]])
    }
    return(granges)
}

#' Internal function
#' 
#' @import magrittr
#' @return NULL Clean-up temporary files generated when running 
#'     getPeriodicityTrack.

cleanUpDirectory <- function() {
    list.files(pattern = '.*tmp.*') %>% file.remove()
    return(TRUE)
}