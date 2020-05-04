#' Core function to generate a track of periodicity strenght of a given 
#' dinucleotide at a given frequency, over a set of chosen GRanges.
#'
#' @param genome DNAStringSet object obtained from
#'     \code{Biostrings::getSeq(BSgenome.xxx)}.
#' @param granges GRanges object (with seqnames overlapping 
#'     the names of genome).
#' @param MOTIF String Oligonucleotide of interest. 
#' @param EXTEND.GRANGES Integer The width the GRanges are going to 
#'     be extended to (default 1000).
#' @param GENOME.WINDOW.SIZE Integer The width of the bins to split the GRanges
#'     objects in (default 100).
#' @param GENOME.WINDOW.SLIDING Integer The increment between bins over GRanges
#'     (default 2).
#' @param BIN.WINDOW.SIZE Integer The width of the bins to split each primary 
#'     bin in (default 60).
#' @param BIN.WINDOW.SLIDING Integer The increment between secondary bins
#'     (default 5).
#' @param RANGE.FOR.SPECTRUM Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform 
#'     (default seq_len(50)).
#' @param FREQ Float The frequence of the dinucleotide to study (default 1/10).
#' @param PROCS Integer Split the workload over several processors 
#'     (default 12).
#' @param bw_file String. The name of the output bigWig track
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @importFrom parallel mclapply
#' @importFrom rtracklayer export.bw
#' @export
#' @return NULL A bigWig track in the working directory. 

generatePeriodicityTrack <- function(
    genome = NULL, 
    granges, 
    MOTIF = 'WW', 
    EXTEND.GRANGES = 1000, 
    GENOME.WINDOW.SIZE = 100, 
    GENOME.WINDOW.SLIDING = 2, 
    BIN.WINDOW.SIZE = 100, 
    BIN.WINDOW.SLIDING = 5, 
    RANGE.FOR.SPECTRUM = seq_len(50), 
    FREQ = 1/10, 
    PROCS = 12, 
    bw_file = NULL
    ) {
    
    if (is.null(genome)) {
        genome <- Biostrings::getSeq(
            BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
        )
    }
    isCircular(seqinfo(granges)) <- NA
    granges.extended <- GenomicRanges::resize(
        granges, EXTEND.GRANGES, fix = 'center'
    )
    granges.extended_large <- GenomicRanges::resize(
        granges, 2*EXTEND.GRANGES, fix = 'center'
    )
    genome.partionned <- partitionGenome(
        genome, granges.extended_large, granges.extended, 
        GENOME.WINDOW.SIZE, GENOME.WINDOW.SLIDING
    )
    genome.partionned.filtered <- IRanges::subsetByOverlaps(
        genome.partionned, granges.extended
    )
    chunks <- getChunks(genome.partionned.filtered, PROCS = PROCS)
    
    message('
    The genome has been divided in ', GENOME.WINDOW.SIZE, '-bp bins
    with a sliding window of ', GENOME.WINDOW.SLIDING, '-bp.
    
    ', length(genome.partionned.filtered), ' windows overlap the input 
    loci (extended to ', EXTEND.GRANGES, '-bp).
    
    The mapping will be split among ', PROCS, ' jobs.
    
    Each job will analyse ', lengths(chunks)[1], ' bins.
    
    Now starting...
    ')
    
    # Process each chunk separately
    list.results <- parallel::mclapply(
        chunks, 
        chunkWrapper,
        genome, 
        MOTIF, 
        GENOME.WINDOW.SLIDING,
        BIN.WINDOW.SIZE, 
        BIN.WINDOW.SLIDING,
        RANGE.FOR.SPECTRUM, 
        FREQ,
        mc.cores = length(chunks)
    )
    
    # Merge all the chunks
    list.results.2 <- unlist(
        GRangesList(lapply(list.files(pattern = 'FULL'), rtracklayer::import))
    )
    
    # Generate final track
    if (is.null(bw_file)) bw_file = paste0(
        MOTIF, 
        '-periodicity_g-', GENOME.WINDOW.SIZE, '^', GENOME.WINDOW.SLIDING, 
        '_b-', BIN.WINDOW.SIZE, '^', BIN.WINDOW.SLIDING ,'.bw'
    )
    rtracklayer::export.bw(
        coverage(list.results.2, weight = list.results.2$score), 
        bw_file
    )
    message(' SUCCESS: ', bw_file, ' has been created!')
    
    # Remove tmp files
    cleanUpDirectory()
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
#' @param GENOME.WINDOW.SIZE An integer. The width of the bins to split 
#'     the GRanges objects in (default 100).
#' @param GENOME.WINDOW.SLIDING An integer. The increment between bins 
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
    GENOME.WINDOW.SIZE, 
    GENOME.WINDOW.SLIDING
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
        GENOME.WINDOW.SIZE, 
        GENOME.WINDOW.SLIDING
    ) %>% 
        GenomicRanges::GRangesList() %>% 
        unlist()
    GenomeInfoDb::seqinfo(granges.partionned) <- genome.Seqinfo
    granges.partionned <- GenomicRanges::trim(granges.partionned)
    granges.partionned <- granges.partionned[
        GenomicRanges::width(granges.partionned) == GENOME.WINDOW.SIZE
    ]
    return(granges.partionned)
}

#' Internal function
#'
#' @param genome.partionned.filtered Output from partitionGenome() filtered 
#'     over GRanges of interest. 
#' @param PROCS An integer Split the workload over several processors 
#'     (default 12).
#' 
#' @return list Chunks list of primary bins to measure.

getChunks <- function(genome.partionned.filtered, PROCS) {
    chunks <- as.numeric(
        cut(
            seq_len(length(genome.partionned.filtered)), 
            breaks = seq(
                1, length(genome.partionned.filtered), length.out = PROCS + 1
            ), 
            include.lowest = TRUE)
        )
    list.granges <- list()
    for (K in seq_len(PROCS)) {
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
#' @param MOTIF String Oligonucleotide of interest. 
#' @param GENOME.WINDOW.SLIDING Integer The increment between bins over GRanges
#'     (default 2).
#' @param BIN.WINDOW.SIZE Integer The width of the bins to split each primary 
#'     bin in (default 60).
#' @param BIN.WINDOW.SLIDING Integer The increment between secondary bins
#'     (default 5).
#' @param RANGE.FOR.SPECTRUM Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform 
#'     (default seq_len(50)).
#' @param FREQ Float The frequence of the dinucleotide to study (default 1/10).
#' 
#' @import GenomicRanges
#' @importFrom rtracklayer export.bw
#' @return GRanges a GRanges object with a score column containing estimated
#'     periodicity

chunkWrapper <- function(
    chunk, 
    genome, 
    MOTIF, 
    GENOME.WINDOW.SLIDING, 
    BIN.WINDOW.SIZE, 
    BIN.WINDOW.SLIDING, 
    RANGE.FOR.SPECTRUM, 
    FREQ
) 
{
    res <- chunk
    res$score <- 0
    res <- GenomicRanges::resize(
        res, width = GENOME.WINDOW.SLIDING, fix = 'center'
    )
    intervals <- seq(1, length(res), length.out = 100)
    for (K in seq_len(length(chunk))) {
        res$score[K] <- binWrapper(
            chunk[K], genome, MOTIF, BIN.WINDOW.SIZE, 
            BIN.WINDOW.SLIDING, RANGE.FOR.SPECTRUM, FREQ
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
#' @param MOTIF String Oligonucleotide of interest. 
#' @param BIN.WINDOW.SIZE Integer The width of the bins to split each primary
#'     bin in (default 60).
#' @param BIN.WINDOW.SLIDING Integer The increment between secondary bins
#'     (default 5).
#' @param RANGE.FOR.SPECTRUM Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform 
#'     (default seq_len(50)).
#' @param FREQ Float The frequence of the dinucleotide to study (default 1/10).
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
    MOTIF,
    BIN.WINDOW.SIZE,
    BIN.WINDOW.SLIDING, 
    RANGE.FOR.SPECTRUM, 
    FREQ
) 
{
    bins <- withSeq(
        partitionBin(grange, genome, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING), 
        genome
    )
    seqs <- bins$seq
    dists <- unlist(lapply(seq_along(seqs), function(k) {
        Biostrings::vmatchPattern(
            MOTIF, seqs[k], max.mismatch = 0, fixed = FALSE
        ) %>% 
        IRanges::start() %>% 
        '[['(1) %>% 
        dist() %>% 
        c()
    }))
    hist <- hist(dists, breaks = seq(0, 1000, 1), plot = FALSE)
    s <- stats::spectrum(hist$counts[RANGE.FOR.SPECTRUM], plot = FALSE)
    spec <- round(s$spec[which.min(abs(s$freq - FREQ))], 4)
    return(spec)
}

#' Internal function
#'
#' @param grange A GRanges object of length 1
#' @param genome A DNAStringSet object. Ideally, the sequence of an entire 
#'     genome
#' @param BIN.WINDOW.SIZE Integer The width of the bins to split each primary
#'     bin in (default 60).
#' @param BIN.WINDOW.SLIDING Integer The increment between secondary bins
#'     (default 5).
#' 
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @return list A GRanges object

partitionBin <- function(grange, genome, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING) {
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
        IRanges::start(grange)+GenomicRanges::width(grange)-BIN.WINDOW.SIZE, 
        BIN.WINDOW.SLIDING
    )
    granges <- GenomicRanges::GRanges(
        seqnames = rep(GenomicRanges::seqnames(grange), length(starts)),
        IRanges::IRanges(start = starts, width = BIN.WINDOW.SIZE)
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
#'     generatePeriodicityTrack.

cleanUpDirectory <- function() {
    list.files(pattern = '.*tmp.*') %>% file.remove()
    return(TRUE)
}