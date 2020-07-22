namedListToLongFormat <- function(x) {
    l <- lapply(names(x), function(NAME) {
        L <- x[[NAME]]
        if (is.null(ncol(L))) {
            data.frame(value = L, name = rep(NAME, length(L)))
        } 
        else {
            data.frame(L, name = rep(NAME, nrow(L)))
        }
    }) 
    res <- do.call(rbind, l)
    return(res)
}

na.replace <- function(x, value) {
    which.na <- is.na(x)
    x[which.na] <- value
    return(x)
}

na.remove <- function(x) {
    which.na <- is.na(x)
    return(x[!which.na])
}

#' A function to shuffle sequence(s)
#' 
#' This function takes a DNAString or DNAStringSet object and simply 
#' shuffles the order individual nucleotides. 
#'
#' @param dna DNAString or DNAStringSet
#' @param order Integer, which order to take into consideration for shuffling
#' (requires Python ushuffle library)
#' @return A DNAString or DNAStringSet
#' 
#' @import Biostrings
#' @export
#' 
#' @examples
#' shuffleSeq('ACGTGGGCTATTAGCTACTGTACGTG')

shuffleSeq <- function(dna, order = 1) {
    if (is(dna, 'DNAString')) {
        dna <- Biostrings::DNAStringSet(dna)
    } 
    else if (is(dna, 'character')) {
        dna <- Biostrings::DNAStringSet(dna)
    }
    if (order == 1) {
        # message('Not using ushuffle')
        shuffled <- Biostrings::DNAStringSet(
            lapply(dna, function(seq) {
                charvec <- strsplit(as.character(seq),"")[[1]]
                shuffled_charvec <- sample(charvec)
                Biostrings::DNAString( paste(shuffled_charvec, collapse="") )
            })
        )
    }
    else {
        # message('Using ushuffle')
        if (!reticulate::py_module_available("ushuffle")) 
            reticulate::py_install("ushuffle")
        ushuffle <- reticulate::import("ushuffle")
        shuffled <- Biostrings::DNAStringSet(
            lapply(dna, function(seq) {
                seq <- ushuffle$shuffle(charToRaw(as.character(seq)), 2)
                Biostrings::DNAString( as.character(seq) )
            })
        )
    }
    return(shuffled)
}

#' A function to sample GRanges from GRanges/DNAStringSet
#'
#' This function takes a given GRanges and returns another GRanges 
#' object. The new GRanges has the same number of ranges and the same
#' chromosome, width and strand distributions than the original 
#' GRanges. 
#'
#' @param x GRanges or DNAStringSet object
#' @param n Integer, number of sampled GRanges
#' @param width Integer, width of sampled GRanges
#' @param exclude Boolean, should the original GRanges be excluded?
#' @param avoid_overlap Boolean, should the sampled GRanges not 
#' be overlapping?
#' @return A GRanges object of length n
#' 
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' sampleGRanges(ce11_proms, 100)

sampleGRanges <- function(
    x, 
    n = NULL, 
    width = NULL,
    exclude = FALSE, 
    avoid_overlap = FALSE
)
{
    UseMethod("sampleGRanges")
}

#' A function to sample GRanges within GRanges
#'
#' This function takes a given GRanges and returns another GRanges 
#' object. The new GRanges has the same number of ranges and the same
#' chromosome, width and strand distributions than the original 
#' GRanges. 
#'
#' @param x GRanges object
#' @param n Integer, number of sampled GRanges
#' @param width Integer, width of sampled GRanges
#' @param exclude Boolean, should the original GRanges be excluded?
#' @param avoid_overlap Boolean, should the sampled GRanges not 
#' be overlapping?
#' @return A GRanges object of length n
#' 
#' @importFrom methods as
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' sampleGRanges(ce11_proms, 100)

sampleGRanges.GRanges <- function(
    x, 
    n = NULL, 
    width = NULL, 
    exclude = FALSE, 
    avoid_overlap = FALSE
)
{
    granges <- x
    maxed_granges <- GenomicRanges::GRanges(
        levels(GenomicRanges::seqnames(granges)), 
        IRanges::IRanges(
            1, 
            width = lengths(GenomicRanges::coverage(granges))
        )
    )
    if (is.null(n)) n <- length(x)
    N <- 2*n 
    g <- GRanges()
    #
    while (length(g) < n) {
        # Sample chrs based on their size
        chrs <- unlist(lapply(
            seq_along(IRanges::width(maxed_granges)), 
            function(K) {
                rep(
                    as.character(GenomicRanges::seqnames(maxed_granges)[K]), 
                    GenomicRanges::end(maxed_granges)[K]
                )
            }
        ))
        chrs <- sample(chrs, N)
        chrs <- as.character(sort(factor(
            chrs, 
            levels = levels(GenomicRanges::seqnames(maxed_granges))
        )))
        # For each chr, get random positions within this chr.
        pos <- unlist(lapply(
            unique(chrs), 
            function(chr) {
                sample(
                    lengths(GenomicRanges::coverage(granges)[chr]), 
                    table(chrs)[chr]
                )
            }
        ))
        # For each chr, get random strands within this chr.
        strands <- unlist(RleList(lapply(
            unique(chrs), 
            function(chr) {
                sample(
                    as.vector(GenomicRanges::strand(
                        granges[GenomicRanges::seqnames(granges) == chr]
                    )), 
                    table(chrs)[chr], 
                    replace = TRUE
                )
            }
        )))
        # For each chr, get widths within this chr.
        if (is.null(width)) {
            widths <- as.vector(unlist(RleList(lapply(
                unique(chrs), 
                function(chr) {
                    sample(
                        IRanges::width(
                            granges[GenomicRanges::seqnames(granges) == chr]
                        ), 
                        table(chrs)[chr], 
                        replace = TRUE
                    )
                }
            ))))
        } 
        else {
            widths <- width
        }
        # Build a GRanges object 
        suppressWarnings({
            g2 <- GenomicRanges::GRanges(
                chrs, 
                IRanges::IRanges(pos, width = widths),
                strand = strands, 
                seqinfo = GenomeInfoDb::seqinfo(granges)
            )
            GenomeInfoDb::seqlengths(g2) <- lengths(maxed_granges)
        })
        if (avoid_overlap) {
            g2 <- reduce(g2)
        }
        # Remove regions overlapping with initial GRanges
        if (exclude) {
            g2 <- g2[!(IRanges::`%over%`(g2, granges))]
        }
        # Remove the extended granges not embedded within the initial granges
        g2 <- IRanges::subsetByOverlaps(
            g2, 
            maxed_granges, 
            type = 'within'
        )
        # g2 <- GenomicRanges::trim(g2)
        if (length(unique(widths)) == 1) {
            g2 <- g2[GenomicRanges::width(g2) == unique(widths)]
        }
        g <- c(g, g2)
    }
    g <- g[sample(seq_len(length(g)), n)]
    g <- sort(g)
    return(g)
}

#' A function to sample GRanges within DNAStringSet
#'
#' This function takes a DNAStringSet (for instance a set of sequences
#' obtained from a BSgenome object) and returns a GRanges object. 
#' The ranges are contained within the DNAStringSet object. This function
#' is useful to sample sequences from a full genome. 
#'
#' @param x a DNAStringSet object
#' @param ... Additional parameters
#' @return A GRanges object of length n
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' library(Biostrings)
#' seqs <- getSeq(BSgenome.Scerevisiae.UCSC.sacCer3::
#'     BSgenome.Scerevisiae.UCSC.sacCer3)
#' sampleGRanges(seqs, n = 100, w = 100)

sampleGRanges.DNAStringSet <- function(
    x, 
    ...
)
{
    seqs <- x
    granges <- DNAStringSet2GRanges(seqs)
    g <- sampleGRanges(granges, ...)
    g$seq <- seqs[g]
    return(g)
}

#' A function to sample GRanges from BSgenome
#'
#' This function takes a BSgenome and returns a GRanges object. 
#' This function is useful to sample sequences from a full genome. 
#'
#' @param x a BSgenome object
#' @param ... Additional parameters
#' @return A GRanges object of length n
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' g <- sampleGRanges(
#'     BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3, 
#'     w = 150, n = 100
#' )
#' g
#' g$seq

sampleGRanges.BSgenome <- function(
    x, 
    ...
)
{
    genome <- x
    granges <- GenomicRanges::GRanges(
        names(genome), 
        lengths(genome)
    )
    GenomeInfoDb::seqlengths(granges) <- lengths(genome)
    g <- sampleGRanges(granges, ...)
    g$seq <- Biostrings::getSeq(genome, g)
    return(g)
}

#' A function to sample GRanges from genome IDs
#'
#' This function takes a BSgenome UCSC ID and returns a GRanges object. 
#' This function is useful to sample sequences from a full genome. 
#' Currently, sacCer3, ce11, dm6, danRer10, mm10 and hg38 are 
#' supported.
#'
#' @param x a genome ID ()
#' @param ... Additional parameters
#' @return A GRanges object of length n
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @import BSgenome
#' @export
#' 
#' @examples
#' sampleGRanges('sacCer3', width = 150, n = 100)

sampleGRanges.character <- function(
    x, 
    ...
)
{
    if (x %in% c(
        'sacCer3', 'ce11', 'dm6', 'mm10', 'hg38', 'danRer10'
    ) | x %in% BSgenome::installed.genomes()) {
        genome <- BSgenome::getBSgenome(x)
    }
    else {
        return(stop(
            'Only sacCer3, ce11, dm6, mm10, hg38 and danRer10 are supported'
        ))
    }
    sampleGRanges(genome, ...)
}

#' A function to sample GRanges from genome IDs
#'
#' This function is an alias for the sampleGRanges.character method.
#' It takes a BSgenome ID and returns a GRanges object. 
#' This function is useful to sample sequences from a full genome. 
#' Currently, sacCer3, ce11, dm6, danRer10, mm10 and hg38 are 
#' supported. 
#'
#' @rdname sampleGRanges.character
#' @export
#' @examples
#' sampleGenome('sacCer3', w = 150, n = 100)

sampleGenome <- sampleGRanges.character

#' A function to generate a GRanges from a DNAStringSet
#'
#' This function takes a DNAStringSet (for instance a set of sequences
#' obtained from a BSgenome object) and converts it into 
#' a GRanges object. 
#'
#' @param seqs A DNAStringSet
#' @return A GRanges object
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' DNAStringSet2GRanges(
#'     BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
#' )

DNAStringSet2GRanges <- function(seqs) {
    if (
        length(names(seqs)) == 0 | 
        length(unique(names(seqs))) != length(names(seqs))
    ) {
        names(seqs) <- paste0('seq_', seq_len(length(seqs)))
    }
    g <- GenomicRanges::GRanges(
        names(seqs), 
        IRanges::IRanges(1, width = lengths(seqs))
    )
    GenomeInfoDb::seqlengths(g) <- lengths(seqs)
    return(g)
}

#' setUpBPPARAM
#' 
#' A function to dynamically select 
#' MulticoreParam or SnowParam (if Windows)
#'
#' @param nproc number of processors
#' @return A BPPARAM object
#' 
#' @import BiocParallel
#' @export
#' 
#' @examples
#' BPPARAM <- setUpBPPARAM(1)

setUpBPPARAM <- function(nproc = 1) {
    if (.Platform$OS.type == "windows") {
        # windows doesn't support multicore, using snow instead
        result <- BiocParallel::SnowParam(workers = nproc)
    } else {
        result <- BiocParallel::MulticoreParam(workers = nproc)
    }
    return(result)
}

#' checkPeriodFromFourier
#' 
#' A function to quickly check whether a period of interest is in 
#' the result of a FFT analysis, based on the length of the vector which 
#' is used for FFT
#'
#' @param len length of the vector used for FFT
#' @param period Float, period of interest
#' @return Float period returned by FFT closest to the period of interest
#' 
#' @importFrom stats spectrum
#' @importFrom stats runif
#' @export
#' 
#' @examples
#' checkPeriodFromFourier(300, 10)

checkPeriodFromFourier <- function(len, period) {
    s <- stats::spectrum(stats::runif(len), plot = FALSE)
    freq <- 1/period
    isFreqThere <- sum(abs(s$freq - freq) < 10e-6)
    if (isFreqThere) {
        return(1/s$freq[which(abs(s$freq - freq) < 10e-6)])
    } 
    else {
        freq2 <- s$freq[which.min(abs(s$freq - freq))]
        message(
            "Frequency closest to ", freq, " returned by FFT is: ", 
            formatC(freq2, digits = 3)
        )
        return(1/s$freq[which.min(abs(s$freq - freq))])
    }
}
