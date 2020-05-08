#' A function to easily coerce a named list into a long data.frame
#'
#' @param x A named list of vectors or data.frames
#' 
#' @return A long data frame

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

#' A tidyverse-compatible function 
#'     to replace NAs in a vector by a given value
#'
#' @param x vector
#' @param value Replace NAs by this variable
#' 
#' @return A vector with NA replaced by value

na.replace <- function(x, value) {
    which.na <- is.na(x)
    x[which.na] <- value
    return(x)
}

#' A tidyverse-compatible function to remove NAs from a vector
#'
#' @param x vector
#' 
#' @return A vector without NAs

na.remove <- function(x) {
    which.na <- is.na(x)
    return(x[!which.na])
}

#' A function to quickly shuffle sequence(s)
#'
#' @param dna DNAString or DNAStringSet
#' 
#' @return A DNAString or DNAStringSet
#' 
#' @import Biostrings
#' @export
#' 
#' @examples
#' \dontrun{
#'     shuffleSeq('ACGTGGGCTATTAGCTACTGTACGTG')
#' }

shuffleSeq <- function(dna) {
    if (is(dna, 'DNAString')) {
        charvec <- strsplit(as.character(dna),"")[[1]]
        shuffled_charvec <- sample(charvec)
        Biostrings::DNAString( paste(shuffled_charvec, collapse="") )
    } 
    else if (is(dna, 'DNAStringSet')) {
        Biostrings::DNAStringSet(lapply(dna, function(seq) {
            charvec <- strsplit(as.character(seq),"")[[1]]
            shuffled_charvec <- sample(charvec)
            Biostrings::DNAString( paste(shuffled_charvec, collapse="") )
        }))
    }
    else if (is(dna, 'character')) {
        dna <- Biostrings::DNAString(dna)
        charvec <- strsplit(as.character(dna),"")[[1]]
        shuffled_charvec <- sample(charvec)
        Biostrings::DNAString( paste(shuffled_charvec, collapse="") )
    }
}

#' A function to sample GRanges from GRanges/DNAStringSet
#'
#' @param x a GRanges or DNAStringSet object
#' @param n Integer number of sampled GRanges
#' @param width Integer width of sampled GRanges
#' @param exclude Boolean should the original GRanges be excluded?
#' @param avoid_overlap Boolean should the sampled GRanges not be overlapping?
#' 
#' @return A GRanges object of length n
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'     data(ce11_proms)
#'     sampleGRanges(ce11_proms, 100)
#'     #
#'     yeast_seq <- getSeq(
#'         BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
#'     )
#'     sacCer3_random_regions <- sampleGRanges(yeast_seq, n = 10000, w = 200)
#' }

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
#' @param x a GRanges object
#' @param n Integer number of sampled GRanges
#' @param width Integer width of sampled GRanges
#' @param exclude Boolean should the original GRanges be excluded?
#' @param avoid_overlap Boolean should the sampled GRanges not be overlapping?
#' 
#' @return A GRanges object of length n
#' 
#' @importFrom methods as
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' \dontrun{
#'     data(ce11_proms)
#'     sampleGRanges(ce11_proms)
#' }

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
#' @param x a DNAStringSet object
#' @param ... Additional parameters
#' 
#' @return A GRanges object of length n
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' 
#' @examples
#' \dontrun{
#'     data(ce11_proms)
#'     sampleGRanges(ce11_proms, 100)
#' }

sampleGRanges.DNAStringSet <- function(
    x, 
    ...
)
{
    seqs <- x
    granges <- DNAStringSet2GRanges(seqs)
    sampleGRanges(granges, ...)
}

#' A function to sample GRanges from BSgenome
#'
#' @param x a BSgenome object
#' @param ... Additional parameters
#' 
#' @return A GRanges object of length n
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' 
#' @examples
#' \dontrun{
#'     sampleGRanges(
#'         (BSgenome.Scerevisiae.UCSC.sacCer3::
#'             BSgenome.Scerevisiae.UCSC.sacCer3), 
#'         100
#'     )
#' }

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
#' @param x a genome ID ()
#' @param ... Additional parameters
#' 
#' @return A GRanges object of length n
#' 
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' 
#' @examples
#' \dontrun{
#'     sampleGRanges(
#'         'ce11', 
#'         100
#'     )
#' }

sampleGRanges.character <- function(
    x, 
    ...
)
{
    if (x %in% c(
        'sacCer3', 'ce11', 'dm6', 'mm10', 'hg38', 'danRer10'
    )) {
        genome <- char2BSgenome(x)
    }
    else {
        return(stop(
            'Only sacCer3, ce11, dm6, mm10, hg38 and danRer10 are supported'
        ))
    }
    sampleGRanges(genome, ...)
}

#' @rdname sampleGRanges.character
#' @export

sampleGenome <- sampleGRanges.character

#' A function to generate a GRanges from a DNAStringSet
#'
#' @param seqs A DNAStringSet
#' 
#' @return A GRanges object
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' \dontrun{
#'     seq <- BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
#'     DNAStringSet2GRanges(seq)
#' }

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

#' A function to get the BSgenome from genome ID
#'
#' @param ID character
#' 
#' @return BSgenome object
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'     genome <- char2BSgenome('ce11')
#'     genome
#' }

char2BSgenome <- function(ID) {
    genome <- switch(
        ID, 
        'sacCer3' = (BSgenome.Scerevisiae.UCSC.sacCer3::
            BSgenome.Scerevisiae.UCSC.sacCer3), 
        'ce11' = (BSgenome.Celegans.UCSC.ce11::
            BSgenome.Celegans.UCSC.ce11), 
        'dm6' = (BSgenome.Dmelanogaster.UCSC.dm6::
            BSgenome.Dmelanogaster.UCSC.dm6), 
        'danRer10' = (BSgenome.Drerio.UCSC.danRer10::
            BSgenome.Drerio.UCSC.danRer10), 
        'mm10' = (BSgenome.Mmusculus.UCSC.mm10::
            BSgenome.Mmusculus.UCSC.mm10), 
        'hg38' = (BSgenome.Hsapiens.UCSC.hg38::
            BSgenome.Hsapiens.UCSC.hg38)
    )
    return(genome)
}