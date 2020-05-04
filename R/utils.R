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

#' A function to sample GRanges within GRanges
#'
#' @param granges a GRanges object
#' @param n Integer
#' 
#' @return A GRanges object of length n
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' \dontrun{
#'     data(ce_proms)
#'     sampleGRanges(ce_proms, 100)
#' }

sampleGRanges <- function(granges, n){
    rand_c <- sample(
        length(granges),
        n, 
        replace = TRUE, 
        prob = GenomicRanges::width(granges)
    )
    rand_ranges <- granges[rand_c]
    rand <- unlist(
        lapply(GenomicRanges::width(rand_ranges), sample, size=1)
    )
    pos <- GenomicRanges::start(rand_ranges)+rand-1
    res <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(rand_ranges), 
        IRanges::IRanges(pos, width=1),
        strand='*', 
        seqinfo = GenomeInfoDb::seqinfo(granges), 
        GenomicRanges::mcols(rand_ranges)
    )
    res <- sort(res)
    return(res)
}

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
    g <- GenomicRanges::GRanges(
        names(seqs), 
        IRanges::IRanges(1, width = lengths(seqs))
    )
    GenomeInfoDb::seqlengths(g) <- lengths(seqs)
    return(g)
}