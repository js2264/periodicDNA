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

