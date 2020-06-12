#' A function to compute fold power increase (FPI)
#' 
#' This function takes a set of sequences and a k-mer of interest, 
#' along with a given period. It calculates the FPI, as introduced
#' by Pich et al., Cell 2018.
#' 
#' The FPI is calculated as the ratio of the observed k-mer PSD 
#' at the chosen perdiod subtracted by the median k-mer PSD
#' in shuffled sequences, divided by the median k-mer PSD 
#' in shuffled sequences. The number of shufflings is specified 
#' by the n_shuffling argument. 
#'
#' @param x DNAStringSet, sequences of interest
#' @param ... Additional arguments
#' @return Several metrics including FPI, observed PSD, etc...
#' 
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' fpi <- getFPI(
#'     ce11_proms_seqs[1:10], 
#'     genome = 'ce11', 
#'     motif = 'TT', 
#'     cores_shuffling = 1
#' )
#' fpi$FPI
#' fpi$observed_PSD
#' fpi$shuffled_PSD
#' plotFPI(fpi)

getFPI <- function(x, ...) {
    UseMethod("getFPI")
}

#' A function to compute fold power increase (FPI)
#' 
#' This function takes a set of sequences and a k-mer of interest, 
#' along with a given period. It calculates the FPI, as introduced
#' by Pich et al., Cell 2018.
#' 
#' The FPI is calculated as the ratio of the observed k-mer PSD 
#' at the chosen perdiod subtracted by the median k-mer PSD
#' in shuffled sequences, divided by the median k-mer PSD 
#' in shuffled sequences. The number of shufflings is specified 
#' by the n_shuffling argument. 
#'
#' @param x DNAStringSet, sequences of interest
#' @param motif character, k-mer of interest
#' @param period integer, Period of interest
#' @param n_shuffling integer, Number of shuffling
#' @param cores_shuffling integer, Number of threads to use to split
#' shuffling
#' @param verbose integer, Should the function be verbose? 
#' @param ... Additional arguments
#' @return Several metrics including FPI, observed PSD, etc...
#' 
#' @import GenomicRanges
#' @import IRanges
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' fpi <- getFPI(
#'     ce11_proms_seqs[1:10], 
#'     motif = 'TT', 
#'     cores_shuffling = 1, 
#'     cores = 1
#' )
#' fpi$FPI
#' fpi$observed_PSD
#' fpi$shuffled_PSD
#' plotFPI(fpi)

getFPI.DNAStringSet <- function(
    x, 
    motif,
    period = 10, 
    n_shuffling = 10,
    cores_shuffling = 10,
    verbose = 1,
    ...
)
{
    seqs <- x
    
    # Calculating observed PSD ---------------------------------------
    if (verbose) message('- Calculating observed PSD')
    obs <- getPeriodicity(
        seqs, 
        motif,
        skip_shuffling = TRUE,
        verbose = 1,
        ...
    )
    obs_PSD <- obs$PSD$PSD[which.min(abs(1/obs$PSD$freq - period))]
    if (verbose) message('>> Measured PSD @ ', period, 'bp is: ', obs_PSD)
    # Shuffling sequences and re-computing ---------------------------
    l_shuff <- parallel::mclapply(
        mc.cores = cores_shuffling, 
        seq_len(n_shuffling), 
        function(k) {
            if (verbose) message('- Shuffling ', k, '/', n_shuffling)
            shuff_seqs <- shuffleSeq(seqs)
            shuff <- getPeriodicity(
                shuff_seqs, 
                motif,
                skip_shuffling = TRUE,
                verbose = 0,
                ...
            )
            return(shuff)
        }
    )
    l_shuff_PSD <- lapply(l_shuff, function(shuff) {
        shuff$PSD$PSD[which.min(abs(1/shuff$PSD$freq - period))]
    }) %>% unlist()
    # Calculate FPI --------------------------------------------------
    fpi <- (obs_PSD - median(l_shuff_PSD)) / median(l_shuff_PSD)
    if (verbose) message('>> Calculated FPI @ ', period, 'bp is: ', fpi)
    # Return FPI -----------------------------------------------------
    res <- list(
        FPI = fpi, 
        observed_PSD = obs_PSD, 
        observed_spectra = obs, 
        shuffled_PSD = l_shuff_PSD, 
        shuffled_spectra = l_shuff, 
        motif = motif, 
        period = period
    )
    return(res)
}

#' A function to compute fold power increase (FPI)
#' 
#' This function takes a GRanges and a k-mer of interest, 
#' along with a given period. It calculates the FPI, as introduced
#' by Pich et al., Cell 2018.
#' 
#' The FPI is calculated as the ratio of the observed k-mer PSD 
#' at the chosen perdiod subtracted by the median k-mer PSD
#' in shuffled sequences, divided by the median k-mer PSD 
#' in shuffled sequences. The number of shufflings is specified 
#' by the n_shuffling argument. 
#'
#' @param x GRanges, GRanges of interest
#' @param genome Genome ID or BSgenome
#' @param ... Additional arguments
#' @return Several metrics including FPI, observed PSD, etc...
#' 
#' @import GenomicRanges
#' @import IRanges
#' @export
#' 
#' @examples
#' data(ce11_TSSs)
#' fpi <- getFPI(
#'     ce11_TSSs[['Ubiq.']][1:10], 
#'     genome = 'ce11', 
#'     motif = 'TT', 
#'     cores_shuffling = 1
#' )
#' fpi$FPI
#' fpi$observed_PSD
#' fpi$shuffled_PSD
#' p1 <- plotFPI(fpi)
#' fpi <- getFPI(
#'     ce11_TSSs[['Muscle']][1:10], 
#'     genome = 'ce11', 
#'     motif = 'TT', 
#'     cores_shuffling = 1
#' )
#' fpi$FPI
#' fpi$observed_PSD
#' fpi$shuffled_PSD
#' p2 <- plotFPI(fpi)
#' cowplot::plot_grid(p1, p2)

getFPI.GRanges <- function(
    x,
    genome, 
    ...
)
{
    granges <- x
    
    if (!is.null(granges$seq)) {
        seqs <- granges$seq
    }
    else {
        seqs <- withSeq(granges, genome)$seq
    }
    getFPI(seqs, ...)
}

