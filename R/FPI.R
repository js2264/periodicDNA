#' A function to compute the fold power increase, as introduced by 
#'   Pich et al., Cell 2018.
#'
#' @param seqs DNAStringSet sequences of interest
#' @param motif character k-mer of interest
#' @param period integer Period of interest
#' @param n_shuffling integer Number of shuffling
#' @param cores_shuffling integer Number of threads to use to split
#'   shuffling
#' @param verbose integer Should the function be verbose? 
#' @param ... Additional arguments
#' 
#' @return list Several metrics included (FPI, observed PSD, etc...)
#' 
#' @import GenomicRanges
#' @import IRanges
#' @export

getFPI <- function(
    seqs, 
    motif,
    period = 10, 
    n_shuffling = 10,
    cores_shuffling = 10,
    verbose = 1,
    ...
) {
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