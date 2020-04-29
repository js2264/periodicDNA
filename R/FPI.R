FPI <- function(
    seqs, 
    motif,
    period = 10, 
    n_shuffle = 10,
    parallel_shuffling = 10,
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
        mc.cores = parallel_shuffling, 
        1:n_shuffle, 
        function(k) {
            if (verbose) message('- Shuffling ', k, '/', n_shuffle)
            shuff_seqs <- shuffleSeq(seqs, k)
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
    # Return FPI -----------------------------------------------------
    return(list(
        FPI = fpi, 
        observed_PSD = obs_PSD, 
        observed_spectra = obs, 
        shuffled_PSD = l_shuff_PSD, 
        shuffled_spectra = l_shuff, 
        motif = motif, 
        period = period
    ))
}