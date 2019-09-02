getPeriodicity <- function(
    granges, 
    motif = 'WW', 
    N_MISMATCH = 0, 
    RANGE_FOR_SPECTRUM = 1:100,
    freq = 0.10, 
    period = seq(2, 20, 1),
    plot = FALSE,
    cores = 2, 
    verbose = TRUE) {
    # Get pairwise distances ---------------------------------------------------
    if (verbose) message("- Getting pairwise distances.")
    seqs <- granges$seq
    dists <- parallel::mclapply(1:length(seqs), function(k) {
        seq <- seqs[k]
        Biostrings::vmatchPattern(motif, seq, max.mismatch = N_MISMATCH, fixed = FALSE)[[1]] %>% 
            GenomicRanges::start() %>% 
            dist() %>% 
            c()
    }, mc.cores = cores) %>% 
        unlist()
    # Fourier ------------------------------------------------------------------
    if (verbose) message("- Applying Fast Fourier Transform to the vector of distances.")
    if (max(dists) < tail(RANGE_FOR_SPECTRUM, 1)) {
        if (verbose) message('Range (', head(RANGE_FOR_SPECTRUM, 1), ':', tail(RANGE_FOR_SPECTRUM, 1), ') is wider than any range in the current distances vector. Shortening the range from ', head(RANGE_FOR_SPECTRUM, 1), ' to ', max(dists), '...')
        RANGE_FOR_SPECTRUM <- head(RANGE_FOR_SPECTRUM, 1):max(dists)
    }
    hist <- hist(dists, breaks = seq(1, max(dists)+1, 1), plot = FALSE) %$%
        counts
    norm_hist <- normalizeHistogram(hist)
    spectra <- norm_hist %>%
            '['(RANGE_FOR_SPECTRUM) %>% 
            stats::spectrum(plot = FALSE)
    if (plot) {
        if (verbose) message("- Plotting results.")
        plot(spectra)
    }
    # Get signal-to-noise ratio ------------------------------------------------
    if (verbose) message("- Computing signal-to-noise ratios.")
    mtm <- spectra$spec[unlist(lapply(c(1/(period)), function (freq) {idx <- which.min(abs(freq - spectra$freq))}))]
    ratios <- data.frame(
        period = period,
        SNR = sapply(1:length(mtm), function(x) {mtm[x]/mean(mtm[-x])})
    )
    # Return all results -------------------------------------------------------
    return(list(
        dists = dists, 
        hist = hist, 
        normalized_hist = data.frame(
            distance = seq(1, max(dists), 1), 
            norm_counts = norm_hist
        ),
        spectra = spectra, 
        PSD = data.frame(
            freq = spectra$freq, 
            PSD = spectra$spec
        ),
        signal_to_noise_ratio = ratios, 
        motif = motif,
        freq = freq
    ))
}
normalizeHistogram <- function(hist) {
    h <- hist
    h <- h / sum(h)
    smoothed.h <- c(zoo::rollmean(h, k = 10), rep(0, 9))
    norm.h <- scale(h - smoothed.h)
    return(norm.h)
}
