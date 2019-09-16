#### ---- getPeriodicity methods ---- ####

#' Core function
#'
#' @param x a DNAStringSet, a DNAString of a GRanges
#' 
#' @return direct to the appropriate method

getPeriodicity <- function(x, ...) {
    UseMethod("getPeriodicity")
}

#' Core function
#'
#' @param x a DNAStringSet
#' 
#' @return List a list containing the results of getPeriodicity function. The 
#' dists vector is the raw vector of all distances between any possible dinucleotide. 
#' The hist data.frame is the distribution of distances over RANGE_FOR_SPECTRUM. 
#' The normalized_hist is the raw hist, normalized for decay over increasing distances.
#' The spectra object is the output of the FFT applied over normalized_hist. 
#' The PSD data frame is the power spectrum density scores over given frequencies.
#' The signal_to_noise_ratio is a data.frame containing enrichment scores of TT
#' periodicity, for the periods in the period vector. The motif object is the 
#' dinucleotide being analysed.

getPeriodicity.default <- function(seqs, ...) {
    getPeriodicity.DNAStringSet(seqs, ...)
}

#' Core function
#'
#' @param seqs a DNAStringSet
#' @param motif a dinucleotide of interest
#' @param RANGE_FOR_SPECTRUM Numeric vector The distances between nucleotides
#' to take into consideration when performing Fast Fourier Transform (default 1:50).
#' @param period Vector a numerical vector of periods to extract.
#' @param plot Boolean Should the FFT results be plotted? 
#' @param cores Integer How many processors should be used to split to work? 
#' @param verbose Boolean
#' 
#' @return The output of getPeriodicity

getPeriodicity.DNAStringSet <- function(
    seqs, 
    motif = 'WW', 
    RANGE_FOR_SPECTRUM = 1:100,
    period = seq(2, 20, 1),
    plot = FALSE,
    cores = 2, 
    verbose = TRUE
)
{
    # Get pairwise distances ---------------------------------------------------
    if (verbose) message("- Getting pairwise distances.")
    dists <- parallel::mclapply(1:length(seqs), function(k) {
        seq <- seqs[k]
        Biostrings::vmatchPattern(
            motif, 
            seq, 
            max.mismatch = 0
            , fixed = FALSE
        )[[1]] %>% 
            GenomicRanges::start() %>% 
            dist() %>% 
            c()
    }, mc.cores = cores) %>% 
        unlist()
    # Fourier ------------------------------------------------------------------
    if (verbose) message("- Applying Fast Fourier Transform to the vector of distances.")
    if (max(dists) < tail(RANGE_FOR_SPECTRUM, 1)) {
        if (verbose) message(
            'Range (', 
            head(RANGE_FOR_SPECTRUM, 1), ':', tail(RANGE_FOR_SPECTRUM, 1), 
            ') is wider than any range in the current distances vector. Shortening the range from ', 
            head(RANGE_FOR_SPECTRUM, 1), ' to ', max(dists), '...'
        )
        RANGE_FOR_SPECTRUM <- head(RANGE_FOR_SPECTRUM, 1):max(dists)
    }
    hist <- hist(dists, breaks = seq(1, max(dists)+1, 1), plot = FALSE)$counts
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
    mtm <- spectra$spec[unlist(lapply(c(1/(period)), function (freq) {
        idx <- which.min(abs(freq - spectra$freq))
    }))]
    ratios <- data.frame(
        period = period,
        SNR = sapply(1:length(mtm), function(x) {mtm[x]/mean(mtm[-x])})
    )
    # Return all results -------------------------------------------------------
    return(list(
        dists = dists, 
        hist = data.frame(
            distance = seq(1, max(dists), 1), 
            counts = hist
        ), 
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
        motif = motif
    ))
}

#' Core function - NOT WORKING. DO NOT USE
#'
#' @param seq a DNAString
#' @param motif a dinucleotide of interest
#' @param subseq_len The length of sub-sequences
#' @param n_occurences The targeted number of pairs of dinucleotides to obtain. 
#' 
#' @return List NULL

getPeriodicity.DNAString <- function(seq, motif = 'TT', subseq_len = NULL, n_occurences = 200, ...) {
    if(is.null(subseq_len)) subseq_len <- nchar(seq) * 0.25
    n_motif <- length(Biostrings::matchPattern(motif, seq, fixed = FALSE))
    repets <- round(n_occurences/n_motif * 10)
    starts <- sample(1:(length(seq)-subseq_len), repets, replace = repets > (nchar(seq) - subseq_len))
    seqs <- lapply(seq_along(starts), function(K) {seq[starts[K]:(starts[K]+subseq_len)]}) %>% as('DNAStringSet')
    spectra <- getPeriodicity(seqs, ...)$spectra
    m <- max(
        spectra$spec[
        (which.min(abs(spectra$freq - freq))-1):(which.min(abs(spectra$freq - freq))+1)
        ]
    )
    return(m)
}

#' Core function
#'
#' @param granges a GRanges
#' @param genome DNAStringSet object. Ideally, the sequence of an entire genome, 
#' obtained for instance by running 
#' \code{Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11)}.
#' @param ... other parameters required in getPeriodicity
#'
#' @return The output of getPeriodicity

getPeriodicity.GRanges <- function(
    granges, 
    genome = Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11), 
    ...
)
{
    if (is.null(granges$seq)) granges <- withSeq(granges, genome)
    seqs <- granges$seq
    getPeriodicity(seqs, ...)
}

#' Internal function
#'
#' @param hist Vector a numeric vector
#' 
#' @return a normalized vector


normalizeHistogram <- function(hist) {
    h <- hist
    h <- h / sum(h)
    smoothed.h <- c(zoo::rollmean(h, k = 10), rep(0, 9))
    norm.h <- scale(h - smoothed.h)
    return(norm.h)
}
