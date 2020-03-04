#' Function to compute the overall periodicity of a motif over sequence(s). 
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
#' 
#' @export

getPeriodicity <- function(x, ...) {
    UseMethod("getPeriodicity")
}

#' Function to compute the overall periodicity of a motif over a set of 
#'   sequences.
#'
#' @param seqs a DNAStringSet
#' @param bg_seqs a DNAStringSet (ideally from random loci) for genome 
#'   background model
#' @param model_background Boolean Should the genome background be taken into 
#'   account during normalization?
#' @param motif a dinucleotide of interest
#' @param RANGE_FOR_SPECTRUM Numeric vector The distances between nucleotides
#' to take into consideration when performing Fast Fourier Transform.
#' @param period Vector a numerical vector of periods to extract.
#' @param plot Boolean Should the FFT results be plotted? 
#' @param cores Integer How many processors should be used to parallelize
#'   the mapping? 
#' @param roll Integer Window to smooth the distribution of pairwise distances
#'   (default: 3, to discard the 3-bp periodicity of dinucleotides which 
#'   can be very strong in vertebrate genomes)
#' @param verbose Boolean
#' @param sample Integer if > 0, will randomly sample this many integers
#'   from the dists vector before normalization. This ensures consistency 
#'   when looking at periodicity in different genomes, since different genomes
#'   will have different GC%
#' @import parallel
#' @import Biostrings
#' @import IRanges
#' @importFrom stats spectrum
#' @export
#' @return List a list containing the results of getPeriodicity function. The 
#' dists vector is the raw vector of all pairwise distances between dinucleotides. 
#' The hist dataframe is the distribution of distances over RANGE_FOR_SPECTRUM. 
#' The normalized_hist is the raw histogram normalized for 
#'   decay over increasing distances between pairs of dinucleotides.
#' The spectra object is the output of the FFT applied over normalized_hist. 
#' The PSD dataframe contains power spectrum density scores over given 
#'   frequencies.
#' The signal_to_noise_ratio is a dataframe containing enrichment scores of
#' dinucleotide periodicity, for each period in the period vector. 
#' The motif character is the analysed dinucleotide.

getPeriodicity.DNAStringSet <- function(
    seqs, 
    bg_seqs = NULL,
    model_background = FALSE, 
    motif = 'WW', 
    RANGE_FOR_SPECTRUM = 1:200,
    period = seq(2, 20, 1),
    plot = FALSE,
    cores = 2, 
    roll = 3,
    verbose = TRUE, 
    sample = 0, 
    doZscore = TRUE
)
{
    # Get pairwise distances ---------------------------------------------------
    if (verbose) message("- Mapping k-mers.")
    dists <- parallel::mclapply(1:length(seqs), function(k) {
        seq <- seqs[k]
        Biostrings::vmatchPattern(
            motif, 
            seq, 
            max.mismatch = 0, 
            fixed = FALSE
        )[[1]] %>% 
            IRanges::start() %>% 
            dist() %>% 
            c()
    }, mc.cores = cores) %>% unlist()
    if (length(dists) < 10) {
        if (verbose) message("- Only ", length(dists), " pairs of k-mers found. Returning null results")
        return(list(
            dists = dists, 
            hist = 0, 
            normalized_hist = 0,
            spectra = 0, 
            PSD = 0,
            signal_to_noise_ratio = 0, 
            motif = motif
        ))
    }
    if (verbose) message("- ", length(dists), " k-mers found.")
    if (sample < length(dists) & sample > 0 & is.null(bg_seqs)) {
        if (verbose) message("- Sampling ", sample, " k-mers.")
        set.seed(222) 
        dists <- dists[sample(1:length(dists), sample)]
    }
    if (verbose) message("- Calculating pairwise distances.")
    hist <- hist(dists, breaks = seq(1, max(dists)+1, 1), plot = FALSE)$counts
    if (!is.null(bg_seqs) & model_background) {
        if (verbose) message("- Mapping k-mers in background.")
        bg_dists <- parallel::mclapply(1:length(bg_seqs), function(k) {
            seq <- bg_seqs[k]
            Biostrings::vmatchPattern(
                motif, 
                seq, 
                max.mismatch = 0, 
                fixed = FALSE
            )[[1]] %>% 
                IRanges::start() %>% 
                dist() %>% 
                c()
        }, mc.cores = cores) %>% unlist()
        if (verbose) message("- Calculating pairwise distances in background.")
        bg_hist <- hist(bg_dists, breaks = seq(1, max(bg_dists)+1, 1), plot = FALSE)$counts
    } 
    else {bg_hist <- NULL}
    if (verbose) message("- Normalizing histogram vector.")
    if (length(hist) > 10) {
        norm_hist <- normalizeHistogram(hist, roll, doZscore)
    } 
    else {
        norm_hist <- hist
    }
    if (!is.null(bg_hist) & model_background) {
        if (verbose) message("- Normalizing histogram vector in background.")
        norm_bg_hist <- normalizeHistogram(bg_hist, roll)
        if (verbose) message("- Substracting background.")
        norm_hist <- norm_hist - norm_bg_hist
    }
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
    spectra <- norm_hist %>%
            '['(RANGE_FOR_SPECTRUM) %>% 
            stats::spectrum(plot = FALSE)
    PSD <- data.frame(
        freq = spectra$freq, 
        PSD = spectra$spec
    )
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
        PSD = PSD,
        signal_to_noise_ratio = ratios, 
        motif = motif
    ))
}

# #' Core function - NOT WORKING. DO NOT USE
# #'
# #' @param seq a DNAString
# #' @param motif a dinucleotide of interest
# #' @param subseq_len The length of sub-sequences
# #' @param n_occurences The targeted number of pairs of dinucleotides to obtain. 
# #' 
# #' @export
# #' @return List NULL
# 
# getPeriodicity.DNAString <- function(seq, motif = 'TT', subseq_len = NULL, n_occurences = 200, force_use = FALSE, ...) {
#     # STOP ERROR ---------------------------------------------------------------
#     if (!force_use) stop('FUNCTION NOT WORKING YET. ABORTING NOW.')
#     #
#     if(is.null(subseq_len)) subseq_len <- nchar(seq) * 0.25
#     n_motif <- length(Biostrings::matchPattern(motif, seq, fixed = FALSE))
#     repets <- round(n_occurences/n_motif * 10)
#     starts <- sample(1:(length(seq)-subseq_len), repets, replace = repets > (nchar(seq) - subseq_len))
#     seqs <- lapply(seq_along(starts), function(K) {seq[starts[K]:(starts[K]+subseq_len)]}) %>% as('DNAStringSet')
#     spectra <- getPeriodicity(seqs, ...)$spectra
#     m <- max(
#         spectra$spec[
#         (which.min(abs(spectra$freq - freq))-1):(which.min(abs(spectra$freq - freq))+1)
#         ]
#     )
#     return(m)
# }

#' Core function
#'
#' @param granges a GRanges
#' @param bg a GRanges to estimate background periodicity. 
#' @param genome DNAStringSet object. The sequence of an entire genome, 
#'   obtained for instance by running 
#'   \code{Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11)}.
#' @param ... other parameters forwarded to getPeriodicity.DNAStringSet()
#'
#' @return List a list containing the results of getPeriodicity function. The 
#' dists vector is the raw vector of all pairwise distances between dinucleotides. 
#' The hist dataframe is the distribution of distances over RANGE_FOR_SPECTRUM. 
#' The normalized_hist is the raw histogram normalized for 
#'   decay over increasing distances between pairs of dinucleotides.
#' The spectra object is the output of the FFT applied over normalized_hist. 
#' The PSD dataframe contains power spectrum density scores over given 
#'   frequencies.
#' The signal_to_noise_ratio is a dataframe containing enrichment scores of
#' dinucleotide periodicity, for each period in the period vector. 
#' The motif character is the analysed dinucleotide.
#' 
#' @export

getPeriodicity.GRanges <- function(
    granges,
    bg = NULL, 
    genome = Biostrings::getSeq(BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11), 
    ...
)
{
    seqs <- withSeq(granges, genome)$seq
    if (!is.null(bg)) {
        bg_seqs <- withSeq(bg, genome)$seq
    } else {
        bg_seqs <- NULL
    }
    getPeriodicity(seqs, bg_seqs, ...)
}

#' Internal function to normalize a pairwise distance 
#'
#' @param hist Vector a numeric vector
#' 
#' @importFrom zoo rollmean
#' @export
#' @return a normalized vector

normalizeHistogram <- function(hist, roll = 1, doZscore = TRUE, roll_smoothed.h = 10) {
    h <- hist
    h <- h / sum(h) # Normalize to total number of pairwise distances
    smoothed.h <- c(
        zoo::rollmean(h, k = roll_smoothed.h), 
        rep(0, roll_smoothed.h-1)
    ) # Smoothed distribution
    if (doZscore) {
        norm.h <- scale(
            h - smoothed.h
        ) # Substract smoothed distribution and then z-score
    } 
    else {
        norm.h <- h - smoothed.h
    }
    norm.h <- zoo::rollmean(
        norm.h, k = roll, na.pad = TRUE, align = 'center'
    ) # Smooth the normalized distribution 
    norm.h[is.na(norm.h)] <- 0
    return(norm.h)
}
