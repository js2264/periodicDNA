#' Function to compute the overall periodicity of a 
#'     motif over sequence(s). 
#'
#' @param x a DNAStringSet, or a GRanges
#' @param ... additional parameters
#' 
#' @return List a list containing the results of getPeriodicity function. 
#'     The dists vector is the raw vector of all distances between any 
#'     possible dinucleotide. The hist data.frame is the distribution of 
#'     distances over RANGE_FOR_SPECTRUM. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spect.  density scores over given 
#'     frequencies. The signal_to_noise_ratio is a data.frame containing 
#'     enrichment scores of TT periodicity, for the periods in the period 
#'     vector. The motif object is the dinucleotide being analysed.
#' 
#' @export

getPeriodicity <- function(x, ...) {
    UseMethod("getPeriodicity")
}

#' Function to compute the overall periodicity of a motif over a set of 
#'     sequences.
#' 
#' @param x a DNAStringSet
#' @param motif a dinucleotide of interest
#' @param RANGE_FOR_SPECTRUM Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform.
#' @param period Vector a numerical vector of periods to extract.
#' @param plot Boolean Should the FFT results be plotted? 
#' @param cores Integer How many processors should be used to parallelize
#'     the mapping? 
#' @param roll Integer Window to smooth the distribution of pairwise distances
#'     (default: 3, to discard the 3-bp periodicity of dinucleotides which 
#'     can be very strong in vertebrate genomes)
#' @param verbose Boolean
#' @param sample Integer if > 0, will randomly sample this many integers
#'     from the dists vector before normalization. This ensures consistency 
#'     when looking at periodicity in different genomes, since different
#'     genomes will have different GC percent
#' @param doZscore Boolean should the normalized dampened signal be z-scored?
#' @param skip_shuffling Boolean should the shuffling sequences be done?
#' @param ... additional parameters
#' 
#' @return List a list containing the results of getPeriodicity function. 
#'     The dists vector is the raw vector of all distances between any 
#'     possible dinucleotide. The hist data.frame is the distribution of 
#'     distances over RANGE_FOR_SPECTRUM. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spect.  density scores over given 
#'     frequencies. The signal_to_noise_ratio is a data.frame containing 
#'     enrichment scores of TT periodicity, for the periods in the period 
#'     vector. The motif object is the dinucleotide being analysed.
#' 
#' @importFrom parallel mclapply
#' @import Biostrings
#' @import IRanges
#' @import magrittr
#' @importFrom stats spectrum
#' 
#' @export

getPeriodicity.DNAStringSet <- function(
    x, 
    motif = 'WW', 
    RANGE_FOR_SPECTRUM = seq_len(200),
    period = seq(2, 20, 1),
    plot = FALSE,
    cores = 2, 
    roll = 3,
    verbose = TRUE, 
    sample = 0, 
    doZscore = FALSE,
    skip_shuffling = TRUE, 
    ...
)
{
    seqs <- x
    
    # Get pairwise distances --------------------------------------------------
    if (verbose) message("- Mapping k-mers.")
    dists <- parallel::mclapply(seq_len(length(seqs)), function(k) {
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
        if (verbose) message("- Only ", length(dists), 
        " pairs of k-mers found. Returning null results")
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
    if (sample < length(dists) & sample > 0) {
        if (verbose) message("- Sampling ", sample, " k-mers.")
        dists <- dists[sample(seq_len(dists), sample)]
    }
    if (verbose) message("- Calculating pairwise distances.")
    hist <- hist(dists, breaks = seq(1, max(dists)+1, 1), plot = FALSE)$counts
    if (verbose) message("- Normalizing histogram vector.")
    if (length(hist) > 10) {
        norm_hist <- normalizeHistogram(hist, roll, doZscore)
    } 
    else {
        norm_hist <- hist
    }
    # Fourier -----------------------------------------------------------------
    if (verbose) message(
        "- Applying Fast Fourier Transform to the vector of distances."
    )
    if (max(dists) < tail(RANGE_FOR_SPECTRUM, 1)) {
        if (verbose) message(
            'Range (', 
            head(RANGE_FOR_SPECTRUM, 1), ':', tail(RANGE_FOR_SPECTRUM, 1), 
            ') is wider than any range in the current 
            distances vector. Shortening the range from ', 
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
    # Get signal-to-noise ratio -----------------------------------------------
    if (verbose) message("- Computing signal-to-noise ratios.")
    mtm <- spectra$spec[unlist(lapply(c(1/(period)), function (freq) {
        idx <- which.min(abs(freq - spectra$freq))
    }))]
    ratios <- data.frame(
        period = period,
        SNR = unlist(
            lapply(seq_along(mtm), function(x) {mtm[x]/mean(mtm[-x])})
        )
    )
    # Return all results ------------------------------------------------------
    if (skip_shuffling) {
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
    #
    # DO THE SAME FOR SHUFFLED SEQUENCES --------------------------------------
    {
        if (verbose) message("- SHUFFLING: suffling sequences.")
        seqs <- shuffleSeq(seqs)
        if (verbose) message("- SHUFFLING: Mapping k-mers.")
        dists_shuffled <- parallel::mclapply(seq_len(length(seqs)), function(k) {
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
        if (length(dists_shuffled) < 10) {
            if (verbose) message(
                "- Only ", 
                length(dists_shuffled), 
                " pairs of k-mers found. Returning null results"
            )
            return(list(
                dists = dists_shuffled, 
                hist = 0, 
                normalized_hist = 0,
                spectra = 0, 
                PSD = 0,
                signal_to_noise_ratio = 0, 
                motif = motif
            ))
        }
        if (verbose) message(
            "- SHUFFLING: ", length(dists_shuffled), " k-mers found."
        )
        if (sample < length(dists_shuffled) & sample > 0) {
            if (verbose) message("- Sampling ", sample, " k-mers.")
            dists_shuffled <- dists_shuffled[
                sample(seq_len(dists_shuffled), sample)
            ]
        }
        if (verbose) message("- SHUFFLING: Calculating pairwise distances.")
        hist_shuffled <- hist(
            dists_shuffled, 
            breaks = seq(1, max(dists_shuffled)+1, 1), 
            plot = FALSE
        )$counts
        if (verbose) message("- SHUFFLING: Normalizing histogram vector.")
        if (length(hist_shuffled) > 10) {
            norm_hist_shuffled <- normalizeHistogram(
                hist_shuffled, roll, doZscore
            )
        } 
        else {
            norm_hist_shuffled <- hist_shuffled
        }
        # Fourier -------------------------------------------------------------
        if (verbose) message(
            "- SHUFFLING: Applying Fast Fourier Transform 
            to the vector of distances."
        )
        if (max(dists_shuffled) < tail(RANGE_FOR_SPECTRUM, 1)) {
            if (verbose) message(
                'Range (', 
                head(RANGE_FOR_SPECTRUM, 1), ':', tail(RANGE_FOR_SPECTRUM, 1), 
                ') is wider than any range in the current distances vector. 
                Shortening the range from ', 
                head(RANGE_FOR_SPECTRUM, 1), ' to ', max(dists_shuffled), '...'
            )
            a <- head(RANGE_FOR_SPECTRUM, 1)
            RANGE_FOR_SPECTRUM <- a:max(dists_shuffled)
        }
        spectra_shuffled <- norm_hist_shuffled %>%
                '['(RANGE_FOR_SPECTRUM) %>% 
                stats::spectrum(plot = FALSE)
        PSD_shuffled <- data.frame(
            freq = spectra_shuffled$freq, 
            PSD = spectra_shuffled$spec
        )
        # Get signal-to-noise ratio -------------------------------------------
        if (verbose) message("- SHUFFLING: Computing signal-to-noise ratios.")
        mtm_shuffled <- spectra_shuffled$spec[
            unlist(lapply(c(1/(period)), function (freq) {
                idx <- which.min(abs(freq - spectra_shuffled$freq))
            }))
        ]
        ratios_shuffled <- data.frame(
            period = period,
            SNR = unlist(
                lapply(
                    seq_along(mtm_shuffled), 
                    function(x) {
                        mtm_shuffled[x]/mean(mtm_shuffled[-x])
                    }
                )
            )
        )
    }
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
        dists_shuffled = dists_shuffled, 
        hist_shuffled = data.frame(
            distance = seq(1, max(dists_shuffled), 1), 
            counts = hist_shuffled
        ), 
        normalized_hist_shuffled = data.frame(
            distance = seq(1, max(dists_shuffled), 1), 
            norm_counts = norm_hist_shuffled
        ),
        spectra_shuffled = spectra_shuffled, 
        PSD_shuffled = PSD_shuffled,
        signal_to_noise_ratio_shuffled = ratios_shuffled, 
        motif = motif
    ))
}

# #' Core function - NOT WORKING. DO NOT USE
# #'
# #' @param seq a DNAString
# #' @param motif a dinucleotide of interest
# #' @param subseq_len The length of sub-sequences
# #' @param n_occurences The targeted number of pairs of dinucleotides
# #' 
# #' @export
# #' @return List NULL
# 
# getPeriodicity.DNAString <- function(
#     seq, motif = 'TT', 
#     subseq_len = NULL,
#     n_occurences = 200,
#     force_use = FALSE, 
#     ...
# ) 
# {
#     # STOP ERROR ------------------------------------------------------------
#     if (!force_use) stop('FUNCTION NOT WORKING YET. ABORTING NOW.')
#     #
#     if(is.null(subseq_len)) subseq_len <- nchar(seq) * 0.25
#     n_motif <- length(Biostrings::matchPattern(motif, seq, fixed = FALSE))
#     repets <- round(n_occurences/n_motif * 10)
#     starts <- sample(
#       seq_len(length(seq)-subseq_len), repets, 
#       replace = repets > (nchar(seq) - subseq_len)
#     )
#     seqs <- lapply(
#       seq_along(starts), 
#       function(K) {seq[starts[K]:(starts[K]+subseq_len)]}
#     ) %>% as('DNAStringSet')
#     spectra <- getPeriodicity(seqs, ...)$spectra
#     a <- (which.min(abs(spectra$freq - freq))-1)
#     b <- (which.min(abs(spectra$freq - freq))+1)
#     m <- max(
#       spectra$spec[
#         a:b
#       ]
#     )
#     return(m)
# }

#' Core function
#' 
#' @param x a GRanges
#' @param genome DNAStringSet object. The sequence of an entire genome, 
#'     obtained for instance by running 
#'     \code{Biostrings::getSeq(
#'         BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
#'     )}.
#' @param ... other parameters forwarded to getPeriodicity.DNAStringSet()
#'
#' @return List a list containing the results of getPeriodicity function. 
#'     The dists vector is the raw vector of all distances between any 
#'     possible dinucleotide. The hist data.frame is the distribution of 
#'     distances over RANGE_FOR_SPECTRUM. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spect.  density scores over given 
#'     frequencies. The signal_to_noise_ratio is a data.frame containing 
#'     enrichment scores of TT periodicity, for the periods in the period 
#'     vector. The motif object is the dinucleotide being analysed.
#' 
#' @importFrom methods is
#' 
#' @export

getPeriodicity.GRanges <- function(
    x,
    genome = 'ce11',
    ...
)
{
    granges <- x
    
    if (methods::is(genome, 'character')) {
        if (genome %in% c(
            'sacCer3', 'ce11', 'dm6', 'mm10', 'hg38', 'danRer10'
        )) {
            genome <- switch(
                genome, 
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
        }
        else {
            return(stop(
                'Only ce11, dm6, mm10, hg38 and danRer10 are supported'
            ))
        }
    }
    if (methods::is(genome, 'BSgenome')) {
        genome <- Biostrings::getSeq(genome)
    }
    seqs <- withSeq(granges, genome)$seq
    getPeriodicity(seqs, ...)
}

#' Internal function to normalize a pairwise distance 
#'
#' @param hist Vector a numeric vector
#' @param roll Integer window used to roll the flatten histogram
#' @param doZscore Boolean should the normalized dampened signal be z-scored?
#' @param roll_smoothed.h Integer window used to flatten the histogram
#' 
#' @return a normalized vector
#' 
#' @importFrom zoo rollmean

normalizeHistogram <- function(
    hist, 
    roll = 1, 
    doZscore = TRUE, 
    roll_smoothed.h = 10
) 
{
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
