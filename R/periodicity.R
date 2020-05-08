#' Function to compute the overall periodicity of a 
#'     motif over sequence(s). 
#'
#' @param x a DNAStringSet, or a GRanges
#' @param ... additional parameters
#' 
#' @return List a list containing the results of getPeriodicity function. 
#'     The dists vector is the raw vector of all distances between any 
#'     possible dinucleotide. The hist data.frame is the distribution of 
#'     distances over range_spectrum. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spectral density scores over given 
#'     frequencies. The motif object is the dinucleotide being analysed.
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
#' @param range_spectrum Numeric vector The distances between nucleotides
#'     to take into consideration when performing Fast Fourier Transform.
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
#'     distances over range_spectrum. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spectral density scores over given 
#'     frequencies. The motif object is the dinucleotide being analysed.
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
    range_spectrum = seq(1, 200),
    cores = 1, 
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
        Biostrings::vmatchPattern(
            motif, 
            seqs[k], 
            max.mismatch = 0, 
            fixed = FALSE
        )[[1]] %>% 
            IRanges::start() %>% 
            dist() %>% 
            c()
    }, mc.cores = cores) %>% unlist()
    max_dist <- max(dists)
    if (length(dists) < 10) {
        if (verbose) message("- Only ", length(dists), 
        " pairs of k-mers found. Returning null results")
        return(list(
            dists = dists, 
            hist = 0, 
            normalized_hist = 0,
            spectra = 0, 
            PSD = 0,
            motif = motif
        ))
    }
    if (verbose) message("- ", length(dists), " pairwise distances measured.")
    if (sample < length(dists) & sample > 0) {
        if (verbose) message("- Subsampling ", sample, " k-mers.")
        dists <- dists[sample(seq_len(dists), sample)]
    }
    if (verbose) message("- Calculating pairwise distance distribution.")
    hist <- hist(dists, breaks = seq(1, max_dist+1, 1), plot = FALSE)$counts
    hist <- zoo::rollmean(
        hist, k = roll, na.pad = TRUE, align = 'center'
    )
    if (verbose) message("- Normalizing histogram vector.")
    if (length(hist) > 10) {
        norm_hist <- normalizeHistogram(hist, roll = 1, doZscore)
    } 
    else {
        norm_hist <- hist
    }
    # Fourier -----------------------------------------------------------------
    if (verbose) message(
        "- Applying Fast Fourier Transform to the vector of distances."
    )
    if (max_dist < tail(range_spectrum, 1)) {
        if (verbose) message(
            'Range (', 
            head(range_spectrum, 1), ':', tail(range_spectrum, 1), 
            ') is wider than any range in the current 
            distances vector. Shortening the range from ', 
            head(range_spectrum, 1), ' to ', max_dist, '...'
        )
        range_spectrum <- head(range_spectrum, 1):max_dist
    }
    spectra <- do.call(
        stats::spectrum, 
        list(norm_hist[range_spectrum], plot = FALSE)
    )
    PSD <- data.frame(
        freq = spectra$freq, 
        period = 1/spectra$freq, 
        PSD = spectra$spec
    )
    # Return all results if skip_shuffling ------------------------------------
    if (skip_shuffling) {
        return(list(
            dists = dists, 
            hist = data.frame(
                distance = seq(1, max_dist, 1), 
                counts = hist
            ), 
            normalized_hist = data.frame(
                distance = seq(1, max_dist, 1), 
                norm_counts = norm_hist
            ),
            spectra = spectra, 
            PSD = PSD,
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
        max_dists_shuffled <- max(dists_shuffled)
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
            breaks = seq(1, max_dists_shuffled+1, 1), 
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
        if (max_dists_shuffled < tail(range_spectrum, 1)) {
            if (verbose) message(
                'Range (', 
                head(range_spectrum, 1), ':', tail(range_spectrum, 1), 
                ') is wider than any range in the current distances vector. 
                Shortening the range from ', 
                head(range_spectrum, 1), ' to ', max_dists_shuffled, '...'
            )
            a <- head(range_spectrum, 1)
            range_spectrum <- a:max_dists_shuffled
        }
        spectra_shuffled <- norm_hist_shuffled %>%
                '['(range_spectrum) %>% 
                stats::spectrum(plot = FALSE)
        PSD_shuffled <- data.frame(
            freq = spectra_shuffled$freq, 
            PSD = spectra_shuffled$spec
        )
    }
    return(list(
        dists = dists, 
        hist = data.frame(
            distance = seq(1, max_dist, 1), 
            counts = hist
        ), 
        normalized_hist = data.frame(
            distance = seq(1, max_dist, 1), 
            norm_counts = norm_hist
        ),
        spectra = spectra, 
        PSD = PSD,
        dists_shuffled = dists_shuffled, 
        hist_shuffled = data.frame(
            distance = seq(1, max_dists_shuffled, 1), 
            counts = hist_shuffled
        ), 
        normalized_hist_shuffled = data.frame(
            distance = seq(1, max_dists_shuffled, 1), 
            norm_counts = norm_hist_shuffled
        ),
        spectra_shuffled = spectra_shuffled, 
        PSD_shuffled = PSD_shuffled,
        motif = motif
    ))
}

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
#'     distances over range_spectrum. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spectral density scores over given 
#'     frequencies. The motif object is the dinucleotide being analysed.
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
    
    if (!is.null(granges$seq)) {
        seqs <- granges$seq
    }
    else {
        if (methods::is(genome, 'character')) {
            if (genome %in% c(
                'sacCer3', 'ce11', 'dm6', 'mm10', 'hg38', 'danRer10'
            )) {
                genome <- char2BSgenome(genome)
            }
            else {
                return(stop(
                    'Only sacCer3, ce11, dm6, mm10, hg38 
                    and danRer10 are supported'
                ))
            }
        }
        if (methods::is(genome, 'BSgenome')) {
            genome <- Biostrings::getSeq(genome)
        }
        seqs <- withSeq(granges, genome)$seq
    }
    getPeriodicity(seqs, ...)
}

#' Core function
#' 
#' @param x a DNAString
#' @param ... other parameters forwarded to getPeriodicity.DNAStringSet()
#'
#' @return List a list containing the results of getPeriodicity function. 
#'     The dists vector is the raw vector of all distances between any 
#'     possible dinucleotide. The hist data.frame is the distribution of 
#'     distances over range_spectrum. The normalized_hist is the raw 
#'     hist, normalized for  decay over increasing distances. The spectra 
#'     object is the output of the FFT applied over normalized_hist. 
#'     The PSD data frame is the power spectral density scores over given 
#'     frequencies. The motif object is the dinucleotide being analysed.
#' 
#' @importFrom Biostrings DNAStringSet
#' 
#' @export

getPeriodicity.DNAString <- function(
    x,
    ...
)
{
    seq <- Biostrings::DNAStringSet(x)
    getPeriodicity(seq, ...)
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
    h <- h / sum(h, na.rm = TRUE) # Normalize to total number of pairwise distances
    # Substract smoothed distribution
    smoothed.h <- c(
        zoo::rollmean(h, k = roll_smoothed.h), 
        rep(0, roll_smoothed.h-1)
    ) 
    norm.h <- h - smoothed.h
    # Z-score
    if (doZscore) {
        norm.h <- scale(norm.h) 
    } 
    # Smooth normalized distribution 
    if (roll > 1) {
        norm.h <- zoo::rollmean(
            norm.h, k = roll, na.pad = TRUE, align = 'center'
        ) 
    }
    norm.h[is.na(norm.h)] <- 0
    return(norm.h)
}
