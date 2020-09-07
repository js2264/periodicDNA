#' A function to compute k-mer periodicity in sequence(s). 
#'
#' This function takes a set of sequences and a k-mer of interest, 
#' map a k-mer of interest in these sequences, computes all the 
#' pairwise distances (distogram), normalize it for distance decay, 
#' and computes the resulting power spectral density of the 
#' normalized distogram.
#' 
#' @param x a DNAString, DNAStringSet or GRanges object. 
#' @param motif a k-mer of interest
#' @param range_spectrum Numeric vector Range of the distogram to use to run 
#' the Fast Fourier Transform on (default: 1:200, i.e. all pairs of k-mers 
#' at a maximum of 200 bp from each other). 
#' @param BPPARAM split the workload over several processors using 
#' BiocParallel
#' @param roll Integer Window to smooth the distribution of pairwise distances
#' (default: 3, to discard the 3-bp periodicity of dinucleotides which 
#' can be very strong in vertebrate genomes)
#' @param verbose Boolean
#' @param sample Integer if > 0, will randomly sample this many integers
#' from the dists vector before normalization. This ensures consistency 
#' when looking at periodicity in different genomes, since different
#' genomes will have different GC percent
#' @param n_shuffling Integer, how many times should the sequences be 
#' shuffled? (default = 0)
#' @param genome genome ID, BSgenome or DNAStringSet object 
#' (optional, if x is a GRanges)
#' @param cores_shuffling integer, Number of cores used for shuffling 
#' (used if n_shuffling > 0)
#' @param cores_computing integer, split the workload over several processors 
#' using BiocParallel (used if n_shuffling > 0)
#' @param order Integer, which order to take into consideration for shuffling
#' (ushuffle python library must be installed for orders > 1) 
#' (used if n_shuffling > 0)
#' @param ... Arguments passed to S3 methods
#' 
#' @return A list containing the results of getPeriodicity function.  
#' \itemize{
#'     \item The dists vector is the raw vector of all distances between 
#'     any possible k-mer. 
#'     \item The hist data.frame is the distribution of distances
#'     over range_spectrum. 
#'     \item The normalized_hist is the raw hist, 
#'     normalized for decay over increasing distances. 
#'     \item The spectra object is the output of 
#'     the FFT applied over normalized_hist. 
#'     \item The PSD data frame is the power 
#'     spectral density scores over given  frequencies. 
#'     \item The motif object is the k-mer being analysed.
#'     \item The final periodicity metrics computed by getPeriodicity()
#' }
#' If getPeriodicity() is ran with n_shuffling > 0, the resulting 
#' list also contains PSD values computed when iterating through shuffled 
#' sequences.
#' 
#' @import BiocParallel
#' @import Biostrings
#' @import IRanges
#' @import magrittr
#' @importFrom stats spectrum
#' @importFrom stats dist
#' @importFrom stats setNames
#' @importFrom methods is
#' @importFrom BSgenome getBSgenome
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' periodicity_result <- getPeriodicity(
#'     ce11_proms_seqs[1:100],
#'     motif = 'TT'
#' )
#' head(periodicity_result$PSD)
#' plotPeriodicityResults(periodicity_result)
#' #
#' data(ce11_TSSs)
#' periodicity_result <- getPeriodicity(
#'     ce11_TSSs[['Ubiq.']][1:10],
#'     motif = 'TT',
#'     genome = 'BSgenome.Celegans.UCSC.ce11'
#' )
#' head(periodicity_result$PSD)
#' plotPeriodicityResults(periodicity_result)
#' #
#' data(ce11_TSSs)
#' periodicity_result <- getPeriodicity(
#'     ce11_TSSs[['Ubiq.']][1:10],
#'     motif = 'TT',
#'     genome = 'BSgenome.Celegans.UCSC.ce11',
#'     n_shuffling = 10
#' )
#' head(periodicity_result$PSD)
#' plotPeriodicityResults(periodicity_result)

getPeriodicity <- function(x, motif, ...) {
    UseMethod("getPeriodicity")
}

#' @export
#' @describeIn getPeriodicity S3 method for DNAStringSet

getPeriodicity.DNAStringSet <- function(
    x,
    motif,
    range_spectrum = seq(1, 200),
    BPPARAM = setUpBPPARAM(1),
    roll = 3,
    verbose = TRUE,
    sample = 0,
    n_shuffling = 0,
    cores_shuffling = 1,
    cores_computing = 1,
    order = 1,
    ...
)
{
    seqs <- x
    if (n_shuffling == 0) {
        # Get pairwise distances ----------------------------------------------
        if (verbose & n_shuffling == 0) message("- Mapping k-mers.")
        dists <- BiocParallel::bplapply(seq_len(length(seqs)), function(k) {
            Biostrings::vmatchPattern(
                motif, 
                seqs[k], 
                max.mismatch = 0, 
                fixed = FALSE
            )[[1]] %>% 
                IRanges::start() %>% 
                stats::dist() %>% 
                c()
        }, BPPARAM = BPPARAM) %>% unlist()
        if (length(dists) == 0) {
            stop(
                'ERROR: No ', motif, '..', motif, ' pair was found in the ',
                'input sequences. Please input other sequences or search ', 
                'for another motif.\n',
                '  ABORTING NOW.'
            )
        }
        max_dist <- max(dists)
        if (max(dists) < tail(range_spectrum, 1)) {
            stop(
                'ERROR: range_spectrum (', 
                head(range_spectrum, 1), ':', tail(range_spectrum, 1), 
                ') is wider than the maximum distance between ', 
                'pairs of k-mers (', max(dists), 
                '). Try shortening range_spectrum to a smaller range.\n',
                '  ABORTING NOW.'
            )
        }
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
        if (verbose & n_shuffling == 0) 
            message("- ", length(dists), " pairwise distances measured.")
        if (sample < length(dists) & sample > 0) {
            if (verbose) 
                message("- Subsampling ", sample, " pairwise distances.")
            dists <- sample(dists, sample)
        }
        if (verbose & n_shuffling == 0) 
            message("- Calculating pairwise distance distribution.")
        hist <- hist(
            dists, breaks = seq(1, max_dist+1, 1), plot = FALSE
        )$counts
        hist <- zoo::rollmean(
            hist, k = roll, na.pad = TRUE, align = 'center'
        )
        if (verbose & n_shuffling == 0) 
            message("- Normalizing distogram vector.")
        if (length(hist) > 10) {
            norm_hist <- normalizeHistogram(hist, roll = 1)
        } 
        else {
            norm_hist <- hist
        }
        # Fourier -------------------------------------------------------------
        if (verbose) message(
            "- Applying Fast Fourier Transform to the normalized distogram."
        )
        spectra <- do.call(
            stats::spectrum, 
            list(norm_hist[range_spectrum], plot = FALSE)
        )
        PSD <- data.frame(
            freq = spectra$freq, 
            period = 1/spectra$freq, 
            PSD = spectra$spec
        )
        l <- list(
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
            motif = motif, 
            periodicityMetrics = stats::setNames(
                PSD, c('Freq', 'Period', 'PSD')
            )
        )
    }
    else {
        res <- getPeriodicityWithIterations(
            x = seqs, 
            motif = motif,
            range_spectrum = range_spectrum,
            roll = roll,
            verbose = verbose,
            sample = sample,
            n_shuffling = n_shuffling, 
            cores_shuffling = cores_shuffling,
            cores_computing = cores_computing,
            order = order,
            ...
        )
        l <- list(
            dists = res$observed_spectra$dists, 
            hist = res$observed_spectra$hist,
            normalized_hist = res$observed_spectra$normalized_hist,
            spectra = res$observed_spectra$spectra, 
            PSD = res$observed_spectra$PSD,
            PSD_withShuffling = res, 
            motif = res$observed_spectra$motif, 
            periodicityMetrics = res$periodicityMetrics
        )
    }
    # Return results ----------------------------------------------------------
    return(l)
}

#' @export
#' @describeIn getPeriodicity S3 method for GRanges

getPeriodicity.GRanges <- function(
    x,
    motif, 
    genome = 'BSgenome.Celegans.UCSC.ce11',
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
            ) | genome %in% BSgenome::installed.genomes()) {
                genome <- BSgenome::getBSgenome(genome)
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
    getPeriodicity(seqs, motif, ...)
}

#' @export
#' @describeIn getPeriodicity S3 method for DNAString

getPeriodicity.DNAString <- function(
    x,
    motif, 
    ...
)
{
    seq <- Biostrings::DNAStringSet(x)
    getPeriodicity(seq, motif, ...)
}

#' @importFrom zoo rollmean

normalizeHistogram <- function(
    hist, 
    roll = 1, 
    roll_smoothed.h = 10
) 
{
    h <- hist
    # Normalize to total number of pairwise distances
    h <- h / sum(h, na.rm = TRUE) 
    # Substract smoothed distribution
    smoothed.h <- c(
        zoo::rollmean(h, k = roll_smoothed.h), 
        rep(0, roll_smoothed.h-1)
    ) 
    norm.h <- h - smoothed.h
    # Smooth normalized distribution 
    if (roll > 1) {
        norm.h <- zoo::rollmean(
            norm.h, k = roll, na.pad = TRUE, align = 'center'
        ) 
    }
    norm.h[is.na(norm.h)] <- 0
    return(norm.h)
}
