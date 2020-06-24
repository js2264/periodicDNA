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
#' @param order Integer, which order to take into consideration for shuffling
#' (currently only 1st order is available)
#' @param n_shuffling integer, Number of shuffling
#' @param cores_shuffling integer, Number of cores used for shuffling
#' @param cores_computing integer, split the workload over several processors 
#' using BiocParallel
#' @param verbose integer, Should the function be verbose? 
#' @param ... Additional arguments
#' @return Several metrics including FPI, observed PSD, etc...
#' 
#' @import GenomicRanges
#' @import BiocParallel
#' @import IRanges
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' BiocParallel::register(BiocParallel::SnowParam(workers = 1))
#' fpi <- getFPI(
#'     ce11_proms_seqs[1:10], 
#'     motif = 'TT'
#' )
#' fpi$FPI
#' fpi$observed_PSD
#' fpi$shuffled_PSD
#' plotFPI(fpi)

getFPI.DNAStringSet <- function(
    x, 
    motif,
    period = 10, 
    order = 1,
    n_shuffling = 10,
    cores_shuffling = 1,
    cores_computing = 1,
    verbose = 1,
    ...
)
{
    seqs <- x
    
    # Calculating observed PSD ---------------------------------------
    if (verbose) message('- Calculating observed PSD')
    obs <- getPeriodicity(
        seqs, 
        motif = motif,
        skip_shuffling = TRUE,
        verbose = verbose,
        BPPARAM = SnowParam(workers = cores_computing),
        ...
    )
    obs_PSD <- obs$PSD$PSD[which.min(abs(1/obs$PSD$freq - period))]
    if (verbose) message('>> Measured PSD @ ', period, 'bp is: ', obs_PSD)
    # Shuffling sequences and re-computing ---------------------------
    l_shuff <- BiocParallel::bplapply(
        BPPARAM = SnowParam(workers = cores_shuffling), 
        seq_len(n_shuffling), 
        function(k) {
            if (verbose) message('- Shuffling ', k, '/', n_shuffling)
            shuff_seqs <- shuffleSeq(seqs, order)
            shuff <- getPeriodicity(
                shuff_seqs, 
                motif,
                skip_shuffling = TRUE,
                verbose = 0,
                BPPARAM = SnowParam(workers = cores_computing),
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
    median_fold <- obs_PSD / median(l_shuff_PSD)
    if (verbose) message('>> Calculated FPI @ ', period, 'bp is: ', fpi)
    # Return FPI -----------------------------------------------------
    res <- list(
        FPI = fpi, 
        median_fold = median_fold,
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
        seqs <- genome[granges]
    }
    getFPI(seqs, ...)
}

#' Internal function
#'
#' @param fpi result of getFPI() function
#' @param threshold Float, p-value used as significance threshold
#' @return periods from the FPI list at which PSD are statistically higher 
#' than those from shuffled sequences
#' @importFrom stats t.test

significantPeriods <- function(fpi, threshold = 0.0001) {
    obsPsds <- fpi$observed_spectra$PSD
    expPsds <- lapply(fpi$shuffled_spectra, '[[', 'PSD') %>% do.call(rbind, .)
    df <- data.frame(
        freq = obsPsds$freq, 
        period = 1/obsPsds$freq,
        pval = unlist(lapply(obsPsds$freq, function(freq) {
            stats::t.test(
                obsPsds$PSD[obsPsds$freq == freq], 
                expPsds$PSD[expPsds$freq == freq],
                var.equal = TRUE
            )$p.value
        }))
    )
    periods <- df$period[which(df$pval < threshold)]
    return(periods)
}
