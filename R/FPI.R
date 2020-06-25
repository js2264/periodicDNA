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
#' (ushuffle python library must be installed for orders > 1)
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
#' @importFrom parallel mclapply
#' @import IRanges
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' BiocParallel::register(setUpBPPARAM(1))
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
        verbose = verbose,
        BPPARAM = setUpBPPARAM(cores_computing),
        ...
    )
    obs_PSD <- obs$PSD$PSD[which.min(abs(1/obs$PSD$freq - period))]
    if (verbose) message('>> Measured PSD @ ', period, 'bp is: ', obs_PSD)
    # Shuffling sequences and re-computing ---------------------------
    l_shuff <- BiocParallel::bplapply(
        BPPARAM = setUpBPPARAM(cores_shuffling), 
        seq_len(n_shuffling), 
        function(k) {
            if (verbose) message('- Shuffling ', k, '/', n_shuffling)
            shuff_seqs <- shuffleSeq(seqs, order)
            shuff <- getPeriodicity(
                shuff_seqs, 
                motif,
                verbose = 0,
                BPPARAM = setUpBPPARAM(cores_computing),
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

#' A function used to compute FPIs and associated p-values
#'
#' @param fpi result of getFPI() function or getPeriodicity() function with
#' n_shuffling > 1
#' @return A table showing FPI and p-value for each observed PSD value 
#' 
#' @importFrom stats t.test
#' @importFrom stats p.adjust
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
#' head(getSignificantPeriods(fpi))

getSignificantPeriods <- function(fpi) {
    if (all(c("FPI", "dists") %in% names(fpi))) {
        fpi <- fpi$FPI
    }
    freqs <- fpi$observed_spectra$PSD$freq
    obsPsds <- fpi$observed_spectra$PSD
    expPsds <- lapply(fpi$shuffled_spectra, '[[', 'PSD') %>% do.call(rbind, .)
    df <- data.frame(
        Freq = obsPsds$freq, 
        Period = 1/obsPsds$freq,
        ObservedPSD = formatC(obsPsds$PSD, format = "e", digits = 2),
        FPI = unlist(lapply(obsPsds$freq, function(freq) {
            obs_PSD <- obsPsds$PSD[obsPsds$freq == freq]
            l_shuff_PSD <- expPsds$PSD[expPsds$freq == freq]
            (obs_PSD - median(l_shuff_PSD)) / median(l_shuff_PSD)
        })), 
        pval = formatC(
            unlist(lapply(obsPsds$freq, function(freq) {
                stats::p.adjust(stats::t.test(
                    obsPsds$PSD[obsPsds$freq == freq], 
                    expPsds$PSD[expPsds$freq == freq],
                    var.equal = TRUE
                )$p.value, "BH")
            })), 
             format = "e", digits = 2
         )
    )
    return(df)
}
