#' A function to compute PSDs with iterations
#' 
#' This function computes PSD values of a given k-mer of interest in 
#' a set of input sequences. It also iterates the PSD calculation process 
#' over shuffled sequences, if n_shuffling is used. 
#' 
#' @param x DNAStringSet, sequences of interest
#' @param motif character, k-mer of interest
#' @param n_shuffling integer, Number of shuffling
#' @param cores_shuffling integer, Number of cores used for shuffling
#' @param cores_computing integer, split the workload over several processors 
#' using BiocParallel
#' @param order Integer, which order to take into consideration for shuffling
#' (ushuffle python library must be installed for orders > 1)
#' @param verbose integer, Should the function be verbose? 
#' @param genome genome ID, BSgenome or DNAStringSet object 
#' (optional, if x is a GRanges)
#' @param ... Arguments passed to S3 methods
#' 
#' @return Several metrics
#' 
#' @import GenomicRanges
#' @import BiocParallel
#' @importFrom parallel mclapply
#' @import IRanges
#' @export
#' 
#' @examples
#' data(ce11_proms_seqs)
#' res <- getPeriodicityWithIterations(
#'     ce11_proms_seqs[1:10], 
#'     genome = 'BSgenome.Celegans.UCSC.ce11', 
#'     motif = 'TT', 
#'     cores_shuffling = 1
#' )
#' res$observed_PSD
#' res$shuffled_PSD

getPeriodicityWithIterations <- function(x, ...) {
    UseMethod("getPeriodicityWithIterations")
}

#' @export
#' @describeIn getPeriodicityWithIterations S3 method for DNAString

getPeriodicityWithIterations.DNAStringSet <- function(
    x, 
    motif,
    n_shuffling = 10,
    cores_shuffling = 1,
    cores_computing = 1,
    order = 1,
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
    # Shuffling sequences and iterate ---------------------------
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
    # Return res -----------------------------------------------------
    res <- list(
        observed_spectra = obs, 
        shuffled_spectra = l_shuff, 
        motif = motif
    )
    res <- getPeriodsMetrics(res)
    return(res)
}

#' @export
#' @describeIn getPeriodicityWithIterations S3 method for GRanges

getPeriodicityWithIterations.GRanges <- function(
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
        seqs <- genome[granges]
    }
    getPeriodicityWithIterations(seqs, ...)
}

#' @importFrom stats p.adjust

getPeriodsMetrics <- function(res) {
    if (all(c("PSD_withShuffling", "dists") %in% names(res))) {
        res <- res$PSD_withShuffling
    }
    freqs <- res$observed_spectra$PSD$freq
    obsPsds <- res$observed_spectra$PSD
    expPsds <- lapply(res$shuffled_spectra, '[[', 'PSD') %>% do.call(rbind, .)
    n_perms <- table(expPsds$freq)[[1]]
    if (n_perms < 100) {
        message(
            "Only ", table(expPsds$freq)[[1]], " shufflings. ", 
            "Cannot compute accurate empirical p-values. ",
            "To compute empirical p-values, ",
            "set up n_shuffling to at least 100. ",
            "Only l2FC values are returned"
        )
        pvals <- rep(as.numeric(NA), length(obsPsds$freq))
    }
    else {
        minpval <- formatC(1/(n_perms+1), format = 'e', digits = 2)
        message(
            table(expPsds$freq)[[1]], " permutations found. ", 
            "Computing empirical p-values ", 
            "(minimum = ", minpval, ")."
        )
        pvals <- unlist(lapply(obsPsds$freq, function(freq) {
            obs <- obsPsds$PSD[obsPsds$freq == freq]
            exp <- expPsds$PSD[expPsds$freq == freq]
            (sum(exp > obs) + 1)/(n_perms + 1)
        }))
        pvals <- as.numeric(formatC(pvals, format = "e", digits = 2))
    }
    df <- data.frame(
        Freq = obsPsds$freq, 
        Period = 1/obsPsds$freq,
        PSD_observed = as.numeric(
            formatC(obsPsds$PSD, format = "e", digits = 2)
        ), 
        l2FC = unlist(lapply(obsPsds$freq, function(freq) {
            obs_PSD <- obsPsds$PSD[obsPsds$freq == freq]
            l_shuff_PSD <- expPsds$PSD[expPsds$freq == freq]
            log2(obs_PSD/median(l_shuff_PSD))
        })),
        pval = pvals,
        fdr = stats::p.adjust(pvals, method = 'fdr')#,
    )
    res$periodicityMetrics <- df
    return(res)
}
