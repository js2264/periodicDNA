generatePeriodicityTrack <- function(
    ce_seq, 
    granges, 
    MOTIF = 'WW', 
    EXTEND.GRANGES = 1000, 
    GENOME.WINDOW.SIZE = 250, 
    GENOME.WINDOW.SLIDING = 2, 
    BIN.WINDOW.SIZE = 100, 
    BIN.WINDOW.SLIDING = 10, 
    RANGE.FOR.SPECTRUM = 1:50, 
    FREQ = 1/10, 
    PROCS = 10, 
    bw.file = NULL
    ) {
    
    granges.extended <- GenomicRanges::resize(granges, EXTEND.GRANGES, fix = 'center')
    genome.partionned <- partitionGenome(ce_seq, GENOME.WINDOW.SIZE, GENOME.WINDOW.SLIDING)
    genome.partionned.filtered <- IRanges::subsetByOverlaps(genome.partionned, granges.extended)
    chunks <- getChunks(genome.partionned.filtered, split.nb = PROCS)
    
    message('
    The genome has been divided in ', GENOME.WINDOW.SIZE, '-bp bins
    with a sliding window of ', GENOME.WINDOW.SLIDING, '-bp.
    
    ', length(genome.partionned.filtered), ' windows overlap the input loci (extended to 
    ', EXTEND.GRANGES, '-bp).
    
    The mapping will be split among ', PROCS, ' jobs.
    
    Each job will analyse ', lengths(chunks)[1], ' bins.
    
    Now starting...
    ')
    
    # Process each chunk separately
    list.results <- parallel::mclapply(
        chunks, 
        chunkWrapper,
        ce_seq, 
        MOTIF, 
        GENOME.WINDOW.SLIDING,
        BIN.WINDOW.SIZE, 
        BIN.WINDOW.SLIDING,
        RANGE.FOR.SPECTRUM, 
        FREQ,
        mc.cores = length(chunks)
    )
    # Recover jobs if it crashed
    #list.results <- recoverTmp()
    # Merge all the chunks
    motif.granges <- unlistResults(list.results)
    # Generate final track
    if (is.null(bw.file)) bw.file = paste0(MOTIF, '-periodicity_g-', GENOME.WINDOW.SIZE, '^', GENOME.WINDOW.SLIDING, '_b-', BIN.WINDOW.SIZE, '^', BIN.WINDOW.SLIDING ,'.bw')
    exportBigWigTrack(motif.granges, bw.file)
    # Remove tmp files
    cleanUpDirectory()
}
partitionGenome <- function(ce_seq, GENOME.WINDOW.SIZE, GENOME.WINDOW.SLIDING) {
    genome.Seqinfo <- GenomeInfoDb::Seqinfo(seqnames = names(ce_seq), seqlengths = lengths(ce_seq), isCircular = rep(F, length(ce_seq)), genome = 'custom')
    genome.granges <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(genome.Seqinfo), IRanges::IRanges(start = rep(1, length(genome.Seqinfo)), width = GenomeInfoDb::seqlengths(genome.Seqinfo)))
    GenomeInfoDb::seqinfo(genome.granges) <- genome.Seqinfo
    granges.partionned <- genome.granges[0,]
    for (CHR in GenomicRanges::seqnames(genome.Seqinfo)) {
        seqnames <- rep(CHR, GenomeInfoDb::seqlengths(genome.Seqinfo)[GenomicRanges::seqnames(genome.Seqinfo) == CHR])
        seqnames <- seqnames[seq(1, length(seqnames), GENOME.WINDOW.SLIDING)]
        seqnames <- factor(seqnames, levels = levels(GenomicRanges::seqnames(genome.granges)))
        seqstarts <- seq(1, GenomeInfoDb::seqlengths(genome.Seqinfo)[GenomicRanges::seqnames(genome.Seqinfo) == CHR], 1)
        seqstarts <- seqstarts[seq(1, length(seqstarts), GENOME.WINDOW.SLIDING)]
        granges.partionned <- c(granges.partionned, GenomicRanges::GRanges(seqnames, IRanges::IRanges(start = seqstarts, width = GENOME.WINDOW.SIZE)))
    }
    
    granges.partionned <- GenomicRanges::trim(granges.partionned)
    granges.partionned <- granges.partionned[GenomicRanges::width(granges.partionned) == GENOME.WINDOW.SIZE]
    
    return(granges.partionned)
}
getChunks <- function(ce11.w250s2.filtered, split.nb = nb.procs) {
    chunks <- as.numeric(cut(1:length(ce11.w250s2.filtered), breaks = seq(1, length(ce11.w250s2.filtered), length.out = split.nb + 1), include.lowest = T))
    list.granges <- list()
    for (K in 1:split.nb) {
        list.granges[[K]] <- ce11.w250s2.filtered[chunks == K]
    }
    return(list.granges)
}
chunkWrapper <- function(chunk, ce_seq, MOTIF, GENOME.WINDOW.SLIDING, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING, RANGE.FOR.SPECTRUM, FREQ) {
    res <- chunk
    res$score <- 0
    res <- GenomicRanges::resize(res, width = GENOME.WINDOW.SLIDING, fix = 'center')
    intervals <- seq(1, length(chunk), length.out = 100)
    for (K in 1:length(chunk)) {
        res$score[K] <- binWrapper(chunk[K], ce_seq, MOTIF, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING, RANGE.FOR.SPECTRUM, FREQ)
        if (K %in% intervals) {
            message('>> ', K, ' bins computed /// TIME: ', Sys.time())
            rtracklayer::export.bw(res[res$score > 0], paste0("tmp.", K, ".", Sys.getpid(), ".bw"))
        }
    }
    rtracklayer::export.bw(res, paste0("tmp.FULL.", Sys.getpid(), ".bw"))
    return(res)
}
binWrapper <- function(grange, ce_seq, MOTIF, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING, RANGE.FOR.SPECTRUM, FREQ) {
    bins <- withSeq(partitionBin(grange, ce_seq, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING), ce_seq)
    seqs <- bins$seq
    dists <- unlist(sapply(1:length(seqs), function(k) {
        Biostrings::vmatchPattern(MOTIF, seqs[k], max.mismatch = 0, fixed = F) %>% 
        GenomicRanges::start() %>% 
        '[['(1) %>% 
        dist() %>% 
        c()
    }))
    hist <- hist(dists, breaks = seq(0, 1000, 1), plot = F)
    s <- spectrum(hist$counts[1:50], plot = F)
    spec <- round(s$spec[which.min(abs(s$freq - FREQ))], 4)
    return(spec)
}
partitionBin <- function(grange, ce_seq, BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING) {
    genome.Seqinfo <- GenomeInfoDb::Seqinfo(seqnames = names(ce_seq), seqlengths = lengths(ce_seq), isCircular = rep(F, length(ce_seq)), genome = 'custom')
    genome.granges <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(genome.Seqinfo), IRanges::IRanges(start = rep(1, length(genome.Seqinfo)), width = GenomeInfoDb::seqlengths(genome.Seqinfo)))
    GenomeInfoDb::seqinfo(genome.granges) <- genome.Seqinfo
    starts <- seq(IRanges::start(grange), IRanges::start(grange) + GenomicRanges::width(grange) - BIN.WINDOW.SIZE, BIN.WINDOW.SLIDING)
    granges <- GenomicRanges::GRanges(seqnames = rep(GenomicRanges::seqnames(grange), length(starts)), IRanges::IRanges(start = starts, width = BIN.WINDOW.SIZE))
    return(granges)
}
unlistResults <- function(list.results) {
    granges <- list.results[[1]][0]
    for (K in 1:length((list.results))) {
        granges <- c(granges, list.results[[K]])
    }
    return(granges)
}
exportBigWigTrack <- function(granges, bw.file) {
    rtracklayer::export.bw(granges, bw.file)
}
cleanUpDirectory <- function() {
    list.files(pattern = '.*tmp.*') %>% file.remove()
    return(TRUE)
}