#### ---- BigWig utils ---- ####
scaleBigWigs <- function(bw.as.rle) {
    l <- IRanges::RleList(lapply(bw.as.rle, function(L) {
        S4Vectors::Rle(scale(L))
    }))
    return(l)
}
