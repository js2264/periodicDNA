#### ---- BigWig utils ---- ####

#' Core function
#'
#' @param bw.as.rle a RleList
#' 
#' @importFrom S4Vectors Rle
#' @export
#' @return a list of scaled Rle (z-score)

scaleBigWigs <- function(bw.as.rle) {
    l <- IRanges::RleList(lapply(bw.as.rle, function(L) {
        S4Vectors::Rle(scale(L))
    }))
    return(l)
}
