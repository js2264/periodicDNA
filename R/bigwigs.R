#' A function to z-score a bigWig track imported by rtracklayer. 
#'
#' @param bigWigs a RleList, usually obtained by rtracklayer::import() 
#' 
#' @return a list of scaled Rle (z-score)
#' 
#' @import S4Vectors
#' @export

scaleBigWigs <- function(bigWigs) {
    if (class(bigWigs) == 'SimpleRleList' | class(bigWigs) == 'CompressedRleList') {
        bigWig <- bigWigs
        ucov <- unlist(bigWig)
        mi <- mean(ucov)
        mu <- S4Vectors::sd(ucov)
        zsc <- (bigWig-mi)/mu
        return(zsc)
    } else {
        l <- lapply(bigWigs, function(bigWig) {
            ucov <- unlist(bigWig)
            mi <- mean(ucov)
            mu <- S4Vectors::sd(ucov)
            zsc <- (bigWig-mi)/mu
            return(zsc)
        })
        return(l)
    }
}
