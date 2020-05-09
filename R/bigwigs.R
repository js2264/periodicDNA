#' A function to z-score a bigWig track imported by rtracklayer. 
#'
#' This function takes a RleList (usually obtained by 
#' rtracklayer::import(..., as = 'Rle')) and scales each chromosome
#' idependently by Z-score. 
#'
#' @param bigWigs a RleList
#' @return a list of scaled Rle (z-score)
#' 
#' @import S4Vectors
#' @importFrom methods is
#' @export
#' 
#' @examples
#' data(ce11_WW_10bp)
#' scaleBigWigs(ce11_WW_10bp)

scaleBigWigs <- function(bigWigs) {
    if (methods::is(bigWigs, 'RleList')) {
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
