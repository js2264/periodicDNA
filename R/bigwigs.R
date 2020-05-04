#' A function to z-score a bigWig track imported by rtracklayer. 
#'
#' @param bigWigs a RleList, usually obtained by rtracklayer::import() 
#' 
#' @return a list of scaled Rle (z-score)
#' 
#' @import S4Vectors
#' @importFrom methods is
#' @export

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
