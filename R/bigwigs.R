#' @import S4Vectors
#' @importFrom methods is

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
