#' A function to easily coerce a named list into a long data.frame
#'
#' @param x A named list.
#' 
#' @return A long data frame
#' 
#' @import magrittr
#' @export

namedListToLongFormat <- function(x) {
    lapply(names(x), function(NAME) {
        L <- x[[NAME]]
        if (is.null(ncol(L))) {
            data.frame(value = L, name = rep(NAME, length(L)))
        } 
        else {
            data.frame(L, name = rep(NAME, nrow(L)))
        }
    }) %>% do.call(rbind, .)
}

#' A tidyverse-compatible function to replace NAs in a vector by a given value
#'
#' @param x vector
#' @param value Replace NAs by this variable
#' 
#' @return A vector with NA replaced by value
#' 
#' @export

na.replace <- function(x, value) {
    which.na <- is.na(x)
    x[which.na] <- value
    return(x)
}

#' A tidyverse-compatible function to remove NAs from a vector
#'
#' @param x vector
#' 
#' @return A vector without NAs
#' 
#' @export

na.remove <- function(x) {
    which.na <- is.na(x)
    return(x[!which.na])
}
