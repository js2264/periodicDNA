#' A function to easily coerce a named list into a long data.frame
#'
#' \code{namedListToLongFormat(x)} returns a data.frame in long format, with 
#' an added 'name' column, containing the names of the input list.  
#'
#' @param x A named list vector.
#' 
#' @export
#' @import magrittr
#' 
#' @return A long data frame

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

#' A function to replace NA by something else in a vector
#'
#' \code{namedListToLongFormat(x)} returns a data.frame in long format, with 
#' an added 'name' column, containing the names of the input list.  
#'
#' @param x vector
#' @param value Replace NA by this variable
#' 
#' @export
#' 
#' @return A vector with NA replaced by value
#' 
na.replace <- function(x, value) {
    which.na <- is.na(x)
    x[which.na] <- value
    return(x)
}
