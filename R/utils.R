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
