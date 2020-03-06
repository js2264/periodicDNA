#' C. elegans promoters
#'
#' Promoters annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce_proms)
#'
#' @format An object of class \code{"GRanges"}.
#'
#' @keywords datasets
#'
#' @references Serizay et al. 2020, 
#' "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#' (\href{https://doi.org/10.1101/2020.02.20.958579}{DOI})
#'
#' @source \href{https://doi.org/10.1101/2020.02.20.958579}{BiorXiv}
#'
#' @examples
#' data(ce_proms)
#' table(ce_proms$regulatory_class)
#' table(ce_proms$which.tissues)
"ce_proms"
