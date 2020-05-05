#' ce_proms
#'
#' Promoters annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce_proms)
#'
#' @format GRanges
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

#' proms
#'
#' Sample of promoters annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(proms)
#'
#' @format GRanges
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
#' data(proms)
#' table(proms$regulatory_class)
#' table(proms$which.tissues)
"proms"

#' ce_proms_seqs
#'
#' Sample of sequences of promoters annotated in C. elegans (ce11) according to
#' Serizay et al. 2020, "Tissue-specific profiling reveals distinctive 
#' regulatory architectures for ubiquitous, germline and somatic genes", 
#' BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce_proms_seqs)
#'
#' @format DNAStringSet
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
#' data(ce_proms_seqs)
#' head(ce_proms_seqs)
"ce_proms_seqs"

