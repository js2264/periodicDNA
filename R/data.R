#' ce11_all_REs
#'
#' Regulatory elements annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce11_all_REs)
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
#' data(ce11_all_REs)
#' table(ce11_all_REs$regulatory_class)
#' table(ce11_all_REs$which.tissues)
"ce11_all_REs"

#' ce11_proms
#'
#' Promoters annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce11_proms)
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
#' data(ce11_proms)
#' table(ce11_proms$which.tissues)
"ce11_proms"

#' ce11_TSSs
#'
#' Coordinates of promoter TSSs annotated in C. elegans (ce11) used in Serizay 
#' et al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce11_TSSs)
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
#' data(ce11_TSSs)
#' lengths(ce11_TSSs)
#' ce11_TSSs[[1]]
"ce11_TSSs"

#' ce11_proms_seqs
#'
#' Sample of sequences of promoters annotated in C. elegans (ce11) according to
#' Serizay et al. 2020, "Tissue-specific profiling reveals distinctive 
#' regulatory architectures for ubiquitous, germline and somatic genes", 
#' BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce11_proms_seqs)
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
#' data(ce11_proms_seqs)
#' head(ce11_proms_seqs)
"ce11_proms_seqs"

