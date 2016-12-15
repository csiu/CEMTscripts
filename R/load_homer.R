#' Load HOMER knownResults.txt
#' @param knownresults_file path to the "knownResults.txt" file
#' @param returncounts boolean; If TRUE, will return
#'                     vector of target & background sequences
#' @format Outputs a dataframe with the following columns:
#'         \code{tf, dnabindingdomain, origin, source, consensus,
#'         pval, pvalLOG, qvalBEN, tagetN, tagetP, bgN, bgP}
#' @source \url{http://homer.salk.edu/homer/motif/motifDatabase.html}
#' @export
load_homer_knownresults <- function(knownresults_file, returncounts=FALSE){
  motifs <- readr::read_tsv(knownresults_file, progress=FALSE)

  ## Adjust column names
  original_columns <- colnames(motifs)
  if (returncounts) {
    ## Number of sequences
    seqcounts <- c(
      "targets" = as.integer(sub(".*?(\\d+).*", "\\1", original_columns[6])),
      "background" = as.integer(sub(".*?(\\d+).*", "\\1", original_columns[8]))
    )
    return(seqcounts)
  }

  ## Tidy column names & split motif name
  colnames(motifs) <-
    c("name", "consensus", "pval", "pvalLOG", "qvalBEN",
      "tagetN", "targetP", "bgN", "bgP")
  motifs %>%
    tidyr::separate("name", into=c("tf", "origin", "source"), sep="/") %>%
    tidyr::separate("tf", into=c("tf", "dnabindingdomain"), sep="\\(") %>%
    dplyr::mutate(dnabindingdomain = sub("\\)$", "", dnabindingdomain))
}
