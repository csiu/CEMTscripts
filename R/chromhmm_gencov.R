#' Pull out column(s) from ChromHMM overlap enrichment file
#' @param enrichment_file
#'   ChromHMM's enrichment overlap file
#' @param col Column(s) to pull out
chromhmm_getenrichmentcol <- function(enrichment_file, col) {
  dat <- readr::read_tsv(enrichment_file)
  dat <- dat[,c(colnames(dat)[1], col)]
  colnames(dat)[1] <- "state"
  dat  %>%
    filter(state != "Base")
}

# -------------------------------------------------------------------------
#' Load ChromHMM genomic coverage
#'
#' Load ChromHMM genomic coverage and find the average genomic coverage
#' if multiple enrichment files are given.
#' @param enrichment_file
#'   Enrichment overlap file containing genomic coverage.
#'   The format of this file should be TSV, have headers,
#'   and column 1 be the states.
#'   When more than 1 file is given, the average is computed
#' @param gencov_col
#'   Name of genomic coverage column in the
#'   \code{enrichment_file}.
#' @import dplyr
#' @export
chromhmm_gencov <- function(enrichment_file, gencov_col="Genome %"){
  ## Load gencovs into list
  dats <- NULL
  for (f in enrichment_file){
    dats[[f]] <-
      CEMTscripts:::chromhmm_getenrichmentcol(f, gencov_col)
  }
  ## Reduce list to single dataframe
  dats <-
    Reduce(function(...){dplyr::left_join(..., by="state")}, dats)
  ## Compute average
  dats %>%
    tidyr::gather(batch, gencov, -state) %>%
    group_by(state) %>%
    summarize(gencov = mean(gencov)) %>%
    ungroup()
}
