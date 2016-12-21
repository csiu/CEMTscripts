#' Load ChromHMM enrichment file
#' @param enrichment_file
#'          Assumes the first column of the TSV is chromatin state
#' @export
chromhmm_loadenrichment <- function(enrichment_file) {
  dat <- readr::read_tsv(enrichment_file, progress = FALSE)
  colnames(dat)[1] <- "state"
  dat  %>%
    filter(state != "Base")
}

# -------------------------------------------------------------------------
#' Pull out column(s) from ChromHMM overlap enrichment file
#' @param enrichment_file
#'   ChromHMM's enrichment overlap file
#' @param col Column(s) to pull out
chromhmm_getenrichmentcol <- function(enrichment_file, col) {
  dat <- CEMTscripts:::chromhmm_loadenrichment(enrichment_file)
  dat[,c(colnames(dat)[1], col)]
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
#' @param isstateorder boolean; treat state as factor and order states?
#' @import dplyr
#' @export
chromhmm_gencov <- function(enrichment_file, gencov_col="Genome %",
                            isstateorder=FALSE){
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
  dats <- dats %>%
    tidyr::gather(batch, gencov, -state) %>%
    group_by(state) %>%
    summarize(gencov = mean(gencov)) %>%
    ungroup()

  if (isstateorder){
    dats$state<-factor(dats$state, levels=1:nrow(dats))
  }
  arrange(dats, state)
}
