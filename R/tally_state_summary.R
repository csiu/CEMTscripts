#' Summarize overlaps in a given state
#'
#' In a given state, tally the number of bins which overlap
#'
#' @param input_file Format is described in \code{load_tallystate()}
tally_state <- function(input_file){
  require(dplyr)
  ## The number of bins which agree among n samples
  dd <- CEMTscripts:::load_tallystate(input_file) %>%
    group_by(samples) %>%
    summarize(nbins = sum(bins))
  setNames(dd$nbins, as.character(dd$samples))
}

#  ------------------------------------------------------------------------
#' Summarize overlaps for all states
#'
#' In a given state, tally the number of bins which overlap.
#' Assumes input files end in "bed" and
#' a state e.g. state 1's filename is named "*.U1.*", "*.E1.*" *.1.*"
#'
#' @param tally_dir directory containing the state files
#' @param verbose boolean; print messages?
#' @param cache
#'          boolean; If TRUE, results will be saved to
#'          \code{tally_dir}/tallystate_summary.txt &/or
#'          data will be loaded from this file
#' @import dplyr
#' @export
tally_state_summary <- function(tally_dir, verbose=FALSE, cache=TRUE){
  cache_file <- file.path(tally_dir, "tallystate_summary.txt")

  if (cache && file.exists(cache_file)) return(readr::read_tsv(cache_file))

  # Input files
  input_files <- CEMTscripts:::getfile_tallystate(tally_dir, state=NULL)

  # Generate summary data
  dat <- NULL
  for (s in 1L:length(input_files)){
    if (verbose) message(paste("Loading state:", s))
    input_file <-
      CEMTscripts:::getfile_tallystate(tally_dir, state=s)
    dat <-
      rbind(
        dat,
        c(state=s, CEMTscripts:::tally_state(input_file))
      )
  }

  if (cache) readr::write_tsv(tbl_df(dat), cache_file)
  tbl_df(dat)
}
