#' Summarize overlaps in a given state
#'
#' In a given state, tally the number of bins which overlap
#'
#' @param input_file
#'
tally_state <- function(input_file){
  require(dplyr)
  ## The number of bins which agree among n samples
  dd <- readr::read_tsv(input_file,
                 col_names = c("chrom", "start", "end",
                               "samples", "label", "bins")) %>%
    group_by(samples) %>%
    summarize(nbins = sum(bins))
  setNames(dd$nbins, as.character(dd$samples))
}

#  ------------------------------------------------------------------------
#' Summarize overlaps in a given state
#'
#' In a given state, tally the number of bins which overlap.
#' Assumes input files end in "bed" and
#' a state i.e. state 1 is named ".U1." or ".1."
#'
#' @param tally_dir
#' @export
tally_state_summary <- function(tally_dir, out_file=NULL){
  require(dplyr)

  # Input files
  input_files <- list.files(path=tally_dir, pattern="*bed")

  # Generate summary data
  dat <- NULL
  for (s in 1L:length(input_files)){
    input_file <-
      file.path(tally_dir,
                grep(sprintf("[\\.U]%s\\.", s), input_files, value=T)
                )
    dat <-
      rbind(
        dat,
        c(state=s, tally_state(input_file))
      )
  }

  if (!is.null(out_file)) readr::write_tsv(dat, out_file)
  tbl_df(dat)
}
