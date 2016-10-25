#' Summarize overlaps in a given state
#'
#' In a given state, tally the number of bins which overlap
#'
#' @import readr
#' @param input_file TSV file with 6 columns + no header:
#'                   chrom, start, stop,
#'                   number of samples, sample representation
#'                   e.g. 0010 (where 0=absent; 1=present),
#'                   number of bins in region i.e. (stop-start)/200bp.
tally_state <- function(input_file){
  require(dplyr)
  ## The number of bins which agree among n samples
  dd <- read_tsv(input_file,
                 col_names = c("chrom", "start", "end",
                               "samples", "label", "bins"),
                 col_types = cols(
                   chrom = col_character(),
                   start = col_integer(),
                   end = col_integer(),
                   samples = col_integer(),
                   label = col_character(),
                   bins = col_integer()
                 )) %>%
    group_by(samples) %>%
    summarize(nbins = sum(bins))
  setNames(dd$nbins, as.character(dd$samples))
}

#  ------------------------------------------------------------------------
#' Summarize overlaps for all states
#'
#' In a given state, tally the number of bins which overlap.
#' Assumes input files end in "bed" and
#' a state e.g. state 1's filename is named "*.U1.*" or "*.1.*"
#'
#' @param tally_dir directory containing the state files
#' @param out_file when the output filename is specified,
#'                 the output will be saved to this file
#' @param verbose boolean; print messages?
#' @export
tally_state_summary <- function(tally_dir, out_file=NULL, verbose=FALSE){
  require(dplyr)

  # Input files
  input_files <- list.files(path=tally_dir, pattern="*bed")

  # Generate summary data
  dat <- NULL
  for (s in 1L:length(input_files)){
    if (verbose) message(paste("Loading state:", s))
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
