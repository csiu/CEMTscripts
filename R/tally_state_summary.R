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

#  ------------------------------------------------------------------------
#' Visualize overlaps for all states
#'
#' Plot output of \code{tally_state_summary()} into a barchart
#'
#' @param input output of \code{tally_state_summary()}
#'              i.e. dataframe where column 1 is state, and the remaining
#'              column are the number of bins which overlap 1, 2, ... samples.
#' @param statepalette character list where each element is a color and
#'                     the name of each element is a \code{state} level.
#' @param facet_wrap boolean
#' @import dplyr
#' @export
visualize_state_summary <- function(input, statepalette=NULL, facet_wrap=FALSE){
  require(ggplot2)
  dat <-
    input %>%
    dplyr::mutate(state = sub("^U", "", state) %>% as.integer()) %>%
    tidyr::gather(samples, bins, -state) %>%
    dplyr::mutate(samples = as.integer(samples))

  # Plot
  plt <-
    dat %>%
    ggplot(aes(x = samples, y = bins)) +
    geom_bar(stat = "identity")
  if (facet_wrap) {
    plt <- plt + facet_wrap(~state, scales = "free")
  } else {
    plt <- plt + facet_grid(~state)
  }

  if (!is.null(statepalette)) {
    CEMTscripts:::plot_colorstrip(dat, plt, statepalette)
  } else {
    plt
  }
}
