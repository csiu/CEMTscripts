#' Get regions
#'
#' Get the list of regions overlaping n samples in a particular state
#'
#' @param tally_dir directory containing the state files ending in "bed"
#' @param state integer; state file in \code{tally_dir} to load
#' @param nsamples Filter: How many samples must region overlap?
#'                 Can be single element e.g. "4" or
#'                 a numeric list e.g. "c(3, 4)"
#' @param regex.in Filter: Which samples must be in the result?
#'                 Character length must equal to the number of samples.
#'                 Only "0", "1", and "." characters are permitted, where
#'                 "1" = present, "0" = absent, "." = either present or absent.
#'                 e.g. ".1001" means in the 5 samples:
#'                 sample 1 can be present or absent;
#'                 samples 2 and 5 must be present; and
#'                 samples 3 and 4 must be absent in the results.
#' @param regex.notin
#'                 Filter: Which samples must not be in the result?
#'                 Requirements are the same as \code{regex.in}
#' @export
regions_get <- function(tally_dir, state, nsamples="all",
                        regex.in="*", regex.notin=NULL){
  validate_regex <- function(regex, l){
    if (!(grepl("^[01\\.]+$", regex) && nchar(regex)==l)) {
      stop(paste("regex is not valid:",
                 "must be [01?] and length be equal to number of samples"))
    } else {
      paste0("^", regex, "$")
    }
  }

  ## Load data
  input_file <- CEMTscripts:::getfile_tallystate(tally_dir, state)
  state_tally <- CEMTscripts:::load_tallystate(input_file)
  ## Load as GR
  gr <- GenomicRanges::makeGRangesFromDataFrame(state_tally)
  GenomicRanges::mcols(gr) <-
    dplyr::select(state_tally, samples, label, bins)

  ## Number of samples
  total_samples <- nchar(state_tally[[1,"label"]])

  ## Filter: How many samples in result?
  if (!(nsamples=="all" || nsamples==total_samples)) {
    gr <-
      gr[GenomicRanges::mcols(gr)$samples %in% nsamples]
  }
  ## Filter: Which samples must be in the result?
  if (regex.in!="*") {
    regex.in <-
      validate_regex(regex.in, total_samples)
    gr <-
      gr[grepl(regex.in, GenomicRanges::mcols(gr)$label)]
  }
  ## Filter: Which samples must not be in the result?
  if (!is.null(regex.notin)) {
    regex.notin <-
      validate_regex(regex.notin, total_samples)
    gr <-
      gr[!grepl(regex.notin, GenomicRanges::mcols(gr)$label)]
  }

  return(gr)
}
