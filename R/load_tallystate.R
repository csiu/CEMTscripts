#' Find tally state file in tally dir
#' @param tally_dir directory containing the state files ending in "bed"
#' @param state Specify state file to be returned.
#'              If \code{NULL} (default) is used,
#'              all files bed files in tally_dir will be returned.
getfile_tallystate <- function(tally_dir, state=NULL){
  # List of all files in tally dir
  # Grep for state of interest
  input_files <- list.files(path=tally_dir, pattern="*bed")

  if (!is.null(state)){
    file.path(tally_dir,
              grep(sprintf("[\\.UE]%s\\.", state), input_files, value=T))
  } else {
    input_files
  }
}

#  ------------------------------------------------------------------------
#' Load unionbedg tally state file
#'
#' @param input_file TSV file with 6 columns + no header:
#'                   chrom, start, stop,
#'                   number of samples, sample representation
#'                   e.g. 0010 (where 0=absent; 1=present),
#'                   number of bins in region i.e. (stop-start)/200bp.
#' @import readr
load_tallystate <- function(input_file){
  read_tsv(input_file,
           col_names = c("chrom", "start", "end",
                         "samples", "label", "bins"),
           col_types = cols(
             chrom = col_character(),
             start = col_integer(),
             end = col_integer(),
             samples = col_integer(),
             label = col_character(),
             bins = col_integer()
           ))
}
