#' Run Bedtools unionbedg
#'
#' @param bgfiles Character list containing the list of BedGraphs to combine.
#'                The names of the list elements will be used as output
#'                header names. If no element names are set,
#'                the file basename will be used instead.
#' @details Invoke System Command: \code{bedtools unionbedg}
bedtools_unionbedg <- function(bgfiles, outfile){
  if (is.null(names(bgfiles))) names(bgfiles)<-basename(bgfiles)

  e <-system2("bedtools",
              args = c("unionbedg", "-i", bgfiles,
                       "-header", "-names", names(bgfiles)),
              stdout = outfile)
  if (e != 0) stop("Failed to execute 'bedtools unionbedg'")
}

#  ------------------------------------------------------------------------
#' State comparison
#'
#' Do the states between the segmentation match?
#'
#' @details
#'   Intended to be used to compare the states between the same bins
#'   of the same sample generated using different models
regions_binmatch <- function(x.filename, y.filename, x.state, y.state){

}
#  ------------------------------------------------------------------------
