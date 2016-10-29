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
#' Do the particular states of a genomic feature
#' between the segmentation of a sample match?
#'
#' @param filenames
#'           Vector containing the paths of 2 ChromHMM segmentations.
#'           Element names will be used as the id.
#' @param state
#'          List specifying the states in a genomic feature
#'          to be treated the same.
#' @details
#'   Intended to be used to compare the states between the same bins
#'   of the same sample generated using different models
#' @format
#'   List object with 3 names
#'   \itemize{
#'     \item{\code{raw}}{
#'        Dataframe produced with bedtools unionbedg}
#'     \item{\code{summary_viz}}{
#'        ggplot summarizing the number of bins in the states of interest
#'        in the two samples}
#'     \item{\code{summary}}{
#'        Dataframe containing the number of bins which have matching states}
#'   }
#' @import ggplot2
#' @import dplyr
#' @export
regions_binmatch <- function(filenames, state){
  tmpfile_unionbedg <- tempfile()

  ## Use system's bedtools to run command
  CEMTscripts:::bedtools_unionbedg(filenames, tmpfile_unionbedg)
  ## Load data
  unionbedg <- CEMTscripts:::load_unionbedg(tmpfile_unionbedg)

  total_bins <- sum(unionbedg$nbins)

  ## Detailed summary
  unionbedg_filtered <-
    unionbedg %>%
    dplyr::rename_(x=names(filenames)[1], y=names(filenames)[2]) %>%
    filter(x%in%state[[1]] | y%in%state[[2]]) %>%
    group_by(x, y) %>%
    summarise(nbins = sum(nbins)) %>%
    ungroup()
  ## Crude summary plot
  gg_sumary <-
    unionbedg_filtered %>%
    ggplot(aes(
      x = reorder(x, -nbins),
      y = reorder(y, -nbins),
      fill=log(nbins))) +
    geom_tile() +
    scale_fill_gradient(
      #midpoint = 9, mid="yellow",
      low="white", high="red") +
    geom_text(aes(label=nbins), angle=45, size=3) +
    xlab(names(filenames)[1]) +
    ylab(names(filenames)[2]) +
    theme_bw() +
    theme(legend.position = "none")

  ## Generic summary
  unionbedg_filtered2 <-
    unionbedg_filtered %>%
    mutate(
      x = x %in% state[[1]],
      y = y %in% state[[2]]
      ) %>%
    group_by(x, y) %>%
    summarise(nbins = sum(nbins)) %>%
    ungroup() %>%
    mutate(proportion = nbins/total_bins) %>%
    arrange(desc(nbins))
  colnames(unionbedg_filtered2)[1:2] <- names(filenames)
  unionbedg_filtered2

  ## Remove temp file
  file.remove(tmpfile_unionbedg)

  list("raw" = unionbedg,
       "summary_viz" = gg_sumary,
       "summary" = unionbedg_filtered2)
}
