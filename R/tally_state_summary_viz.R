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
#' @param scale_y_log10
#'          boolean; when this is TRUE, the bin counts will be
#'          log10 transformed and the plot will use points instead of bars.
#' @param shape boolean; when this is TRUE, the bin counts are
#'              rescaled to range 0-1 and the
#'              \code{scale_y_log10} argument is not used
#' @import dplyr
#' @import ggplot2
#' @export
tally_state_summary_viz <- function(input, statepalette=NULL,
                                    facet_wrap=FALSE, scale_y_log10=FALSE,
                                    shape=FALSE){
  dat <-
    input %>%
    dplyr::mutate(state = sub("^[UE]", "", state) %>% as.integer()) %>%
    tidyr::gather(samples, bins, -state) %>%
    dplyr::mutate(samples = as.integer(samples))

  # Plot
  if (shape) {
    dat <-
      dat %>%
      mutate(bins = log10(bins)) %>%
      group_by(state) %>%
      mutate(bins = scales::rescale(bins, to=c(0,1)))
    scale_y_log10 <- TRUE
  }
  if (scale_y_log10) {
    plt <-
      dat %>%
      ggplot(aes(x = samples, y = bins, group = state)) +
      geom_point(size = .75) +
      geom_line() +
      scale_y_log10() +
      ylab("")
  } else {
    plt <-
      dat %>%
      ggplot(aes(x = samples, y = bins)) +
      geom_bar(stat = "identity")
  }
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
