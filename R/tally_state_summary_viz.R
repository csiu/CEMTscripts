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
#' @import dplyr
#' @export
tally_state_summary_viz <- function(input, statepalette=NULL,
                                    facet_wrap=FALSE, scale_y_log10=FALSE){
  require(ggplot2)
  dat <-
    input %>%
    dplyr::mutate(state = sub("^[UE]", "", state) %>% as.integer()) %>%
    tidyr::gather(samples, bins, -state) %>%
    dplyr::mutate(samples = as.integer(samples))

  # Plot
  if (scale_y_log10) {
    plt <-
      dat %>%
      ggplot(aes(x = samples, y = bins, group = state)) +
      geom_point(size = .75) +
      geom_line() +
      scale_y_log10()
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
