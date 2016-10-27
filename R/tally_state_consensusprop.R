#' Consensus proportion
#'
#' Proportion of bins marked-in-all-samples among bins marked-in-at-least-1-sample
#'
#' @param tallystatesummary Output of \code{tally_state_summary()}
#' @details
#'   The number of bins marked as state \code{i} in the consensus divided by
#'   the number of bins marked as state \code{i} in at least 1 sample
#'   (...which gives somewhat a measure of variability).
#' @format Dataframe with 4 columns sorted by increasing "proportion" size
#'   \itemize{
#'    \item{\code{state}}
#'    \item{\code{bins_consensus}}
#'    \item{\code{bins_atleast1}}
#'    \item{\code{proportion}}
#'   }
#' @export
tally_state_consensusprob <- function(tallystatesummary){
  ## Reformat
  d <-
    tallystatesummary %>%
    as.data.frame() %>%
    gather(nsamples, bins, -state)

  ## Number of bins marked in at least 1 sample
  d_atleast1 <-
    d %>%
    group_by(state) %>%
    summarize(bins_atleast1 = sum(bins))

  ## Consensus size
  N <- max(as.integer(d$nsamples))

  ## Output
  d %>%
    filter_(paste("nsamples ==", N)) %>%
    dplyr::left_join(d_atleast1, by="state") %>%
    mutate(
      proportion = bins/bins_atleast1
    ) %>%
    arrange(proportion) %>%
    select(-nsamples) %>%
    dplyr::rename(bins_consensus = bins)
}
