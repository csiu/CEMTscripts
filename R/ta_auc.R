#' Compute middle Riemann sum
#' @param transcriptabundance output of \code{transcript_abundance}
#' @param ngenes integer; how many top expr genes?
#' @export
ta_auc <- function(transcriptabundance, ngenes=10000){
  cumsum_list <- transcriptabundance$cumsum[1:ngenes]

  ## Set up widths and heights so that width * height = 1
  x <- (1:length(cumsum_list)) / length(cumsum_list)
  y <- cumsum_list / max(cumsum_list)

  id <- order(x)
  ## middle Riemann sum
  ## sum of (mini widths * rollover mean heights)
  sum(diff(x[id]) * zoo::rollmean(y[id],2))
}
