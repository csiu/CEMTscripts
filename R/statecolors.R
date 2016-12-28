#' State colors
#'
#' @format List of state colors; produced by: \code{
#' statecolors <- c(
#'   "tssactive"="#FF0000",
#'   "tss"="#FF4500",
#'   "txn"="#008000",
#'   "txnweak"="#006400",
#'   "znf"="#66CDAA",
#'   "enhgenic"="#C2E105",
#'   "enhactive"="#FFC34D",
#'   "enh"="#FFFF00",
#'   "bivtss"="#CD5C5C",
#'   "bivenh"="#BDB76B",
#'   "heterochromatin"="#8A91D0",
#'   "repeats"="#C0C0C0",
#'   "polycomb"="#808080",
#'   "low"="#FFFFFF"
#' );
#' save(statecolors, file="data/statecolors.RData")
#' }
#'
#' @source \url{http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html}
"statecolors"

#  ------------------------------------------------------------------------
#' Make statepalette
#' @param x List of states. Valid states are described in \code{?statecolors}
#' @details
#'   For use as \code{statepalette} in \code{tally_state_summary_viz()}
#' @export
statecolors_makepalette <- function(x){
  s <- CEMTscripts::statecolors[x]
  names(s) <- 1L:length(x)
  s
}
