#' Total unique elements in embedded list
#' @param x a list whereby each element is a character vector/set and
#'          the name of the element is the name of the set
#' @param getunique boolean; if TRUE, will instead return the
#'                  list of unique elements
set_total <- function(x, getunique=FALSE) {
  elements <- Reduce(union, x)
  if (getunique) elements else length(elements)
}

#  ------------------------------------------------------------------------
#' Plot sets
#'
#' @param x a list whereby each element is a character vector/set and
#'          the name of the element is the name of the set
#' @param type what kind of plot to produce: one of "upset" (default)
#'             or "venn"
#' @param ... additional arguments to \code{UpSetR::upset(...)} or
#'            \code{VennDiagram::venn.diagram(...)}
#' @import dplyr
#' @export
#' @examples
#' x <-
#'   list(
#'     set1 = letters[1:10],
#'     set2 = c("a", "e", "i", "o", "u"),
#'     set3 = letters[1:3]
#'   )
#' plot_sets(x)
#' plot_sets(x, type="venn")
plot_sets <- function(x, type="upset", ...){
  if (type=="upset") {
    UpSetR::fromList(x) %>%
      UpSetR::upset(...)
  } else if (type=="venn") {
    ## ignore VennDiagram*log
    futile.logger::flog.threshold(futile.logger::ERROR,
                                  name = "VennDiagramLogger")

    grid::grid.newpage()
    VennDiagram::venn.diagram(x, filename=NULL, ...) %>%
      grid::grid.draw()
  }
}
