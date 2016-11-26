#' Create a list of flanking regions
#' @param gr gr object
#' @param width size of flanking region
#' @export
gr_createFlank <- function(gr, width=2000) {
  upstream_flank <-
    GenomicRanges::flank(gr, width, start=TRUE) %>%
    GenomicRanges::reduce()
  downstream_flank <-
    GenomicRanges::flank(gr, width, start=FALSE) %>%
    GenomicRanges::reduce()
  total_flank <-
    c(upstream_flank, downstream_flank) %>%
    GenomicRanges::reduce()

  list(up = upstream_flank,
       down = downstream_flank,
       both = total_flank)
}

#  ------------------------------------------------------------------------
#' Create a list of boundary regions
#' @param gr gr object
#' @param width size of boundary region
#' @export
gr_createBoundary <- function(gr, width=2000) {
  upstream_border <-
    GenomicRanges::flank(gr, width/2, start=TRUE, both=TRUE) %>%
    GenomicRanges::reduce()
  downstream_border <-
    GenomicRanges::flank(gr, width/2, start=FALSE, both=TRUE) %>%
    GenomicRanges::reduce()
  total_border <-
    c(upstream_border, downstream_border) %>%
    GenomicRanges::reduce()

  list(up = upstream_border,
       down = downstream_border,
       both = total_border)
}
