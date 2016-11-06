#' Cumulative sum of a (desc) sorted vector
#'
#' @param x vector of elements
#' @examples
#' x <- c("a"=1, "b"=3, "c"=NA, "d"=6, "e"=2, "f"=NA)
#' CEMTscripts:::transcript_cumsum(x)
transcript_cumsum <- function(x){
  ## remove NA
  x <- x[!is.na(x)]
  ## sort
  x <- sort(x, decreasing=TRUE)
  ## cumsum
  cumsum(x)
}

#  ------------------------------------------------------------------------
#' Transcript abundance
#'
#' @param txi_df one element (i.e. "abundance") from \code{transcript_togene}
#' @param expr_column column containing transcript information
#' @format Dataframe whose rows are sorted in descending order based on \code{expr_column}.
#'         There are 5 columns: gene_id, expr, cumsum, proportion, and count
#' @export
transcript_abundance <- function(txi_df, expr_column){
  x <- setNames(txi_df[[expr_column]], txi_df[["gene_id"]])
  x <- x[!is.na(x)]

  cs <- CEMTscripts:::transcript_cumsum(x)
  total <- sum(x)

  data.frame(
    gene_id = names(cs),
    expr = x[names(cs)],
    cumsum = cs,
    proportion = cs/total,
    count = 1:length(cs)
  ) %>%
    as.data.frame() %>%
    tbl_df()
}

#  ------------------------------------------------------------------------
#' Plot Transcript abundance
#'
#' Use a dot plot to plot the transcript abundance
#' @param transcriptabundance output of \code{transcript_abundance}
#' @param ngenes integer; how many top expr genes to plot?
#' @import ggplot2
#' @export
transcript_abundance_viz <- function(transcriptabundance, ngenes=10000){
  transcriptabundance %>%
    head(ngenes) %>%
    ggplot(aes(x=count, y=proportion)) +
    geom_point() +
    ylab("Proportion of transcripts") +
    xlab("Number of genes") +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "white")
      )
}

#  ------------------------------------------------------------------------
#' Return the list of top genes
#' @param transcriptabundance output of \code{transcript_abundance}
#' @param ngenes integer; how many top expr genes?
ta_getgenes <- function(transcriptabundance, ngenes){
  transcriptabundance %>%
    filter(count <= ngenes) %>%
    {.$gene_id} %>%
    as.character()
}

#  ------------------------------------------------------------------------
#' Proportion of expressed transcripts that make up the top ngenes
#' @param transcriptabundance output of \code{transcript_abundance}
#' @param ngenes integer; how many top expr genes?
#' @param getgenes
#'          boolean; if TRUE, will instead return the list of top ngenes
#' @export
ta_ngenes2proportion <- function(transcriptabundance, ngenes,
                                 getgenes=FALSE){
  if (getgenes){
    ## the genes
    CEMTscripts:::ta_getgenes(transcriptabundance, ngenes)
  } else {
    ## the proportion
    transcriptabundance %>%
      filter(count == ngenes) %>%
      {.$proportion}
  }
}

#  ------------------------------------------------------------------------
#' Number of genes which make up a proportion of expressed transcript
#' @param transcriptabundance output of \code{transcript_abundance}
#' @param p proportion. Value should be between 0 - 1.
#' @param gt boolean. If FALSE (default), the function will return the
#'           number of genes equal to or below the specified proportion.
#'           If TRUE, the function will return +1 gene which makes
#'           the proportion slightly larger than the specified value
#' @param getgenes
#'          boolean; if TRUE, will instead return the list of top ngenes
#' @export
ta_proportion2ngenes <- function(transcriptabundance, p, gt=FALSE,
                                 getgenes=FALSE){
  the_ngenes <-
    transcriptabundance %>%
    filter(proportion <= p) %>%
    {.$count} %>%
    tail(1)
  if (gt) the_ngenes<-the_ngenes+1

  if (getgenes) {
    CEMTscripts:::ta_getgenes(transcriptabundance, the_ngenes)
  } else {
    the_ngenes
  }
}
