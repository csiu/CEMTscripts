#' make Expression counts for UpSetR plots
#'
#' Count groups for for use in
#' UpSetR::fromExpression(.) %>% UpSetR::upset()
#'
#' @param input either a dataframe or gr object
#' @param isgr boolean; is \code{input} a gr object?
#'             (i.e. output of \code{regions_get()} or
#'             \code{regions_addgenes()})
#' @param column.group Name of column in mcols(gr) to group by;
#'                     should contain only "0" and "1",
#'                     length be equal to number of samples,
#'                     "0" = sample is absent and
#'                     "1" = sample is present.
#' @param column.tally Name of column in mcols(gr) to add up
#'                     when we group by \code{column.group}.
#'                     Default NULL, we count each row in gr as 1.
#' @param group.name character vector of column.group sample names
#' @export
makeExpression_fromLabelcol <- function(input, isgr=TRUE,
                                        column.group="label",
                                        column.tally=NULL,
                                        group.names=NULL){
  ## Set defaults
  ## - Make input to dataframe
  ## - Ensure there is tally column
  ## - Ensure there is group name
  if (isgr) {
    x <-
      GenomicRanges::mcols(input) %>%
      as.data.frame() %>%
      tbl_df()
  } else {
    x <- input
  }
  if (is.null(column.tally)) {
    x <- mutate(x, columntally=1)
    column.tally <- "columntally"
  }
  if (is.null(group.names)) {
    group.names <-
      paste("group", 1:nchar(x[[column.group]][1]), sep="_")
  }

  ## Create list and
  ## Update the names
  x <-
    x %>%
    group_by_(column.group) %>%
    summarise_(count = sprintf("sum(%s)", column.tally)) %>%
    {setNames(.$count, .[[column.group]])}
  names(x) <-
    strsplit(names(x), "") %>%
    sapply(function(s){paste(group.names[s==1], collapse="&")})
  x
}
