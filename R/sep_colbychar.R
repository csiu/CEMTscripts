#' Separate column by each character
#' @param dat dataframe
#' @param column The name of the column you want to split.
#'               This column should be character; if not,
#'               this column will be converted to character.
#' @param column.names The new column names after the split
#' @param keep boolean; keep original column in dataframe?
#' @examples
#' dat <- head(iris)
#' sep_colbychar(dat, column="Petal.Length")
#' @export
sep_colbychar <- function(dat, column="label", column.names=NULL, keep=TRUE){
  ## Create some column IDs
  if (is.null(column.names)) {
    column.names <-
      paste("tmpID", 1:nchar(dat[[column]][1]), sep="_")
  }

  ## Split column & update colnames
  dat[[column]] <- as.character(dat[[column]])
  mat <-
    do.call(rbind, strsplit(dat[[column]], ""))
  colnames(mat) <- column.names

  output <- cbind(dat, mat)
  if (!keep) select_(output, paste0("-",column)) else output
}
