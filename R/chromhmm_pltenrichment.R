#' Combine matrices
#' @param matrix_list a list where each element is a matrix
#' @param method How should we combine the matrices?
#'               One of "mean" (default),
#'               standard deviation "sd", or
#'               coefficient of variation "cov"
combine_matrix <- function(matrix_list, method="mean"){
  if (method=="mean") {
    apply(simplify2array(matrix_list), 1:2, mean)
  } else if (method=="sd") {
    apply(simplify2array(matrix_list), 1:2, sd)
  } else if (method=="cov") {
    apply(simplify2array(df_list), 1:2, function(x){sd(x)/mean(x)})
  }
}

#  ------------------------------------------------------------------------
#' Combine a list of matrices and plot heatmap
#' @param mat either a list where each element is a matrix,
#'            or a matrix object
#' @param is_list boolean; is \code{mat} a list?
#' @param method when \code{mat} is a list of matrices,
#'               how should we combine the matrices?
#'               "mean" (default), "sd", or "cov" (coefficient of variation)
#' @param use_manual_opt boolean; choose your own pheatmap parameters?
#' @param ... additional arguments to pheatmap::pheatmap
#' @export
plot_combinematrix <- function(mat,
                               is_list=TRUE, method="mean",
                               use_manual_opt=FALSE, ...){
  if (is_list) {
    mat <- CEMTscripts:::combine_matrix(mat, method)
  }

  if (!use_manual_opt) {
    if (method=="mean"){
      scale <- "column"
      color <-
        colorRampPalette(c("white", "#e6e6ff", "blue"))(100)
      use_manual_opt <- FALSE
    } else if (method=="cov" || method == "sd") {
      scale <- "none"
      color <-
        colorRampPalette(c("#DDF3DF", "white", "#E34234"))(100)
      use_manual_opt <- FALSE
    }
  }

  if (use_manual_opt) {
    pheatmap::pheatmap(mat, ...)
  } else {
    pheatmap::pheatmap(
      mat,
      scale = scale,
      color = color,
      ...)
  }
}
