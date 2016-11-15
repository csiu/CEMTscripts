#' Plot overlap enrichment
#' @param input Output of \code{chromhmm_loadenrichment}
#' @param ... arguments to pheatmap::pheatmap(...)
#' @import ggplot2
#' @export
chromhmm_pltenrichment <- function(input, ...){
  input %>%
    as.data.frame() %>%
    tibble::column_to_rownames("state") %>%
    data.matrix() %>%
    pheatmap::pheatmap(
      scale = "column",
      ## colors from: viridis::viridis(5, option="B")
      color = colorRampPalette(
        c("#000004FF", "#56106EFF", "#BB3754FF",
          "#F98C0AFF", "#FCFFA4FF"))(100),
      ...)
}
