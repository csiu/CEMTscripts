#' Split histone emission probabilities by state
#' @param emissions dataframe produced by
#'                  \code{chromhmm_loadmodel(...)$emissions}
#' @examples
#' library(dplyr)
#' system.file("extdata", "model_18_core_K27ac.txt", package = "CEMTscripts") %>%
#'  {chromhmm_loadmodel(.)$emissions} %>%
#'  split_by_stateprofile()
split_by_stateprofile <- function(emissions){
  split(setNames(emissions$prob, emissions$markname), emissions$state)
}

#' Spread emissions
#' @param emissions dataframe produced by
#'                  \code{chromhmm_loadmodel(...)$emissions}
chromhmm_spreademissions <- function(emissions){
  mat <-
    emissions %>%
    select(-mark) %>%
    tidyr::spread(markname, prob)
  rownames(mat) <- NULL
  mat <- tibble::column_to_rownames(as.data.frame(mat), "state")
  mat
}
# -------------------------------------------------------------------------
#' Matrix pearson correlation
#' @param x \code{n * p} matrix
#' @param y \code{p * m} matrix
#' @details Computes the pearson correlation of
#'          the rows of x with the columns of y
#' @import dplyr
#' @export
`%cor%` <- function(x,y){
  ## ensure element order is the same
  y <- y[colnames(x),]
  ## compute correlation of each row of x
  ## by each columns of y
  mat <-
    apply(x, 1, function(eachrow){
      cor(eachrow, y, method="pearson")
    }) %>%
    t()
  colnames(mat) <- colnames(y)
  mat
}
# -------------------------------------------------------------------------

library(dplyr)

## Defile model files
# roadmap_core5 <- system.file("extdata", "model_15_coreMarks.txt", package = "CEMTscripts")
roadmap_model_file <- system.file("extdata", "model_18_core_K27ac.txt", package = "CEMTscripts")
model_file <- "test/model_15_unordered.txt"
model_file <- "test/model_15.txt"

# Roadmap labels
rlabel <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx",
            "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF", "9_Het",
            "10_TssBiv", "11_BivFlnk", "12_EnhBiv",
            "13_ReprPC", "14_ReprPCWk", "15_Quies")
rlabel <- c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
            "5_Tx", "6_TxWk",
            "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk",
            "12_ZNF", "13_Het", "14_TssBiv", "15_EnhBiv",
            "16_ReprPC", "17_ReprPCWk", "18_Quies")

## Load models
model <- chromhmm_loadmodel(model_file)
model_roadmap <- chromhmm_loadmodel(roadmap_model_file)

# ## Split each state
# m <- split_by_stateprofile(model$emissions)
# r <- split_by_stateprofile(model_roadmap$emissions)

## Spread emissions
M <- chromhmm_spreademissions(model$emissions)
R <- chromhmm_spreademissions(model_roadmap$emissions)

mat <- R %cor% t(M)
rownames(mat) <- rlabel
pheatmap::pheatmap(mat, cluster_rows=F, cluster_cols=T, row.names=rlabel)
