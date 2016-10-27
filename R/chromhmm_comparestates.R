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
#' Compute pearson correlation between the emission probabilities of different models
#' @param model_file ChromHMM model file
#' @param reference_file ChromHMM model file for comparison
#' @import dplyr
#' @export
chromhmm_comparestates <- function(model_file, reference_file){
  ## Load 2 emissions matrices
  M <-
    chromhmm_loadmodel(model_file) %>%
    {chromhmm_spreademissions(.$emissions)}
  R <-
    chromhmm_loadmodel(reference_file) %>%
    {chromhmm_spreademissions(.$emissions)}

  ## Compute pearson correlation
  R %cor% t(M)
}

# -------------------------------------------------------------------------
#' Compute pearson correlation between the emission probabilities of different models
#' @param chromhmmcomparestates
#'          Output of \code{chromhmm_comparestates()}
#' @param ref
#'          Used to label the states of the reference model if the reference
#'          model is Roadmap's 15 state core model "\code(roadmap5marks)"
#'          or 18 state expanded model "\code(roadmap6marks)"
#' @param ...
#'          Additional argument to \code{pheatmap(...)}
#' @export
chromhmm_comparestates_vizH <- function(chromhmmcomparestates, ref=NULL, ...){
  if (ref=="roadmap5marks") {
    rlabel <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx",
                "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF", "9_Het",
                "10_TssBiv", "11_BivFlnk", "12_EnhBiv",
                "13_ReprPC", "14_ReprPCWk", "15_Quies")
  } else if (ref=="roadmap6marks") {
    rlabel <- c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
                "5_Tx", "6_TxWk",
                "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk",
                "12_ZNF", "13_Het", "14_TssBiv", "15_EnhBiv",
                "16_ReprPC", "17_ReprPCWk", "18_Quies")
  }
  if (ref %in% c("roadmap5marks", "roadmap6marks")) rownames(chromhmmcomparestates)<-rlabel
  pheatmap::pheatmap(chromhmmcomparestates, cluster_rows=F, cluster_cols=T)
}

# -------------------------------------------------------------------------
library(dplyr)

## Defile model files
# roadmap_core5 <- system.file("extdata", "model_15_coreMarks.txt", package = "CEMTscripts")
roadmap_model_file <- system.file("extdata", "model_18_core_K27ac.txt", package = "CEMTscripts")
#model_file <- "test/model_15_unordered.txt"
model_file <- "test/model_15.txt"
reference_file <- roadmap_model_file

mat <- chromhmm_comparestates(model_file, reference_file)
chromhmm_comparestates_vizH(mat, ref="roadmap6marks")
