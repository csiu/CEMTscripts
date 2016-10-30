#' Protein coding transcripts
#'
#' List of protein coding transcripts from Homo sapiens
#' on chromosomes 1-22, X and Y. Data is
#' obtained from Ensembl gene set annotations GR37.75
#' @format
#'   An object of class GRanges of length 81732.
#'   This number refers to the number of
#'   protein_coding transcripts on chromosome 1-22,X,Y.
#'   These transcripts also map to 20154 ensembl gene ids.
#' @source Downloaded from \url{http://www.ensembl.org/info/data/ftp/index.html}
#' @examples
# ensembl_gtf <- "Homo_sapiens.GRCh37.75.gtf"
#' ## Object produce by importing the data & applying filters
#' ensemblgtf <- rtracklayer::import(ensembl_gtf)
#' ensemblgtf <-
#'   ensemblgtf[ensemblgtf$source=="protein_coding" &
#'                ensemblgtf$type=="transcript"]
#' ensemblgtf <-
#'   keepSeqlevels(ensemblgtf, value = c(1:22,"X","Y"))
#'
#' save(ensemblgtf, file="data/ensemblgtf.RData")
"ensemblgtf"
