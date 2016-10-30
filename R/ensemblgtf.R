#' Protein coding transcripts
#'
#' List of protein coding transcripts from Homo sapiens
#' obtained from Ensembl gene set annotations GR37.75
#' @source Downloaded from \url{http://www.ensembl.org/info/data/ftp/index.html}
#' @examples
#' ## Object produce by the following:
#'
#' # ensembl_gtf <- "Homo_sapiens.GRCh37.75.gtf"
#' ensemblgtf <- rtracklayer::import(ensembl_gtf)
#' ensemblgtf <-
#'   ensemblgtf[ensemblgtf$source=="protein_coding" &
#'                ensemblgtf$type=="transcript"]
#' save(ensemblgtf, file="data/ensemblgtf.RData")
"ensemblgtf"

