#' Ensembl gene info
#'
#' @format Dataframe with 3 columns: gene_id, gene_name, gene_type
#'
#' @source Ensembl BioMart download of "Homo sapiens genes (GRCh37.p13)"
#'         with attributes: "Ensembl Gene ID", "Associated Gene Name", and "Gene type"
#'         \url{http://grch37.ensembl.org/biomart}
#' @examples
#' ensemblgeneinfo <-
#'   system.file("extdata", "grch37p13.geneinfo.tsv", package="CEMTscripts") %>%
#'   readr::read_tsv()
#' save(ensemblgeneinfo, file="data/ensemblgeneinfo.RData")
"ensemblgeneinfo"
