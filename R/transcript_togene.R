#' Create tx2gene for use in tximport
#' @format a two-column data.frame linking
#'         transcript id (column 1) to gene id (column 2)
#' @source Ensembl BioMart download of "Homo sapiens genes (GRCh37.p13)"
#'         with attributes: "Ensembl Gene ID" and "Ensembl Transcript ID"
#'         \url{http://grch37.ensembl.org/biomart}
#tx2gene <- transcript_tx2gene()
transcript_tx2gene <- function(){
  system.file("extdata", "grch37p13.tx2gene.tsv",
              package="CEMTscripts") %>%
    readr::read_tsv(col_types = readr::cols(.default = "c"))
}

#  ------------------------------------------------------------------------
#' Map transcripts to genes
#'
#' Uses tximport to map Salmon transcripts to genes
#' @param files a character vector of filenames
#'              for the transcript-level abundances
transcript_tximport <- function(files, ...){
  tximport::tximport(files,
                     type = "salmon",
                     tx2gene = CEMTscripts:::transcript_tx2gene(),
                     reader = readr::read_tsv,
                     ...)
}

#  ------------------------------------------------------------------------
#' Summarize genes
#'
#' Uses tximport to map Salmon transcripts to genes.
#' Column names/sample ids are added at this step.
#' @param files a character vector of filenames
#'              for the transcript-level abundances
#' @param file_ids
#'          a character vector of ids in the same order as \code{files}.
#'          By default, the directory name of the base file is used.
#' @export
transcript_togene <- function(files, file_ids=NULL){
  txi <- CEMTscripts:::transcript_tximport(files)

  ## Add column name
  if (is.null(file_ids)) file_ids<-basename(dirname(files))
  for (i in c("abundance", "counts", "length"))
    colnames(txi[[i]]) <- file_ids

  txi
}
