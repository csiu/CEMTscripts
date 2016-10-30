#' Label regions with gene ids
#' @param regions GRanges object
#' @param regiontype
#' @format
#'   Format is the same as \code{regions}, but with an extra "gene_id" column.
#'   The length of the output depends on the number of
#'   matching regions with gene_ids.
#' @import GenomicRanges
regions_addgenes <- function(regions, regiontype="tss",
                            tss.upstream=2000, tss.downstream=2000){
  if (regiontype == "tss") {
    ## From gene list, infer a list of promoters
    gene_db <-
      GenomicRanges::promoters(
        ensemblgtf, upstream=tss.upstream, downstream=tss.downstream)
    if (grepl("^chr", seqlevels(regions)[1]) &&
        grepl("^(?!chr)", seqlevels(gene_db)[1], perl=T)) {
      gene_db <-
        renameSeqlevels(
          gene_db,
          sub("^", "chr\\1", seqlevels(gene_db)))
    }
  } else {
    stop("regiontype must be one of: 'tss'")
  }
  # ## Overlap genes/promoters with regions of interest
  # GenomicRanges::subsetByOverlaps(query=gene_db, subject=regions)

  ## Find overlaps
  hits <-
    GenomicRanges::findOverlaps(query=gene_db, subject=regions)
  ## Merge subject & query hits
  dat.hits <- regions[subjectHits(hits),]
  dat.hits$gene_id <- mcols(gene_db)[queryHits(hits),"gene_id"]
  dat.hits
}

