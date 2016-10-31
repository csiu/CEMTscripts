#' Add chr
add_chr <- function(gr){
  if (grepl("^(?!chr)", GenomeInfoDb::seqlevels(gr)[1], perl=T)) {
    GenomeInfoDb::renameSeqlevels(
      gr,
      sub("^", "chr\\1", GenomeInfoDb::seqlevels(gr)))
  } else { gr }
}

# ------------------------------------------------------------------------
#' Label regions with gene ids
#' @param regions GRanges object
#' @param regiontype
#' @param tss.upstream
#'          When \code{regiontype} is "tss", this argument is used
#'          to specify the upstream GenomicRanges::promoter size
#' @param tss.downstream See \code{tss.upstream}
#' @format
#'   Format is the same as \code{regions}, but with an extra "gene_id" column.
#'   The length of the output depends on the number of
#'   matching regions with gene_ids.
#' @export
regions_addgenes <- function(regions, regiontype="tss",
                            tss.upstream=2000, tss.downstream=2000){
  ensure_correct_chr_prefix <-
    grepl("^chr", GenomeInfoDb::seqlevels(regions)[1]) &&
    grepl("^(?!chr)", GenomeInfoDb::seqlevels(ensemblgtf)[1], perl=T)

  if (regiontype == "tss") {
    ## From gene list, infer a list of promoters
    gene_db <-
      GenomicRanges::promoters(
        ensemblgtf, upstream=tss.upstream, downstream=tss.downstream)
    if (ensure_correct_chr_prefix) gene_db<-CEMTscripts:::add_chr(gene_db)

    # ## Overlap genes/promoters with regions of interest
    # GenomicRanges::subsetByOverlaps(query=gene_db, subject=regions)

    ## Find overlaps
    hits <-
      GenomicRanges::findOverlaps(query=gene_db, subject=regions)
    ## Merge subject & query hits
    output <- regions[GenomicRanges::subjectHits(hits),]
    output$gene_id <-
      GenomicRanges::mcols(gene_db)[GenomicRanges::queryHits(hits),"gene_id"]

  } else if (regiontype == "enh") {
    }
  } else {
    stop("regiontype must be one of: 'tss'")
  }
  output
}

#  ------------------------------------------------------------------------
#' Convert regions to genes
#' @param regionsaddgenes
#'          Output of \code{regions_addgenes()}.
#'          Assumes metadata contain columns:
#'          "samples" (number of overlapping samples) and "gene_id"
#' @format
#'   Dataframe with 4 columns: "gene_id", "samples_maxoverlap",
#'   "samples_neighbor", and "label_neighbor".
#'   samples_maxoverlap refers to the maximum number of samples
#'   overlapping the same bin of a gene; samples_neighbor refers
#'   to the samples in the total number of unique samples in the
#'   gene neighborhood, identity of which is given in label_neighborhood.
#' @export
regions_togenes <- function(regionsaddgenes){

  dat_maxoverlap <-
    GenomicRanges::mcols(regionsaddgenes) %>%
    ## Transform to dataframe
    as.data.frame() %>%
    tbl_df() %>%
    ## For each gene, determine max sample overlap
    ## Here we yet to consider the "label"
    select(samples, gene_id, label) %>%
    group_by(gene_id) %>%
    filter(samples == max(samples)) %>%
    unique() %>%
    ungroup()

  ## Create some tmp IDS
  tmp_ids <-
    paste("tmpID", 1:nchar(dat_maxoverlap$label[1]), sep="_")

  ## Here we account for samples labelled to the same gene
  ## but not necessarily overlapping each other
  ## i.e. samples in the same neighborhood
  mat <-
    do.call(rbind, strsplit(dat_maxoverlap$label, ""))
  colnames(mat) <- tmp_ids
  mat <-
    cbind(gene_id=dat_maxoverlap$gene_id, mat)
  mat <-
    mat %>%
    ## Reformatting
    as.data.frame() %>%
    tbl_df() %>%
    mutate_each_(funs(as.character(.)), tmp_ids) %>%

    ## Do the counting of samples
    tidyr::gather(tmpid, call, starts_with("tmpID_")) %>%
    filter(call == 1) %>%
    unique() %>%

    tidyr::spread(tmpid, call)
  mat[is.na(mat)] <- 0
  ## Gather output
  mat <-
    data.frame(
      gene_id = mat$gene_id,
      label_neighbor = do.call(paste0, c(mat[tmp_ids])),
      samples_neighbor = rowSums(data.matrix(mat[tmp_ids]))
    ) %>%
    mutate_all(funs(as.character(.))) %>%
    tbl_df()

  unique(select(dat_maxoverlap, -label)) %>%
    left_join(mat, by="gene_id") %>%
    select(gene_id, samples_maxoverlap=samples, samples_neighbor, label_neighbor)
}
