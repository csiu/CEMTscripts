#' Add chr prefix to GRanges object
#'
#' Add "chr" prefix to a GRanges object that is lacking "chr"
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
#'          One of "tss"  or "enh".
#'          "tss" will make associations between regions
#'          overlapping gene promoter neighborhoods;
#'          "enh" will make associations between regions
#'          upstream (i.e "precede") of genes
#'          for both positive and negative strands
#'          if regions strand information is unspecified.
#' @param tss.upstream
#'          When \code{regiontype} is "tss", this argument is used
#'          to specify the upstream GenomicRanges::promoter size.
#' @param tss.downstream See \code{tss.upstream}
#' @param enh.distmax
#'          When \code{regiontype} is "enh", this argument is used
#'          to specify the maximum region to gene distance.
#'          Default (NULL) is to consider all region to gene pairings.
#' @format
#'   GRanges object containing all columns from \code{regions}.
#'   A "gene_id" column is also added.
#'   The length of the output depends on the number of
#'   regions matching with gene_ids.
#'   For \code{regiontype="enh"},
#'   "gene_distance" (distance between the region and the gene) and
#'   "ensemblgtf_id" (\code{ensemblgtf} dataframe row number of gene)
#'   columns are also added.
#' @export
regions_addgenes <- function(regions, regiontype="tss",
                            tss.upstream=2000, tss.downstream=2000,
                            enh.distmax=NULL){
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
    gene_db <- ensemblgtf
    if (ensure_correct_chr_prefix) gene_db<-CEMTscripts:::add_chr(gene_db)

    ## Ensure there is strand info on regions
    if (all(GenomicRanges::strand(regions)=="*")) {
      GenomicRanges::strand(regions) <- "+"
      r <- regions
      GenomicRanges::strand(regions) <- "-"
      r <- append(r, regions)

      regions <- r
      rm(r)
    }

    ## Find upstream enhancers
    hits <-
      GenomicRanges::precede(regions, gene_db, ignore.strand=FALSE)
    ## Add gene id + subject (ensemblgtf) row id
    regions$gene_id <- gene_db$gene_id[hits]
    regions$ensemblgtf_id <- hits

    ## Remove NAs
    regions <- regions[!is.na(regions$gene_id)]

    ## Add distances
    regions$gene_distance <-
      GenomicRanges::distance(regions, gene_db[regions$ensemblgtf_id])

    output <- regions
    if (!is.null(enh.distmax)) {
      output <- output[output$gene_distance<=enh.distmax]
    }

  } else {
    stop("regiontype must be one of: 'tss', 'enh'")
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
