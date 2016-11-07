#' Load MACS2 NAME_peaks.narrowPeak
#' @source \url{https://github.com/taoliu/MACS}
#' @export
load_macs2_narrowpeak <- function(narrowPeak_file) {
  dat <-
    readr::read_tsv(narrowPeak_file,
                    col_names = c(
                      "chr", "start", "end",
                      "name", "score", "strand",
                      "fc", "nlogp", "nlogq", "relsummit"),
                    col_types = readr::cols(
                      chr = readr::col_character()))
  gr <-
    GenomicRanges::makeGRangesFromDataFrame(dat)
  GenomicRanges::mcols(gr) <-
    select(dat, -chr, -start, -end, -strand)
  gr
}
