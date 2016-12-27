#' Load ChromHMM emission probabilities
#' @param filename The ChromHMM model emissions file
#' @import dplyr
# filename <- file.path(chromhmmdir, "LearnModel_ordered/model/emissions_14.txt")
# (emissions_probs <- CEMTscripts:::load_emissionsdb(filename))
load_emissionsdb <- function(filename){
  emissions_probs <- readr::read_tsv(filename, progress=FALSE)
  colnames(emissions_probs)[1] <- "state"
  emissions_probs <-
    emissions_probs %>%
    tidyr::gather(mark, prob, -state)
  emissions_probs
}

#  ------------------------------------------------------------------------
#' Calculate the decisiveness of 1 state
#' @param state_emissions list of state emissions from histone modifications
#' @details Calculated by \code{d_k = sum[h=1 to H]{min(1-E_h, E_h)}}
#'          want to minimize this value.
# state_emissions <- c(0.25, 0.5)
# CEMTscripts:::decisiveness(state_emissions)
decisiveness <- function(state_emissions){
  sum(sapply(state_emissions, function(e){min(1-e, e)}))
}

#  ------------------------------------------------------------------------
#' Calculate the decisiveness of model
#' @param emissions_probs matrix containing the state emssions
#' @param genomic_coverage
#'          list containing the genomic coverage of the states.
#'          Used to scale the states.
#'          By default, there is no scaling
#' @param returntype specify what value is returned:
#'                   "d_k" for vector of d_k;
#'                   "dc_k" for scaled vector of d_k * genomic coverage;
#'                   anything else (default) for the metric itself
#' @details Calculated by \code{D = sum[k=1 to K]{d_k * c_k}};
#'          want to minimize this value.
#' @examples
#' emissions_probs <- rbind(c(1, 0.5), c(1,0))
#' genomic_coverage <- c(.95, .05)
#' CEMTscripts:::decisiveness_model(emissions_probs, genomic_coverage)
decisiveness_model <- function(emissions_probs, genomic_coverage=1,
                               returntype=""){
  d_k  <- apply(emissions_probs, 1, CEMTscripts:::decisiveness)
  dc_k <-  d_k * genomic_coverage

  if (returntype=="d_k") {
    d_k
  } else if (returntype=="dc_k") {
    dc_k
  } else {
    sum(dc_k)
  }
}

#  ------------------------------------------------------------------------
#' Calculate state decisiveness
#'
#' Calculate, based on state emissions, the average decisiveness of a model
#' @param filename.emissions The ChromHMM model emissions file
#' @param genomic_coverage
#'          list containing the genomic coverage of the states.
#'          Used to scale the states.
#' @param returntype specify what value is returned:
#'                   "d_k" for vector of d_k;
#'                   "dc_k" for scaled vector of d_k * genomic coverage;
#'                   anything else (default) for the metric itself
#' @details Calculated by \code{D = sum[k=1 to K]{d_k * c_k}} where
#'          \code{d_k = sum[h=1 to H]{min(1-E_h, E_h)}};
#'          want to minimize this value.
#' @export
state_decisiveness <- function(filename.emissions, genomic_coverage,
                               returntype=""){
  ## Load emissions
  CEMTscripts:::load_emissionsdb(filename.emissions) %>%

    ## Tidy emissions
    tidyr::spread(mark, prob) %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("state") %>%

    ## Calculate the average decisiveness across multiple states
    CEMTscripts:::decisiveness_model(genomic_coverage, returntype)
}
