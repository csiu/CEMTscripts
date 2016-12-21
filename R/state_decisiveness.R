#' Load ChromHMM emission probabilities
#' @param filename The ChromHMM model emissions file
#' @import dplyr
# filename <- file.path(chromhmmdir, "LearnModel_ordered/model/emissions_14.txt")
# (emissions_probs <- load_emissionsdb(filename))
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
#' @details Want to minimize this value
#' @examples
#' state_emissions <- c(1, 0.5)
#' decisiveness(state_emissions)
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
#' @examples
#' emissions_probs <- rbind(c(1, 0.5), c(1,0))
#' genomic_coverage <- c(.95, .05)
#' decisiveness_model(emissions_probs, genomic_coverage)
decisiveness_model <- function(emissions_probs, genomic_coverage=1){
  sum(apply(emissions_probs, 1, CEMTscripts:::decisiveness) * genomic_coverage)
}

#  ------------------------------------------------------------------------
#' Calculate state decisiveness
#'
#' Calculate, based on state emissions, the average decisiveness of a model
#' @param filename.emissions The ChromHMM model emissions file
#' @param genomic_coverage
#'          list containing the genomic coverage of the states.
#'          Used to scale the states.
#' @details Want to minimize this value
#' @export
state_decisiveness <- function(filename.emissions, genomic_coverage){
  ## Load emissions
  CEMTscripts:::load_emissionsdb(filename.emissions) %>%

    ## Tidy emissions
    tidyr::spread(mark, prob) %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("state") %>%

    ## Calculate the average decisiveness across multiple states
    CEMTscripts:::decisiveness_model(genomic_coverage)
}
