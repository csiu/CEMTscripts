#' Load ChromHMM emission probabilities
#' @param filename The ChromHMM model emissions file
#' @import dplyr
# filename <- file.path(chromhmmdir, "LearnModel_ordered/model/emissions_19.txt")
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
#' Calculate the unsynchronization of 1 state
#' @param state_emissions list of state emissions from histone modifications
#' @details Calculated by \code{d_k = sum[h=1 to H]{min(1-E_h, E_h)}}
#'          want to minimize this value.
# state_emissions <- c(0.25, 0.5)
# CEMTscripts:::unsync(state_emissions)
unsync <- function(state_emissions){
  sum(sapply(state_emissions, function(e){min(1-e, e)}))
}

#  ------------------------------------------------------------------------
#' Calculate the unsynchronization of the model
#' @param emissions_probs matrix containing the state emssions
#' @param scale
#'          list used to scale the states;
#'          by default, there is no scaling.
#' @param returntotal
#'          boolean; return the total value (TRUE) or
#'          the list of values making up the total (FALSE)
#' @details Calculated by \code{D = sum[k=1 to K]{d_k * s_k}};
#'          want to minimize this value. s_k represent the scaled term.
#' @examples
#' emissions_probs <- rbind(c(1, 0.5), c(1,0))
#' scale <- c(.95, .05)
#' CEMTscripts:::unsynchronization(emissions_probs, scale)
unsynchronization <- function(emissions_probs, scale=1, returntotal=TRUE){
  d  <- apply(emissions_probs, 1, CEMTscripts:::unsync) * scale
  if (returntotal) sum(d) else d
}

#  ------------------------------------------------------------------------
#' Calculate model unsynchronization from chromhmm emissions file
#'
#' Calculate, based on state emissions, the unsynchronization of a model
#' @param filename.emissions The ChromHMM model emissions file
#' @param usescale
#'          boolean; use scale (1/k) terms in calculation?
#' @param returntotal
#'          boolean; return the total value (TRUE) or
#'          the list of values making up the total (FALSE)
#' @details Calculated by \code{D = sum[k=1 to K]{d_k * s_k}} where
#'          \code{d_k = sum[h=1 to H]{min(1-E_h, E_h)}} and s_k is the
#'          scale term. Note: want to choose a model which minimizes the
#'          unsynchronization (D) term.
#' @export
unsynchronization_score <- function(filename.emissions, usescale=TRUE,
                                    returntotal=TRUE){
  emissions_probs <-
    ## Load emissions
    CEMTscripts:::load_emissionsdb(filename.emissions) %>%

    ## Tidy emissions
    tidyr::spread(mark, prob) %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("state")

  ## Calculate the unsynchronization of the model
  if (usescale) scale<-1/nrow(emissions_probs) else scale<-1
  CEMTscripts:::unsynchronization(emissions_probs, scale, returntotal)
}
