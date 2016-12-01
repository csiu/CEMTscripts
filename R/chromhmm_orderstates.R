#' Load ChromHMM model
#'
#' @param model_file
#'          ChromHMM model file.
#'          This file usually starts with "model_" and ends with ".txt"
#' @param histone_order
#'          list of histone order. Example
#'          \code{c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3",
#'                  "H3K9me3", "H3K27me3")}
#' @return Returns a list object.
#' \itemize{
#'  \item{"\code{num_states}"}{
#'    for the number of states in the model}
#'  \item{"\code{emissions}"}{
#'    for a dataframe with the emission probabilities}
#'  \item{"\code{transitions}"}{
#'    for a dataframe with the transition probabilities}
#' }
#' @examples
#' model_file <-
#'   system.file("extdata", "model_15_coreMarks.txt",
#'               package="CEMTscripts")
#' model_file <-
#'   system.file("extdata", "model_18_core_K27ac.txt",
#'               package="CEMTscripts")
#'
#' chromhmm_loadmodel(model_file)
#'
#' @import readr
#' @export
chromhmm_loadmodel <- function(model_file, histone_order=NULL){
  ## Load data
  model <-
    read_tsv(model_file, col_names = paste0("X", 1:6), col_types = cols(.default = "c"))

  num_states <- as.numeric(model[1,1])

  ## Extract transitions
  transitions <-
    model %>%
    filter(X1 == "transitionprobs") %>%
    select(
      from = X2,
      to = X3,
      prob = X4
    ) %>%
    mutate(
      prob = as.numeric(prob)
    )
  transitions$from <-
    factor(transitions$from, levels = 1:num_states)
  transitions$to <-
    factor(transitions$to, levels = 1:num_states)

  ## Extract emissions
  emissions <-
    model %>%
    filter(
      X1 == "emissionprobs",
      X5 == "1"
    ) %>%
    select(
      state = X2,
      mark = X3,
      markname = X4,
      prob = X6
    ) %>%
    mutate(
      prob = as.numeric(prob)
    )
  emissions$state <-
    factor(emissions$state, levels = 1:num_states)
  if (!is.null(histone_order)) {
    emissions$markname <-
      factor(emissions$markname, levels = histone_order)
  }

  list(
    "num_states" = num_states,
    "emissions" = emissions,
    "transitions" = transitions)
}

# -------------------------------------------------------------------------
#' Order states
#'
#' For a given mark and emission probabilties above a threshold,
#' order the states according to decreasing emission probabilties.
markstate_order <- function(emissions, mark, threshold=0.05){
  emissions %>%
    filter_(
      sprintf("markname == '%s'", mark),
      sprintf("prob > '%s'", threshold)
      ) %>%
    arrange(desc(prob)) %>%
    {.$state} %>%
    as.character()
}

# -------------------------------------------------------------------------
#' Crude ChomHMM state ordering

#' Crude ChromHMM state ordering according to Roadmap
#' @param model output of \code{load_model()}
#' @export
chromhmm_orderstates <- function(model){
  emissions <- model$emissions
  num_states <- model$num_states

  histone_marks <-
    c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3", "H3K9me3", "H3K27me3")
  threshold <- 0.05

  ## What states are emitted in what marks above a theshold
  marks_in_state <- list()
  for (mark in histone_marks) {
    marks_in_state[[mark]] <- markstate_order(emissions, mark, threshold)
  }

  ## States with following feature:

  # - TSS
  # - Active TSS (need to look at TSS neighborhood)
  # (will remove bivlalent states later)
  states_tss <- marks_in_state$H3K4me3

  ## Obtain list of enhancer states
  states_enh <-
    setdiff(c(marks_in_state$H3K27ac,
              marks_in_state$H3K4me1) %>%
              unique(),
            marks_in_state$H3K4me3)

  # - Bivalent
  states_biv <-
    intersect(
      c(marks_in_state$H3K9me3, marks_in_state$H3K27me3),
      c(states_tss, states_enh))

  ## Remove bivalent states from tss and enh list
  states_tss <-
    setdiff(states_tss, states_biv)
  states_enh <-
    setdiff(states_enh, states_biv)

  # - Enh genic
  # - Enh active
  # - Enh
  states_enhG <-
    intersect(states_enh, marks_in_state$H3K36me3)
  states_enhA <-
    intersect(marks_in_state$H3K27ac, marks_in_state$H3K4me1) %>%
    {intersect(states_enh, .)} %>%
    {setdiff(., states_enhG)}
  states_enhW <-
    setdiff(states_enh, c(states_enhG, states_enhA))

  # - Txn
  # - ZNF (need to look at enrichment overlap)
  states_txn <-
    setdiff(marks_in_state$H3K36me3,
            c(states_tss, states_enh, states_biv))

  ## Repressed states
  # - Heterochromatin
  # - PolyComb
  states_het <-
    setdiff(marks_in_state$H3K9me3, c(states_txn, states_biv))
  states_pc <-
    setdiff(marks_in_state$H3K27me3, c(states_txn, states_biv))
  states_het2 <-
    intersect(states_het, states_pc)

  states_repr <-
    c(
      setdiff(states_het, states_het2),
      states_het2,
      setdiff(states_pc, states_het2)
    )

  # - low
  states_low <-
    setdiff(1:num_states, unlist(marks_in_state))

  (stateorder <-
    c(states_tss, states_txn, states_enhG,
      states_enhA, states_enhW, states_biv,
      states_repr, states_low)
  )
}

#  ------------------------------------------------------------------------
#' Output state order as ChromHMM Reorder's stateorderingfile
#' @export
chromhmm_makestateorderingfile <- function(stateorder, filename=NULL){
  num_states <- length(stateorder)
  d <- data.frame(old_state=stateorder, new_state=1:num_states)
  if (!is.null(filename)) readr::write_tsv(d, filename, col_names=F) else d
}

#  ------------------------------------------------------------------------
#' Recreate ChromHMM emissions & transitions plot
#' @param input Emissions or transitions output of \code{chromhmm_loadmodel}
#' @param type Specify what type of data the input is:
#'             \code{e} for emissions or
#'             \code{t} for transitions
#' @param stateorder character list specifying the state order
#' @import ggplot2
#' @export
chromhmm_plt <- function(input, type, stateorder=NULL){
  if (type=="e") {
    plt <-
      ggplot(input, aes(x=markname,
                        y=state,
                        fill=prob)) +
      geom_tile() +
      scale_fill_continuous(low="white", high="blue") +
      xlab("") +
      ylab("") +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

  } else if (type=="t"){
    plt <-
      ggplot(input, aes(x=to, y=from, fill=prob)) +
      geom_tile() +
      scale_fill_continuous(low="white", high="blue") +
      xlab("State To") +
      ylab("State From")
  } else {
    stop("'type' must be one of: 'e' or 't'")
  }

  if (!is.null(stateorder)) plt + ylim(stateorder) else {
    num_states <- dplyr::last(input[[1]])
    plt + ylim(as.character(num_states:1))
  }
}
