#' Make sample pairs
#'
#' From a character list of samples,
#' create a dataframe specifying the pairs of samples
#' to compare for sample concordance.
#'
#' @param the_samples character list of sample ids
#' @examples
#' the_samples <- c("CEMT_33", "CEMT_34", "CEMT_50")
#' (sample_pairs <- make_sample_pairs(the_samples))
make_sample_pairs <- function(the_samples){
  sample_pairs <-
    combn(the_samples, 2) %>%
    t()
  colnames(sample_pairs) <- paste0("s", 1:2)
  (sample_pairs <-
     sample_pairs %>%
     tbl_df() %>%
     group_by(s1,s2)
  )
}

#  ------------------------------------------------------------------------
#' Load ChromHMM transition probabilities
#'
#' @param filename The ChromHMM model transitions file
#' @import dplyr
# filename <- file.path(chromhmmdir, "LearnModel_ordered/model/transitions_14.txt")
# (transition_probs <- load_transitiondb(filename))
load_transitiondb <- function(filename){
  transition_probs <- readr::read_tsv(filename)
  colnames(transition_probs)[1] <- "state_from"
  transition_probs <-
    transition_probs %>%
    tidyr::gather(state_to, prob, -state_from)
  transition_probs
}

#  ------------------------------------------------------------------------
#' Load unionbedg
#'
#' @param filename TSV file containing the unionbedg of ChromHMM sample segments.
#'                 The file contains the following columns:
#'                 \code{chrom}, \code{start}, \code{end}, <sampleID ...>;
#'                 Note: header is included in the file.
#' @param binsize The integer bin size used in ChromHMM (default = 200).
#' @details \code{filename} is created using bedtools \code{unionbedg}.
#' @import readr
#' @export
# filename <- file.path(chromhmmdir, "makeconsensus/CEMTcolonNormal/CEMTcolonNormal_unionbedg.txt")
# (unionbedg <- load_unionbedg(filename))
load_unionbedg <- function(filename, binsize=200){
  dat <-
    readr::read_tsv(
      filename,
      col_types = cols(chrom = col_character()),
      progress = FALSE
    ) %>%
    mutate(nbins = (end-start)/binsize)
}

#  ------------------------------------------------------------------------
#' Create frequency matrix
#'
#' To calculate the pairwise sample concordance,
#' we first create a frequency matrix
#' by comparing/counting the number of bins in each pairs of states
#' (one from sample1, the other from sample2).
#'
#' @details
#' After we have the frequency matrix,
#' we then use transition probabilities to account for transitions between states.
create_freqmatrix <- function(unionbedg, sample1, sample2){
  x <-
    unionbedg %>%
    # for each (sample1:stateA, sample2:stateB) pair
    # count the number of state pairs
    group_by_(sample1, sample2) %>%
    summarize(freq = sum(nbins)) %>%

    # then produce a frequency matrix where
    # sample1 states are the rows and
    # sample2 states are the columns
    tidyr::spread_(sample2, "freq")
  rownames(x) <- NULL
  x <-
    x %>%
    as.data.frame() %>%
    tibble::column_to_rownames(sample1) %>%
    as.matrix()
  # replace NA with 0
  x[is.na(x)] <- 0
  x
}

#  ------------------------------------------------------------------------
#' Distance between 2 states
#'
#' Get the distance score from the transition probabilities between 2 states
#'
#' @param transition_probs output of \code{load_transitiondb(...)}
#' @param x state from sample1
#' @param y state from sample2
#' @param type Specify which probability to use:
#'             "\code{max}" (default) or "\code{avg}".
#'             Specifying a type not in the list will give you "max".
#'             See details.
#' @details
#' Let A and B be the states of sample1 and sample2 in a given bin.
#' Let P(.) be the transition probability of a state.
#' When A==B, the score returned is 0.
#' When type=="avg", the score returned is 1-mean(P(A->B), P(B->A)).
#' When type=="max", the score returned is 1-max(P(A->B), P(B->A)).
state_distance <- function(transition_probs, x, y, type="max"){
  if (x==y) {
    ## if states are equal, return 0
    0
  } else {
    ## else return 1 - transition_score
    ## where transition_score is
    ## mean(P(A->B), P(B->A)) when type=="avg" or
    ## max(P(A->B), P(B->A)) when type=="max"
    x <- sub("^U|^E", "", x)
    y <- sub("^U|^E", "", y)
    t_probs <-
      transition_probs %>%
      filter(
        (state_from==x&state_to==y) | (state_from==y&state_to==x)
      ) %>% {.$prob}
    if (type == "avg") {
      1-mean(t_probs)
    } else {
      1-max(t_probs)
    }
  }
}

#  ------------------------------------------------------------------------
#' Distance between 2 states
#'
#' Get the distance score from the transition probabilities between 2 samples
sample_distance <- function(freqmatrix, transition_probs, type="max"){
  freqmatrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample1") %>%
    tidyr::gather(sample2, count, -sample1) %>%
    group_by(sample1, sample2) %>%
    mutate(
      score = state_distance(transition_probs, sample1, sample2, type),
      tscore = count * score
    ) %>% {.$tscore} %>%
    sum()
}

#  ------------------------------------------------------------------------
#' Calculate sample concordance
#'
#' @param the_samples
#'           character list of sample ids used for the pairwise comparisons
#' @param filename.transitions
#'           The ChromHMM model transitions file
#' @param filename.unionbedg
#'           TSV file containing the unionbedg of ChromHMM sample segments.
#'           The file contains the following columns:
#'           \code{chrom}, \code{start}, \code{end}, <sampleID ...>;
#'           Note: header is included in the file.
#'@param verbose boolean; print messages?
#'
#' @details The output will be a dataframe containing 4 columns:
#'          \code{s1, s2, dist_max, dist_avg}. s1 and s2 refers to
#'          all pairwise combinations of \code{the_samples}.
#'          dist_max refers to the distance score using 1-max transition prob;
#'          dist_avg refers to the distance score using 1-mean transition prob.
#' @export
sample_concordance <- function(the_samples,
                               filename.transitions, filename.unionbedg,
                               verbose=FALSE){
  sample_pairs <- make_sample_pairs(the_samples)

  ## Load transition probs
  transition_probs <- load_transitiondb(filename.transitions)

  ## Load data
  unionbedg <- load_unionbedg(filename.unionbedg)

  ## For each sample pairs,
  ## create frequency matrix and
  ## compute distamce metrics per sample pair
  (dat.dist <-
    sample_pairs %>%
    mutate(
      dist_max = create_freqmatrix(unionbedg, s1,s2) %>%
        sample_distance(transition_probs, "max"),
      dist_avg = create_freqmatrix(unionbedg, s1,s2) %>%
        sample_distance(transition_probs, "avg")
    )
  )
}

#  ------------------------------------------------------------------------
#' Sample concordance summary
#'
#' Computes the mean, median, min, max, sd
#' of the sample pairs from \code{sample_concordance()}
#'
#' @param sampleconcordance the output of \code{sample_concordance()}
#' @param use_scale Boolean; Scale the data with \code{fun_scale}
#' @param fun_scale See \code{sample_concordance_vizH()}
#' @export
sample_concordance_summary <- function(
  sampleconcordance, use_scale = FALSE,
  fun_scale = function(x){1-(x/15181508)}){

  ## Prepare data
  dat <-
    sampleconcordance %>%
    tidyr::gather(stat, val, -starts_with("s"))
  if (use_scale) dat<-mutate(dat, val = fun_scale(val))

  ## Calculate summary statistics
  dat %>%
    group_by(stat) %>%
    summarize(
      mean = mean(val),
      median = median(val),
      min = min(val),
      max = max(val),
      sd = sd(val)
    )
}

#  ------------------------------------------------------------------------
#' Sample concordance visualization
#'
#' @param sampleconcordance the output of \code{sample_concordance()}
#' @param type See \code{state_distance()}
#' @param fun_scale Function to scale y values.
#'                  See \code{sample_concordance_vizH()}.
#'                  Specify "\code{NULL}" for no scaling
#' @param fun_scale.label y-axis label associated with \code{fun_scale}
#' @export
sample_concordance_viz <- function(sampleconcordance, type="max",
                                   fun_scale = function(x){1-(x/15181508)},
                                   fun_scale.label = "Similarity"){
  require(ggplot2)
  if (type=="avg") dmetric<-"dist_avg" else dmetric<-"dist_max"

  if (!is.null(fun_scale)) {
    sampleconcordance[[dmetric]] <- fun_scale(sampleconcordance[[dmetric]])
  } else if (is.null(fun_scale) && fun_scale.label == "Similarity") {
    fun_scale.label <- "Distance"
  }

  sampleconcordance %>%
    mutate(pair = paste(s1, s2, sep="-")) %>%
    ggplot() +
    aes_string(
      x=sprintf("reorder(pair, %s)", dmetric),
      y=dmetric
    ) +
    geom_point() +
    xlab("Pair of samples") +
    ylab(fun_scale.label) +
    theme(axis.text.x = element_text(angle=90, vjust=1, hjust=.5))
}

#  ------------------------------------------------------------------------
#' Sample concordance visualization (heatmap)
#'
#' Visualize sample concordance in terms of a heatmap
#'
#' @param sampleconcordance the output of \code{sample_concordance()}
#' @param type One of "\code{max}" (default) or "\code{avg}";
#'             see \code{state_distance()} for details.
#'             This parameter chooses which distance metric to use.
#' @param fun_scale function to scale matrix;
#'                  use "\code{NULL}" for no scaling
#'                  (i.e use the raw sample distances.
#'                  See details.
#' @param ... additional arguments to \code{pheatmap(...)}
#'
#' @details
#'   The value 15181508 is the number of 200 bp bins as segmented by ChromHMM.
#'   x/15181508 accounts for the total distance (i.e. error) between samples.
#'   1-(x/15181508) transforms the distance to a similarity measure.
#'   For scaling the distance between 0 and 1, use
#'   \code{function(x){(x-min(x,na.rm=T)) / (max(x,na.rm=T)-min(x,na.rm=T))}}
#' @export
sample_concordance_vizH <- function(sampleconcordance, type="max",
                                    fun_scale = function(x){1-(x/15181508)},
                                    ...){
  if (type=="avg") dmetric<-"dist_avg" else dmetric<-"dist_max"

  ## Create matrix of pairwise symmetric distances
  mat <-
    rbind(
      sampleconcordance,
      dplyr::rename(sampleconcordance, s1=s2, s2=s1)
    ) %>%
    select_("s1", "s2", dist=dmetric) %>%
    tidyr::spread(s2, dist)
  rownames(mat) <- NULL
  mat <-
    as.data.frame(mat) %>%
    tibble::column_to_rownames(var = "s1")

  ## Scale matrix
  if (!is.null(fun_scale)) mat<-fun_scale(mat)

  ## Create heatmap
  pheatmap::pheatmap(mat, ...)

}
