% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_colorstrip.R
\name{plot_colorstrip}
\alias{plot_colorstrip}
\title{Add color to ggplot's facet strip.background}
\usage{
plot_colorstrip(d, p, statepalette, side = "top")
}
\arguments{
\item{d}{dataframe of the data used to create \code{p}.
There is at least 3 columns: x, y, and "\code{state}"
\code{state} is used to determine the number of strips
and should be equal to the number of facets of \code{p}.}

\item{p}{ggplot's base plot}

\item{statepalette}{character list where each element is a color and
the name of each element is a \code{state} level.}

\item{side}{Where is the facet strip? One of "top" or "right".}
}
\description{
Add color to ggplot's facet strip.background
}
\examples{
library(ggplot2)
library(dplyr)
d <- iris \%>\%
  mutate(
    state = Species,
    samples = Sepal.Length,
    bins = Sepal.Width)
p <- d \%>\%
  ggplot(aes(x=samples, y=bins)) +
  geom_point() +
  facet_wrap(~state)
statepalette <- c("red", "yellow", "blue")
statepalette <- setNames(statepalette, levels(d$state))
plot_colorstrip(d, p, statepalette)
}

