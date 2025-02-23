% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_weight.R
\name{get_time_wts}
\alias{get_time_wts}
\title{Calculate temporally-biased re-sampling weights without phi look-up}
\usage{
get_time_wts(times, phi, rel.thr)
}
\arguments{
\item{times}{a vector of time-stamps}

\item{phi}{a vector of the relevance values of 
\code{df}'s target variable}

\item{rel.thr}{a relevance threshold above which an 
observation is considered relevant}
}
\value{
A vector of temporally-biased re-sampling weights, scaled
to fit within range [0,1].
}
\description{
Calculate weights for re-sampling with a temporal bias.
For each station, obsrvations that have been preceded by
extreme values the most recently, have higher weights
(meaning that observations where extreme values were recently
measured are more likely to be selected).
Most recent extreme observation in the past at station L: w -> 1; 
no rare observations in the past at station L: w -> 0.
}
\author{
Mariana Oliveira
}
