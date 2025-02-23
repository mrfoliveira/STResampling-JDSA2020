% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resampling_utils.R
\name{calculate_relev}
\alias{calculate_relev}
\title{Calculate phi control}
\usage{
calculate_relev(form, dat, rel, cf = 1.5)
}
\arguments{
\item{form}{A formula describing the prediction problem.}

\item{dat}{A data frame containing the original imbalanced data set.}

\item{rel}{The relevance function which can be automatically ("auto") 
determined (the default) or may be provided by the user through a matrix
with the interpolating points.}

\item{cf}{Parameter needed if rel = 'auto'. The default is 1.5.}
}
\value{
A list containing a vector s.y containing sorted target values,
a vector y.relev containing the relevance values of s.y, and phi.control
}
\description{
Calculate phi control
}
