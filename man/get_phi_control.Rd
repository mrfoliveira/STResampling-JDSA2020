% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resampling_utils.R
\name{get_phi_control}
\alias{get_phi_control}
\title{Get phi control}
\usage{
get_phi_control(y, rel, cf = 1.5)
}
\arguments{
\item{y}{the target variable}

\item{rel}{The relevance function which can be automatically ("auto") 
determined (the default) or may be provided by the user through a matrix
with the interpolating points.}

\item{cf}{Parameter needed if rel = 'auto'. The default is 1.5.}
}
\value{
phi.control object
}
\description{
Get phi control
}
