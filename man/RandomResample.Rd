% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resampling_method.R
\name{RandomResample}
\alias{RandomResample}
\title{Random re-sampling for imbalanced regression problems (with spatio-temporal bias)}
\usage{
RandomResample(
  form,
  dat,
  type,
  C.perc,
  thr.rel,
  rel = "auto",
  cf = 1.5,
  repl = ifelse(type == "under", FALSE, TRUE),
  pert = 0.1,
  time = NULL,
  site_id = NULL,
  bias = FALSE,
  ...
)
}
\arguments{
\item{form}{A formula describing the prediction problem.}

\item{dat}{A data frame containing the original imbalanced data set.}

\item{type}{Type of re-sampling to apply. Can be one of "under", "over", 
and "gauss", depending on whether the user wants to under-sample normal cases,
over-sample extreme cases or add Gaussian noise to replicated extreme cases.}

\item{C.perc}{Vector containing percentage values (or a single  value that 
will be used for all bumps. In under-sampling, C.perc of the size of each 
bump of normal values will be kept in the final data set. In the case of 
over-sampling and Gaussian noise, C.perc of the size of each bump of extreme
values will be added to the final data set. Bumps are ordered in ascending 
order of the target value.}

\item{thr.rel}{A number indicating the relevance threshold below which a 
case is considered as belonging to the normal "class".}

\item{rel}{The relevance function which can be automatically ("auto") 
determined (the default) or may be provided by the user through a matrix
with the interpolating points.}

\item{cf}{Parameter needed if rel = 'auto'. The default is 1.5.}

\item{repl}{A Boolean value controlling whether replication is allowed
when re-sampling observations. Defaults to FALSE when under-sampling
and to TRUE when over-sampling or adding Gaussian noise.}

\item{pert}{Standard deviation of gaussian noise as a percentage of
of each variable original standard deviation. Only necessary if type = "gauss"}

\item{time}{Column name of the time-stamp (if available). 
Only necessary if bias = TRUE or type = "gauss"}

\item{site_id}{Column containing location IDs (if available).
Only necessary if bias = TRUE or type = "gauss"}

\item{bias}{Boolean indicating whether spatio-temporal bias should be
factored in while re-sampling}

\item{...}{Parameters to feed to \code{sample_wts} in case bias = TRUE.}
}
\description{
Random re-sampling for imbalanced regression problems (with spatio-temporal bias)
}
\references{
Paula Branco, Rita P. Ribeiro, Luis Torgo (2016)., 
UBL: an R Package for Utility-Based Learning, 
CoRR abs/1604.08079 [cs.MS], URL: http://arxiv.org/abs/1604.08079
}
\seealso{
\code{\link[UBL]{RandUnderRegress}}, \code{link{sample_wts}}.
}
