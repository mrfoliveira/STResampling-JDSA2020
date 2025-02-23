% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_framework.R
\name{estimates}
\alias{estimates}
\title{Estimate error using a chosen method}
\usage{
estimates(
  data,
  form,
  estimator = "kf_xval",
  est.pars = list(nfolds = 10, fold.alloc.proc = "Trand_SPrand"),
  workflow = "simple_workflow",
  wf.pars = NULL,
  evaluator = "evaluate",
  eval.pars = NULL,
  seed = 1234
)
}
\arguments{
\item{data}{a data frame}

\item{form}{a formula for learning}

\item{estimator}{the name of an error estimator function}

\item{est.pars}{a named list of arguments to feed to \code{estimator}}

\item{workflow}{the name of the workflow to use for making predictions}

\item{wf.pars}{a named list of arguments to feed to \code{workflow}}

\item{evaluator}{the name of the function to use to calculate evaluation results}

\item{eval.pars}{a named list of arguments to feed to \code{evaluator}}

\item{seed}{a seed to set before performing estimates}
}
\value{
The results of \code{evaluator} after applying \code{estimator} to the
learning task
}
\description{
Estimate error using a chosen method
}
