% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_weight.R
\name{get_phi_space_wts}
\alias{get_phi_space_wts}
\title{Calculate spatially-biased re-sampling weights with phi look-up}
\usage{
get_phi_space_wts(
  df,
  phi,
  rel.thr,
  sites_sf = NULL,
  lon = NULL,
  lat = NULL,
  crs = NULL,
  site_id,
  time
)
}
\arguments{
\item{df}{a data frame}

\item{phi}{a vector of the relevance values of 
\code{df}'s target variable}

\item{rel.thr}{a relevance threshold above which an 
observation is considered relevant}

\item{sites_sf}{An sf obejct containing station and IDs and 
geometry points of the locations. As an alternative, provide
\code{lon}, \code{lat}, and \code{crs}}

\item{lon}{the name of the column containing the location's longitude}

\item{lat}{the name of the column containing the location's latitude}

\item{crs}{the code for the Coordinate Reference System}

\item{site_id}{the name of the column containing location IDs}

\item{time}{the column name of the time-stamp}
}
\value{
A vector of spatially-biased re-sampling weights, scaled
to fit within range [0,1].
}
\description{
Calculate weights for re-sampling with a spatial bias with phi look-up.
Observations have a distance that tends to 1 as 
they are farther away from the closest relevant case (besides itself)
at time slice \code{t} (meaning they are more likely to be kept).
Farthest away from relevant cases at time slice t: d -> 1.
}
