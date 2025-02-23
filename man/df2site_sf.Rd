% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{df2site_sf}
\alias{df2site_sf}
\title{Create an sf object of available sites}
\usage{
df2site_sf(df, site_id, lon, lat, crs)
}
\arguments{
\item{df}{a data frame of the data set}

\item{site_id}{the name of the column containing location IDs}

\item{lon}{the name of the column containing the location's longitude}

\item{lat}{the name of the column containing the location's latitude}

\item{crs}{the code for the Coordinate Reference System}
}
\value{
a sf object, containing the geographic information for
each location in \code{df}
}
\description{
Extracts the location information from a data frame
and transforms into a \code{sf} object.
}
\seealso{
\code{\link[sf]{st_as_sf}}
}
