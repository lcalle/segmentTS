% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentTS.4segdist.R
\name{segmentTS.4segdist}
\alias{segmentTS.4segdist}
\title{Segment Distance}
\usage{
segmentTS.4segdist(obs.seg, sim.seg)
}
\arguments{
\item{obs.seg}{data.frame object with variables derived from segmentTS.1matchsignal}

\item{sim.seg}{data.frame object, variables as in obs.seg, but for simulated data}
}
\value{
list object with 4 outputs: time-series of the matching times and values (poly_t,poly) and the distance statistics (dist_tdiff,dist_vdiff)
}
\description{
This function takes in the full time-series data for observed and simulated from segmentTS.1matchsignal.
The boundaries of the segments in the time-series are defined using the
positions of the major peaks and troughs from fn segmentTS.3eqsignal(). Each segment
is then passed to segmentTS.4segdist to calculate the point-by-point segment statistics.
Matches similar signals in the two time-series. Distance statistics are simulation minus observation.
}
