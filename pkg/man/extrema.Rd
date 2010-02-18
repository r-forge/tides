\name{extrema}
\alias{extrema}
\title{Extrema}
\description{Calculate maxima and minima of quasi periodic time series. The function works good when the time series is 'smooth enough'. The function will fail when multiple local extrema exist, with a magnitude similar to the extremum extremorum and in an interval larger than [t-T2, t+T2] around the extremum extremorum at time t. When this is the case, this might be solved by adding extra wrappers around the function.}
\usage{extrema(h, h0, T2 = 5*60*60, hoffset = 3, filtconst = 1)}
\arguments{
  \item{h }{Water level time series. data frame with time and h column}
  \item{h0 }{Reference level, either single valued vector with dimension corresponding to h}
  \item{T2 }{'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h}
  \item{hoffset}{Offset level, to prevent spurious maxima generation due to small fluctuations}
  \item{filtconst }{Filtering constant for smoothing the time series}
}
\value{
  a list containing:
  \item{HL }{Data frame with extrema}
  \item{h }{original water level data frame with additional attributes}
}
\author{Tom Cox <tom.cox@ua.ac.be>}
\keyword{utilities}

