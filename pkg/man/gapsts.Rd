\name{gapsts}
\alias{gapsts}
\title{Gaps in time series}
\description{Determine gaps (time intervals bigger than a certain maximum dtMax) in a time series}
\usage{gapsts(ts, dtMax,unit = "mins")}
\arguments{
  \item{ts }{array of times, consisting of different continuous subseries seperated by large gaps}
  \item{dtMax }{maximum time interval in a continuous series}
  \item{unit }{unit of dtMax}
}
\value{
  Data frame with the initial time, end time and time difference (unit = unit) of each interval > dtMax
}
\author{Tom Cox <tom.cox@ua.ac.be>}
\keyword{utilities}

