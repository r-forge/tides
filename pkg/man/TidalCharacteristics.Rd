\name{TidalCharacteristics}
\alias{TidalCharacteristics}
\title{Calculate tidal characteristics}
\description{Calculates the characteristics of observed tidal water levels. Wrapper of the functions \code{\link{extrema}}, \code{\link{IT}} and \code{\link{IF}}}. Also works on time series with gaps.
\usage{TidalCharacteristics(h,h0 = h$h0,T2 = 5*60*60, hoffset = 3, filtconst = 1, dtMax = 15,unit = "mins", Tavg = 12.4*60 ) 
TC <- TidalCharacteristics(h,...)}
\arguments{
  \item{h }{Water level time series. data frame with time and h column}
  \item{h0}{Reference level, either single valued or vector with dimension corresponding to h}
  \item{T2 }{'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h}
  \item{hoffset}{Offset level, to prevent spurious maxima generation due to small fluctuations}
  \item{filtconst }{Filtering constant for smoothing the time series}
  \item{dtMax}{Maximum accepted time interval in a continuous series. Bigger time intervals are considered to be gaps}
  \item{unit}{Unit of dtMax, Tavg}
  \item{Tavg}{Average period of time series}
}
\value{
  An object of class \code{Tides}, i.e. a list containing:
  \item{HL }{Data frame with extrema}
  \item{h }{original water level data frame with additional attributes}
  \item{gaps}{a data frame containing start and end times of gaps in the series}
  \item{IF}{inundation frequency of the reference level}
  \item{ITs}{inundation times at the reference level}
  \item{DTs}{dry times at the reference level}
  \item{h0}{reference level}
  \item{N}{Total number of cycles in time span}
}
\seealso{
\code{\link{extrema}}, \code{\link{IT}}, \code{\link{plot.Tides}}
}
\author{Tom Cox <tom.cox@ua.ac.be>}
\keyword{utilities}
\examples{
TC <- TidalCharacteristics(waterlevels)
TC
plot(TC)
}

