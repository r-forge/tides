\name{TidalCharacteristics}
\alias{TidalCharacteristics}
\title{Calculate tidal characteristics}
\description{Calculates the characteristics of observed tidal water levels. Wrapper of the functions \code{\link{extrema}}, \code{\link{IT}} and \code{\link{IF}}}
\usage{TidalCharacteristics(h,h0 = h$h0,h0marg = 0.3, T2 = 5*60*60,filtconst = 50, dtMax = 15,unit = "mins", Tavg = 12.4*60 ) 
TC <- TidalCharacteristics(h,...)}
\arguments{
  \item{h }{Water level time series. data frame with time and h column}
  \item{h0 }{Reference level, either single valued vector with dimension corresponding to h}
  \item{h0marg }{Margin on reference level, to cope with small fluctuations in the Water level time series}
  \item{T2 }{'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h}
  \item{filtconst }{Filtering constant for smoothing the time series}
  \item{dtMax}{maximum accepted time interval in a continuous series. Bigger time intervals are considered to be gaps}
  \item{unit}{unit of dtMax, Tavg}
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
TC <- TidalCharacteristics(waterlevels,h0=291)
TC
plot(TC)
}

