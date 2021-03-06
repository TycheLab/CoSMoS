% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratifyData.R
\name{stratifySeasonData}
\alias{stratifySeasonData}
\title{Stratify timeseries by season}
\usage{
stratifySeasonData(TS, season)
}
\arguments{
\item{TS}{time series}

\item{season}{name of the season}
}
\description{
Stratify timeseries by season
}
\examples{

x <- data.frame(date = seq(Sys.Date(), by = 'day', length.out = 1000),
                value = rnorm(1000))

stratifySeasonData(x, 'month')
}
\keyword{internal}
