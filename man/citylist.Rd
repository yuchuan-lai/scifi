\name{citylist}
\alias{citylist}

\title{
Retrive a list of US cities with available climate data
}
\description{
Function for obtaining the list of cities with historical long-term temperature and precipitation records (annual or daily).
}
\usage{
citylist()
}
%- maybe also 'usage' for other objects documented here.
\details{
Information about the historical records can be found at: Lai, Y., & Dzombak, D. A. (2019). Use of Historical Data to Assess Regional Climate Change. Journal of Climate, 32(14), 4299-4320. Additionally, data can be directly accessed at: https://kilthub.cmu.edu/projects/Use_of_historical_data_to_assess_regional_climate_change/61538
}
\value{
A dataframe object that contains the information about the available cities, locations (latitude and longitude), period of record, and name of the weather stations (the raw data was collected).
}
\references{
Lai, Y., and D. A. Dzombak, 2019: Use of Historical Data to Assess Regional Climate Change. J. Clim., 32, 4299–4320, doi:10.1175/JCLI-D-18-0630.1. http://journals.ametsoc.org/doi/10.1175/JCLI-D-18-0630.1.

Owen, T. W., K. Eggleston, A. Degaetano, and R. Leffler, 2006: Accessing NOAA daily temperature and precipitation extremes based on combined/threaded station records. ftp://ftp.ncdc.noaa.gov/pub/data/papers/200686ams12.4tofree.pdf.
}
\author{
Yuchuan Lai
}
\note{
Daily maximum and minimum temperature and daily precipitation were acquired from Applied Climate Information System (ACIS), developed by the NOAA Northeast Regional Climate Center (NRCC).

The historical observations from ACIS belong to Global Historical Climatological Network - daily (GHCN-D) datasets.

The included stations were based on NRCC’s “ThreadEx” project, which combined daily temperature and precipitation extremes at 255 NOAA Local Climatological Locations, representing all large and medium size cities in U.S. (Owen et al. 2006).
}
\examples{
cities.info <- citylist()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% __ONLY ONE__ keyword per line
