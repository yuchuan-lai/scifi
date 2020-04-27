\name{download}
\alias{download}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download the available historical climate data for the U.S. cities
}
\description{
Function for downloading historical long-term temperature and precipitation records (annual or daily) for a U.S. city. The list of cities with available climate records can be obtained by \code{\link{citylist}} function.
}
\usage{
download(city.name, record.type = "annual")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{city.name}{
Character input for the name of the city. The name of city follows the names used in the list of cities.
}
  \item{record.type}{
`annual` or `daily` as records of annual temperature and precipitation indices or daily temperature and precipitation observations for the selected cities
}
}
\details{
Annual average of daily mean temperature, Tmax, and Tmin are included for annual average temperature calculations. Calculations were based on the compiled daily temperature and precipitation records at individual cities.

Temperature and precipitation extreme indices include: warmest daily Tmax and Tmin, coldest daily Tmax and Tmin, warm days and nights, cold days and nights, maximum 1-day precipitation, maximum consecutive 5-day precipitation, precipitation amounts from very wet days.

Number of missing daily Tmax, Tmin, and precipitation values are included for each city.

Further information about the historical records can be found at: Lai, Y., & Dzombak, D. A. (2019). Use of Historical Data to Assess Regional Climate Change. Journal of Climate, 32(14), 4299-4320. Additionally, data can be directly accessed at: https://kilthub.cmu.edu/projects/Use_of_historical_data_to_assess_regional_climate_change/61538
}
\value{
A dataframe object that contains the years, missing daily data, and annual 13 temperature and precipitation index calculations for the selected city (with `annual` argument) or dates and daily maximum and minimum temperature and precipitaiton in the selected city (with `daily` argument)
}
\references{
Lai, Y., and D. A. Dzombak, 2019: Use of Historical Data to Assess Regional Climate Change. J. Clim., 32, 4299–4320, doi:10.1175/JCLI-D-18-0630.1. http://journals.ametsoc.org/doi/10.1175/JCLI-D-18-0630.1.

Owen, T. W., K. Eggleston, A. Degaetano, and R. Leffler, 2006: Accessing NOAA daily temperature and precipitation extremes based on combined/threaded station records. ftp://ftp.ncdc.noaa.gov/pub/data/papers/200686ams12.4tofree.pdf.
}
\author{
Yuchuan Lai
}
\note{
The historical observations from ACIS belong to Global Historical Climatological Network - daily (GHCN-D) datasets. The included stations were based on NRCC’s “ThreadEx” project, which combined daily temperature and precipitation extremes at 255 NOAA Local Climatological Locations, representing all large and medium size cities in U.S. (Owen et al. 2006).
}

\examples{
# Obtain historical annual data for Pittsburgh
pit.annual <- download("Pittsburgh", "annual")

# Obtain historical daily data for Pittsburgh
pit.daily <- download("Pittsburgh", "daily")
}

\keyword{datasets}
