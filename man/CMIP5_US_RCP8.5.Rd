\name{CMIP5.RCP8.5.monthly.US.proj}
\alias{CMIP5.RCP8.5.monthly.US.proj}
\docType{data}
\title{
Monthly temperature and precipitation projections of the U.S. from 28 selected GCMs of CMIP5 (RCP8.5)
}
\description{
Obtained monthly GCM temperature and precipitation projections averaged among the grid locations across the contiguous United States. The GCMs belong to CMIP5 and a total of 28 GCMs are included. The projecitons are under RCP8.5. Note that the projections were obtained for the grid locations within latitude: 25.12993 to 49.38323 and longitude: -124.6813 to -67.00742, which is slightly different from the projections of the contiguous United States.
}
\usage{data("CMIP5.RCP8.5.monthly.US.proj")}
\format{
  A list object with 4 subsets of dataframes as monthly averaged daily Tmax, daily Tmin, and daily Tmean (as temp), and monthly total precipitation.
  \describe{
    \item{\code{tmax, tmin, temp}}{a dataframe contains the dates (the first date of each month is used) and temperature projections (in unit of Fahrenheit) from each GCM (shown as column names)}
    \item{\code{prcp}}{a dataframe contains the dates (the first date of each month is used) and precipitation projections (in unit of inch) from each GCM (shown as column names)}
  }
}

\source{
The GCM simulations of CMIP5 were acquired from the Working Group on Coupled Modelling (WGCM) of World Climate Research Programme (WCRP) at: https://esgf-node.llnl.gov/projects/cmip5/.
}
\references{
Taylor, K. E., Stouffer, R. J., and Meehl, G. A., 2012. An overview of CMIP5 and the experiment design. Bulletin of the American Meteorological Society, 93(4), 485–498.
Lai, Y., and Dzombak, D. A., 2021: Use of Integrated Global Climate Model Simulations and Statistical Time Series Forecasting to Project Regional Temperature and Precipitation. Journal of Applied Meteorology and Climatology.
}
\examples{
data(CMIP5.RCP8.5.monthly.US.proj)
}
\keyword{datasets}
