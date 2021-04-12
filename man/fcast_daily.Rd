\name{fcast_daily}
\alias{fcast_daily}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain daily temperature and precipitation simulations with the ARIMA model
}
\description{
Use the ARIMA model to simulate the future daily maximum and minimum temperature and precipitaiton series for a city based on the compiled long-term location-specific records. The simulations are incorporated with the seasonal forecasts obtained from the SARIMA model and can be seen as a realization of future daily temperature and precipitation.
}
\usage{
fcast_daily(hist.obs.daily, num.fcast.yr = 20, fcast.starting.yr = 2020, num.sets = 1, extra.days.for.bt = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.obs.daily}{
A dataframe object that includes the daily dates, maximum temperature, minimum temperature, and precipitation as four columns. It can be a dataframe object that is obtained using \code{\link{download}} function and with `daily` argument.
}
  \item{num.fcast.yr}{
A number of forecasting years. Number of forecasting years is recommended to be less than 20 years as long-term forecasts are unreliable because of climate variability.
}
  \item{fcast.starting.yr}{
A number specifying the starting year for forecasting (e.g., 2006). The forecasting starting year can be within the period of record to evaluate forecasting accuracy with observations.
}
  \item{num.sets}{
A number specifying the number of daily simulation sets. One set of simulation represents one realization of future daily temperature and precipitation. A larger number for simulation sets can lead to a longer time for running.
}
  \item{extra.days.for.bt}{
A number specifying the number of extra days selected from previous and succeeding months when bootstrapping. A larger number will increase the randomness as a result.
}
}
\details{
Background information about the statistical forecasting models applied can be found at: Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

Additional informaiton about the time series forecasting and bootstrap can be found at: Hyndman, R. J., and G. Athanasopoulos, 2018: Forecasting: principles and practice. OTexts,. ——, and Coauthors, 2018: forecast: Forecasting functions for time series and linear models.
}
\value{

\item{hist.obs}{Historical observations of daily temperature and precipitation}
\item{fcast.tmax}{Simulated future daily maximum temperature}
\item{fcast.tmin}{Simulated future daily maximum temperature}
\item{fcast.temp}{Simulated future daily mean temperature}
\item{fcast.prcp}{Simulated future daily precipitation}

}
\references{
Hyndman, R. J., and G. Athanasopoulos, 2018: Forecasting: principles and practice. OTexts,. ——, and Coauthors, 2018: forecast: Forecasting functions for time series and linear models.

Lai, Y., and Dzombak, D. A, 2020. Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.
}
\author{
Yuchuan Lai
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{fcast_annual}} and \code{\link{fcast_annual_cl}}
}
\examples{
# Download the historical daily data for Pittsburgh
pit.daily <- download("Pittsburgh", "daily")

# Obtain 1 set of 20-year simulations of daily temperature and precipitation
# starting from 2020 in Pittsburgh
pit.daily.simu <- fcast_daily(pit.daily)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
