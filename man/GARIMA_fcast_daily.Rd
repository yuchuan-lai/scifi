\name{GARIMA_fcast_daily}
\alias{GARIMA_fcast_daily}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain daily temperature and precipitation simulations with the G-ARIMA model
}
\description{
Use the G-ARIMA model to simulate the future daily maximum and minimum temperature and precipitaiton series for a city based on the compiled long-term location-specific records and the future simulations from global climate models (GCMs). This is an integrative technique combining the ARIMA forecasting (as \code{\link{fcast_annual}}) and GCM climate projections.
}
\usage{
GARIMA_fcast_daily(hist.obs.daily, model.monthly.proj = "default", num.fcast.yr = "default", fcast.starting.yr = 2020, num.sets = 1, extra.days.for.bt = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.obs.daily}{
A dataframe object that includes the daily dates, maximum temperature, minimum temperature, and precipitation as four columns. It can be a dataframe object that is obtained using \code{\link{download}} function and with `daily` argument.
}
  \item{model.monthly.proj}{
A list object that includes the monthly average temperature (averged daily Tmax, Tmin, and Tmean) and monthly total precipitation projections from different GCMs. Two example datasets CMIP5.RCP4.5.monthly.US.proj and CMIP5.RCP8.5.monthly.US.proj are included in the package as an example. By default, these example model projections, which represent the averaged projections across the contiguous United States, are used. For improved projecting accuracy, the location specific climate projections should be used.
}
  \item{num.fcast.yr}{
A number of forecasting years. By default, the number of forecasting years is determined by macthing the provided climate projections (or the example projections).
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
Background information about the G-ARIMA model: Lai, Y., & Dzombak, D. A, 2021. Use of Integrated Global Climate Model Simulations and Statistical Time Series Forecasting to Project Regional Temperature and Precipitation. Journal of Applied Meteorology and Climatology.

Background information about the statistical forecasting models applied can be found at: Lai, Y., & Dzombak, D. A, 2020. Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

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

Lai, Y., and Dzombak, D. A., 2021: Use of Integrated Global Climate Model Simulations and Statistical Time Series Forecasting to Project Regional Temperature and Precipitation. Journal of Applied Meteorology and Climatology.
}
\author{
Yuchuan Lai
}

\seealso{
See also \code{\link{fcast_daily} as the ARIMA forecasting}
}
\examples{
# Download the historical daily data for Pittsburgh
pit.daily <- download("Pittsburgh", "daily")

# Obtain 1 set of G-ARIMA simulations of daily temperature and precipitation
# starting from 2020 in Pittsburgh (projection period: 2020-2099)
pit.daily.simu <- GARIMA_fcast_daily(pit.daily)

}
