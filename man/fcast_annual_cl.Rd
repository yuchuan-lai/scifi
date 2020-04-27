\name{fcast_annual_cl}
\alias{fcast_annual_cl}

\title{
Obtain ARIMA climate forecasts of annual extremes with estimated confidence levels for return periods
}
\description{
Estimate the confidence levels for different return periods when forecasting annual temperature and precipitaiton index series using the ARIMA model.
}
\usage{
fcast_annual_cl(hist.obs, num.fcast.yrs = 20, fcast.starting.yr = 2020, return.period = 5, num.bt = 100, confidence.level = 0.9)
}

\arguments{
  \item{hist.obs}{
A dataframe object that includes the years and one climate variable as two columns. It can be a datafram object that is obtained using \code{\link{select}} function. Use of annual maxima or minima series is intended, as confidence levels for different return periods for these extreme series are important climatic design values in engineering design.
}
  \item{num.fcast.yrs}{
A number of forecasting years. Number of forecasting years is recommended to be less than 20 years as long-term forecasts are unreliable because of climate variability.
}
  \item{fcast.starting.yr}{
A number specifying the starting year for forecasting (e.g., 2006). The forecasting starting year can be within the period of record to evaluate forecasting accuracy with observations.
}
  \item{return.period}{
A number specifying the desired return period for the estimation (e.g., 5).
}
  \item{num.bt}{
A number specifying the number of times for bootstrapping when estimating confidence levels. A relative large number is recommended for a more accuate estimation of confidence levels.
}
  \item{confidence.level}{
A number specifying specifying the level of statistical confidence. The number has to be bigger than 0 and smaller than 1. A default value of 0.9 indicates the 90 percent confidence level.
}
}
\details{
Background information about the statistical forecasting models applied can be found at: Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

Further information about the use of bootstrapping technique can be found at: 1) Pascual, L., J. Romo, and R. Esther, 2004: Bootstrap predictive inference for ARIMA processes. J. Time Ser. Anal., 25, 449–465. 2) Pascual, L., J. Romo, and E. Ruiz, 2005: Bootstrap prediction intervals for power-transformed time series. Int. J. Forecast., 21, 219–235, doi:10.1016/j.ijforecast.2004.09.006.
}
\value{
A dataframe object that contains the information about the estimated historical and forecasted return period with the given confidence levels. If the forecasting starting year is selected within the historical period of record, observations that can be used for validation of forecasts are also provided.
}
\references{
Hyndman, R. J., and G. Athanasopoulos, 2018: Forecasting: principles and practice. OTexts,. ——, and Coauthors, 2018: forecast: Forecasting functions for time series and linear models.

Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

Pascual, L., J. Romo, and R. Esther, 2004: Bootstrap predictive inference for ARIMA processes. J. Time Ser. Anal., 25, 449–465.

Pascual, L., J. Romo, and E. Ruiz, 2005: Bootstrap prediction intervals for power-transformed time series. Int. J. Forecast., 21, 219–235, doi:10.1016/j.ijforecast.2004.09.006.
}
\author{
Yuchuan Lai
}

\seealso{
See Also as \code{\link{fcast_annual}} and \code{\link{fcast_daily}}
}
\examples{
# Estimate the historical and forecasted 5-year return period in
# annual maximum 1-day precipitation at Chicago, with forecasting
# starting from 2020 and 20-year forecasting period

# Download the historical annual data for Chicago
chi.annual <- download("Chicago", "annual")

# Select the annual maximum 1-day precipitation series
chi.prcp.ex <- select(chi.annual, "Max.1.day.P")

# Forecast the annual maxima series and estimate the confidence levels
# with 5-year return period
chi.fcast.prcp.5yr <- fcast_annual_cl(chi.prcp.ex)

# Alternatively, estimat the 10-year return period with 99 percent
# confidence level starting from 2000
chi.fcast.prcp.10yr <- fcast_annual_cl(chi.prcp.ex, fcast.starting.yr = 2020, return.period = 10, confidence.level = 0.99)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
