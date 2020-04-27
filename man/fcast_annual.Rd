\name{fcast_annual}
\alias{fcast_annual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain regional climate forecasts of annual indices with the statistical forecasting models
}
\description{
Use the ARIMA model or other common statistical techniques to forecast the annual temperature and precipitaiton index series, based on the compiled long-term location-specific records. Other statistical techniques include baseline, baseline (recent past 30-years), linear trend, and non stationary GEV methods.
}
\usage{
fcast_annual(hist.obs, num.fcast.yrs = 20, fcast.starting.yr = 2020, forecast.model = "ARIMA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.obs}{
A dataframe object that includes the years and one climate variable as two columns. It can be a datafram object that is obtained using \code{\link{select}} function.
}
  \item{num.fcast.yrs}{
A number of forecasting years. Number of forecasting years is recommended to be less than 20 years as long-term forecasts are unreliable because of climate variability.
}
  \item{fcast.starting.yr}{
A number specifying the starting year for forecasting (e.g., 2006). The forecasting starting year can be within the period of record to evaluate forecasting accuracy with observations.
}
  \item{forecast.model}{
A character object specifying the statistical model used for forecasting. `ARIMA`, `linear`, `baselineall`, `baseline30`, and `gev` are the acceptable input arguments for the forecasting model, corresponding to the ARIMA, linear trend, baseline, baseline (recent 30-year), and non-stationary GEV methods.
}
}
\details{
Background information about the statistical forecasting models applied can be found at: Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

Further informaiton about the time series forecasting can be found at: Hyndman, R. J., and G. Athanasopoulos, 2018: Forecasting: principles and practice. OTexts,. ——, and Coauthors, 2018: forecast: Forecasting functions for time series and linear models.
}
\value{
A dataframe object that contains the information about the historical fitting, residuals, point forecasts, and interval forecasts (80th and 95th intervals) obtained from the forecasting model. If the forecasting starting year is selected within the historical period of record, observations that can be used for validation of forecasts are also provided.
}
\references{
Hyndman, R. J., and G. Athanasopoulos, 2018: Forecasting: principles and practice. OTexts,. ——, and Coauthors, 2018: forecast: Forecasting functions for time series and linear models.

Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.
}
\author{
Yuchuan Lai
}

\seealso{
See Also as \code{\link{fcast_annual_cl}} and \code{\link{fcast_daily}}
}

\examples{
# Download the annual data for Pittsburgh
pit.annual <- download("Pittsburgh", "annual")

# Process the annual data of Pittsburgh to obtain annual average temperature
pit.temp <- select(pit.annual, "Avg.Temp")

# Obtain 20-year forecasts of annual average temperature
# starting from 2020 in Pittsburgh using the ARIMA model
pit.fcast.temp <- fcast_annual(pit.temp)

# Obtain 15-year forecasts of annual maximum 1-day precipitation
# starting from 2006 in Chicago using the non-stationary GEV model

chi.annual <- download("Chicago", "annual")
chi.prcp.ex <- select(chi.annual, "Max.1.day.P")
chi.fcast.prcp.ex <- fcast_annual(chi.prcp.ex, 15, 2006, "gev")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
