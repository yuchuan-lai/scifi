\name{plot_annual}
\alias{plot_annual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot regional climate forecasts of annual indices from statistical forecasting model
}
\description{
Plot the temperature or precipitation forecasts of annual indices obtained using statistical forecasting models and \code{\link{fcast_annual}}.
}
\usage{
plot_annual(fcast.output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcast.output}{
A dataframe object that is produced by the \code{\link{fcast_annual}} as the output of statistical forecasting in annual resolution.
}
}
\details{
Historical observations are presented as points, while historical fitting is shown as the line during the historical period of record. Point forecasts and interval forecasts (80th and 95th intervals) are presented as the line and shaded are in the forecasting period.
}
\value{
A ggplot object.
}
\references{
Hyndman, R. J., and G. Athanasopoulos, 2018: Forecasting: principles and practice. OTexts,. ——, and Coauthors, 2018: forecast: Forecasting functions for time series and linear models.

Lai, Y., and Dzombak, D. A., 2020: Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

Lai, Y., and Dzombak, D. A., 2021: Use of Integrated Global Climate Model Simulations and Statistical Time Series Forecasting to Project Regional Temperature and Precipitation. Journal of Applied Meteorology and Climatology.
}
\author{
Yuchuan Lai
}
\seealso{
See Also as \code{\link{fcast_annual}}, \code{\link{plot_annual_cl}}, and \code{\link{plot_daily}}
}
\examples{
# Download the annual data for Pittsburgh
pit.annual <- download("Pittsburgh", "annual")

# Process the annual data of Pittsburgh to obtain annual average temperature
pit.temp <- select(pit.annual, "Avg.Temp")

# Obtain 20-year forecasts of annual average temperature
# starting from 2020 in Pittsburgh using the ARIMA model
pit.fcast.temp <- fcast_annual(pit.temp)

# Plot the obtained 20-year forecasts of annual average temperature
plot_annual(pit.fcast.temp)

# Obtain 15-year forecasts of annual maximum 1-day precipitation
# starting from 2006 in Chicago using the non-stationary GEV model

chi.annual <- download("Chicago", "annual")
chi.prcp.ex <- select(chi.annual, "Max.1.day.P")
chi.fcast.prcp.ex <- fcast_annual(chi.prcp.ex, 15, 2006, "gev")
plot_annual(chi.fcast.prcp.ex)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
