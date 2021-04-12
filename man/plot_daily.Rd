\name{plot_daily}
\alias{plot_daily}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot daily temperature and precipitation forecasts from the ARIMA daily simulation model
}
\description{
Plot the daily temperature or precipitation forecasts obtained using the ARIMA daily simulation model in \code{\link{fcast_daily}}. Probability density functions of daily series and annual average calculations from the daily simulations are plotted.
}
\usage{
plot_daily(fcast.output, plot.var.name = "temp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcast.output}{
A list object that is produced by the \code{\link{fcast_daily}} as the output of ARIMA simulations of future temperature and precipitation.
}
  \item{plot.var.name}{
A character specifying the variable name for plotting. `temp`, `tmax`, `tmin`, and `prcp` are the acceptable arguments.
}
}
\details{
Blue colored lines are the individual simulated series while black line in pdf plot and black points in annual plot represent the historical observations.

Background information about the statistical forecasting models applied can be found at: Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

}
\value{
A ggplot object, including a probability density function plot and aggregated annual average of daily simulaitons plot. Blue colored lines are the individual simulated series while black line in pdf plot and black points in annual plot represent the historical observations.
}
\references{
Lai, Y., and Dzombak, D. A., 2021: Use of Integrated Global Climate Model Simulations and Statistical Time Series Forecasting to Project Regional Temperature and Precipitation. Journal of Applied Meteorology and Climatology.

Lai, Y., and Dzombak, D. A., 2020: Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.
}
\author{
Yuchuan Lai
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{fcast_daily}}, \code{\link{plot_annual}}, and \code{\link{plot_annual_cl}}
}
\examples{
# Download the historical daily data for Pittsburgh
pit.daily <- download("Pittsburgh", "daily")

# Obtain 1 set of 20-year simulations of daily temperature and precipitation
# starting from 2020 in Pittsburgh
pit.daily.simu <- fcast_daily(pit.daily)

# Plot the obtained 1 set of simulations of daily maximum temperature
plot_daily(pit.daily.simu, "tmax")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
