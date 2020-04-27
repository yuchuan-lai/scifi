\name{plot_annual_cl}
\alias{plot_annual_cl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot annual forecasts with estimated confidence levels of return periods
}
\description{
Plot the forecasts of temperature or precipitation extremes with estimated return periods obtained using the ARIMA forecasting model from \code{\link{fcast_annual_cl}}.
}
\usage{
plot_annual_cl(fcast.output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcast.output}{
A dataframe object that is produced by the \code{\link{fcast_annual_cl}} as the output of forecasting annual extremes and estimating of return periods.
}
}
\details{
Historical observations are presented as points, the initial estimates of historical and forecasted return periods are presented as lines, while shaded area indicates the confidence levels for the given return period.

Further bckground information about the statistical forecasting models applied can be found at: Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.
}
\value{
A ggplot object.
}
\references{
Lai, Y., & Dzombak, D. A. (in press). Use of the Autoregressive Integrated Moving Average (ARIMA) Model to Forecast Near-term Regional Temperature and Precipitation. Weather and Forecasting.

Pascual, L., J. Romo, and R. Esther, 2004: Bootstrap predictive inference for ARIMA processes. J. Time Ser. Anal., 25, 449–465.

Pascual, L., J. Romo, and E. Ruiz, 2005: Bootstrap prediction intervals for power-transformed time series. Int. J. Forecast., 21, 219–235, doi:10.1016/j.ijforecast.2004.09.006.
}
\author{
Yuchuan Lai
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{fcast_annual_cl}}, \code{\link{plot_annual_cl}}, and \code{\link{plot_daily}}
}
\examples{
# Plot the historical and forecasted 5-year return period in
# annual maximum 1-day precipitation at Chicago, with forecasting
# starting from 2020 and 20-year forecasting period

chi.annual <- download("Chicago", "annual")
chi.prcp.ex <- select(chi.annual, "Max.1.day.P")
chi.fcast.prcp.5yr <- fcast_annual_cl(chi.prcp.ex)
plot_annual_cl(chi.fcast.prcp.5yr)

# Alternatively, plot the 10-year return period with 99 percent confidence level starting from 2000
chi.fcast.prcp.10yr <- fcast_annual_cl(chi.prcp.ex, fcast.starting.yr = 2000, return.period = 10, confidence.level = 0.99)
plot_annual_cl(chi.fcast.prcp.10yr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
