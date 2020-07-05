# ARIMA daily temperature and precipitation simulation model, i.e., forecasts in daily resolutions:
dcast_data <- function(daily.data) {
  daily.data <- data.frame("Year" = format(daily.data$Date, "%Y"), daily.data)
  daily.data$Date <- format(daily.data$Date, "%m-%d")
  daily.data.dcast <- reshape2::dcast(daily.data, Date ~ Year, value.var = "Value")
  return(daily.data.dcast)
}
sampling_date <- function(current.fcast.yr, bt.hist.yr, jj, extra.days.for.bt) {
  select.year <- bt.hist.yr[jj]
  adj.days <- sample(c(- extra.days.for.bt:extra.days.for.bt),1)
  ini.select.date <- as.Date(sprintf("%02d-%02d-01", select.year, jj))
  select.st.date <- ini.select.date - adj.days
  if (lubridate::leap_year(current.fcast.yr) == TRUE & jj ==2) {
    select.days <- lubridate::days_in_month(jj) + 1
  } else {
    select.days <- lubridate::days_in_month(jj)
  }
  select.date <- data.frame("select.date" = seq(select.st.date, (select.st.date + select.days - 1), by ="day"),
                            "baseline.date" = seq(ini.select.date, (ini.select.date + select.days - 1), by ="day"))
  return(select.date)
}

#' @export
fcast_daily <- function(hist.obs.daily, num.fcast.yr = 20, fcast.starting.yr = 2020, num.sets = 1, extra.days.for.bt = 7) {

  #' @import forecast
  #' @import zoo

  #' @importFrom xts xts
  #' @importFrom xts apply.yearly
  #' @importFrom xts apply.monthly
  #' @importFrom lubridate day
  #' @importFrom lubridate month
  #' @importFrom lubridate year
  #' @importFrom lubridate days_in_month
  #' @importFrom stats ts
  #' @importFrom reshape2 dcast

  `%>%` <- dplyr::`%>%`

  # Note: Because of the possible missing daily data in the historical observations (i.e., NA values),
  # and historical residuals are used to bootstrap and create future simulation series,
  # the future simulations can have missing daily values (corresponding to the resampled NA values).
  # While this study has been using quality-checked historical data (no more than 10 daily data missing for any year),
  # interpolation of these missing values in the historical sereis prior to running the codes here
  # can be used to avoid missing data in future series,
  # or interpolation of missing values in the future series after simulations can be applied as well.

  # Custom input for the forecast starting year and number of forecasting years (years of simulation)
  # extra.days.for.bt <- 7
  colnames(hist.obs.daily) <- c("Date", "tmax", "tmin", "prcp")

  end.record.yr <- lubridate::year(hist.obs.daily$Date[nrow(hist.obs.daily)])
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yr - 1

  # Evaluate the number of missing data in the daily records:
  tmax.na.count <- xts::xts(is.na(hist.obs.daily$tmax), as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
  tmin.na.count <- xts::xts(is.na(hist.obs.daily$tmin), as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
  prcp.na.count <- xts::xts(is.na(hist.obs.daily$prcp), as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))

  tmax.na.count.annual <- xts::apply.yearly(tmax.na.count, FUN = sum)
  tmin.na.count.annual <- xts::apply.yearly(tmin.na.count, FUN = sum)
  prcp.na.count.annual <- xts::apply.yearly(prcp.na.count, FUN = sum)

  # Estimate the number of historical years that will be used
  city.hist.years.all <- lubridate::year(tmax.na.count.annual)[tmax.na.count.annual <= 10 & tmin.na.count.annual <= 10 & prcp.na.count.annual <= 10]
  st.yr <- city.hist.years.all[1]
  if (end.record.yr == (fcast.starting.yr - 1)) {
    hist.years <- city.hist.years.all
  } else {
    hist.years <- city.hist.years.all[! city.hist.years.all %in% c(fcast.starting.yr:end.record.yr)]
  }
  num.yr <- length(hist.years)

  # Extract the historical daily observations
  tmax.daily.data <- hist.obs.daily$tmax
  tmin.daily.data <- hist.obs.daily$tmin
  prcp.daily.data <- hist.obs.daily$prcp
  daily.dates.all <- as.Date(hist.obs.daily$Date)

  # Convert daily data into xts time series and extract the years of historical observaitons
  # (leave the years that will be used for cross-validation)
  # Also becasue of QA, there may be some years of record should be discarded
  tmax.daily.xts <- xts::xts(as.numeric(as.character(tmax.daily.data)), as.Date(daily.dates.all, "%Y-%m-%d"))
  tmax.daily.xts <- window(tmax.daily.xts, start = paste0(st.yr,"-01-01"), end = paste0((fcast.starting.yr - 1),"-12-31"))
  tmax.daily.xts <- tmax.daily.xts[paste0(city.hist.years.all)]
  tmin.daily.xts <- xts::xts(as.numeric(as.character(tmin.daily.data)), as.Date(daily.dates.all, "%Y-%m-%d"))
  tmin.daily.xts <- window(tmin.daily.xts, start = paste0(st.yr,"-01-01"), end = paste0((fcast.starting.yr - 1),"-12-31"))
  tmin.daily.xts <- tmin.daily.xts[paste0(city.hist.years.all)]
  prcp.daily.xts <- xts::xts(as.numeric(as.character(prcp.daily.data)), as.Date(daily.dates.all, "%Y-%m-%d"))
  prcp.daily.xts <- window(prcp.daily.xts, start = paste0(st.yr,"-01-01"), end = paste0((fcast.starting.yr - 1),"-12-31"))
  prcp.daily.xts <- prcp.daily.xts[paste0(city.hist.years.all)]

  # daily.dates store the dates for the extracted historical obsvervations (tmax is used here)
  daily.dates <- zoo::index(tmax.daily.xts)
  # daily.monthly.dates is variable for assigning the monthly baseline with the corresponding historical month
  daily.monthly.dates <- seq(as.Date(paste0(hist.years[1], "-01-01")), as.Date(paste0((fcast.starting.yr - 1), "-12-01")), by="month")
  daily.monthly.dates <- daily.monthly.dates[lubridate::year(daily.monthly.dates) %in% hist.years]

  # Daily temperature can also be calcualted
  temp.daily.xts <- (tmax.daily.xts + tmin.daily.xts) / 2

  # Create dataframes to store the historical daily data (log-transformation of precipitation data is performed)
  # The referred percipitation in the following steps are transfromed series unless specified otherwise
  # The "+ 1" in the log-transformation is to yield 0 for no precipitation events
  tmax.daily.df <- data.frame(Date = daily.dates, Value = as.numeric(as.character(tmax.daily.xts)))
  tmin.daily.df <- data.frame(Date = daily.dates, Value = as.numeric(as.character(tmin.daily.xts)))
  prcp.daily.log.xts <- log(prcp.daily.xts + 1)
  prcp.daily.log.df <- data.frame(Date = daily.dates, Value = as.numeric(as.character(prcp.daily.log.xts)))

  # Generate dataframes for all daily data for later use
  tmax.daily.all.df <- data.frame(Date = as.Date(hist.obs.daily$Date), Value = as.numeric(as.character(hist.obs.daily$tmax)))
  tmin.daily.all.df <- data.frame(Date = as.Date(hist.obs.daily$Date), Value = as.numeric(as.character(hist.obs.daily$tmin)))
  prcp.daily.all.df <- data.frame(Date = as.Date(hist.obs.daily$Date), Value = log(hist.obs.daily$prcp + 1))

  # Calculate the missing daily data for later use:
  tmax.daily.na.df <- data.frame(Date = as.Date(hist.obs.daily$Date), Value = is.na(tmax.daily.all.df$Value))
  tmin.daily.na.df <- data.frame(Date = as.Date(hist.obs.daily$Date), Value = is.na(tmin.daily.all.df$Value))
  prcp.daily.na.df <- data.frame(Date = as.Date(hist.obs.daily$Date), Value = is.na(prcp.daily.all.df$Value))
  tmax.daily.na.dcast <- dcast_data(tmax.daily.na.df)
  tmin.daily.na.dcast <- dcast_data(tmin.daily.na.df)
  prcp.daily.na.dcast <- dcast_data(prcp.daily.na.df)

  # Calculate the monthly and annual average series of temperature and preciptiation
  tmax.monthly.xts <- xts::apply.monthly(tmax.daily.xts, FUN = mean, na.rm = TRUE)
  tmin.monthly.xts <- xts::apply.monthly(tmin.daily.xts, FUN = mean, na.rm = TRUE)
  prcp.monthly.log.xts <- xts::apply.monthly(prcp.daily.log.xts, FUN = mean, na.rm = TRUE)

  # Becasue ARIMA function expects a ts object, needs to transfer the xts objects to ts objects
  # The acutal month/dates corresponding to the ts object (becasue of possible missing data in the middle of record) will be
  # assigned later
  tmax.monthly.ts <- stats::ts(tmax.monthly.xts, frequency = 12, start = c(hist.years[1], 1))
  tmin.monthly.ts <- stats::ts(tmin.monthly.xts, frequency = 12, start = c(hist.years[1], 1))
  prcp.monthly.log.ts <- stats::ts(prcp.monthly.log.xts, frequency = 12, start = c(hist.years[1], 1))

  # Annual series are used to determine the ARIMA orders for non-seasonal component
  tmax.annual.xts <- xts::apply.yearly(tmax.daily.xts, FUN = mean, na.rm = TRUE)
  tmax.annual.ts <- stats::ts(tmax.annual.xts, frequency = 1, start = hist.years[1])
  tmin.annual.xts <- xts::apply.yearly(tmin.daily.xts, FUN = mean, na.rm = TRUE)
  tmin.annual.ts <- stats::ts(tmin.annual.xts, frequency = 1, start = hist.years[1])
  prcp.annual.log.xts <- xts::apply.yearly(prcp.daily.log.xts, FUN = mean, na.rm = TRUE)
  prcp.annual.log.ts <- stats::ts(prcp.annual.log.xts, frequency = 1, start = hist.years[1])


  # auto.arima can be used to selecte the ARIMA orders for the non-seasonal components or visual inspection can be used
  # monthly precipitation is further applied with Box-Cox transformation to consider the non-normal distribution
  lambda <- forecast::BoxCox.lambda(prcp.monthly.log.ts[prcp.monthly.log.ts > 0], method = "loglik")
  prcp.monthly.log.trans.ts <- forecast::BoxCox(prcp.monthly.log.ts, lambda)

  # Determine the orders for the non-seasonal component
  # auto.arima can be used to find the orders for non-seasonal component.
  # However, mannually select orders is important, it can lead to substantial variations in the output
  # p order for the AR term of non-seasonal component is suggested to be at least one
  # For example, Chicago is selected as (1, 1, 1) orders for the non-seasonal component

  # arima.order.annual.tmax <- arimaorder(auto.arima(tmax.annual.ts, allowdrift = TRUE))
  # arima.order.annual.tmin <- arimaorder(auto.arima(tmin.annual.ts, allowdrift = TRUE))
  # arima.order.annual.prcp <- arimaorder(auto.arima(prcp.annual.log.ts, allowdrift = TRUE))

  # Now the monthly series can be fitted with SARIMA
  # Depended on the orders from the auto.arima, allow.drift term can be added as well

  arima.fit.tmax <- forecast::Arima(tmax.monthly.ts, order = c(1, 1, 1), seasonal = c(0, 1, 1))
  arima.fit.tmin <- forecast::Arima(tmin.monthly.ts, order = c(2, 1, 1), seasonal = c(0, 1, 1))
  arima.fit.prcp <- forecast::Arima(prcp.monthly.log.trans.ts, order = c(1, 1, 1), seasonal = c(0, 1, 1))

  # auto.arima can also be used to select the orders for seasonal component but it's inefficient and can be erroneous
  # arima.fit.tmax <- auto.arima(tmax.monthly.ts, seasonal = TRUE, allowdrift = TRUE)
  # arima.fit.tmin <- auto.arima(tmin.monthly.ts, seasonal = TRUE, allowdrift = TRUE)
  # arima.fit.prcp <- auto.arima(prcp.monthly.trans.ts, seasonal = TRUE, allowdrift = TRUE)

  # ARIMA forecasts can then be made for the next year
  # Note the dates shown in these arima forecast object can be off (because of missing data in the middle of records)
  arima.fcast.tmax <- forecast::forecast(arima.fit.tmax, h = 12)
  arima.fcast.tmin <- forecast::forecast(arima.fit.tmin, h = 12)
  arima.fcast.prcp <- forecast::forecast(arima.fit.prcp, h = 12)


  # Historical monthly baseline series are stored in dataframes (days for a same month have the same values)
  fitted.daily.tmax <- data.frame(Date = daily.dates,
                                  Value = stats::setNames(stats::fitted(arima.fit.tmax),
                                                   as.Date(daily.monthly.dates))[format(daily.dates, format = "%Y-%m-01")])
  fitted.daily.tmin <- data.frame(Date = daily.dates,
                                  Value = stats::setNames(stats::fitted(arima.fit.tmin),
                                                   as.Date(daily.monthly.dates))[format(daily.dates, format = "%Y-%m-01")])
  fitted.daily.prcp <- data.frame(Date = daily.dates,
                             Value = stats::setNames(BoxCox.rev(stats::fitted(arima.fit.prcp), lambda),
                                              as.Date(daily.monthly.dates))[format(daily.dates, format = "%Y-%m-01")])

  # Forecast daily dates are created here
  fcast.dates <- seq(as.Date(paste0(fcast.starting.yr,"-01-01")), as.Date(paste0(fcast.ed.yr,"-12-31")), by="days")
  fcast.index <- data.frame(Date = fcast.dates, index = format(fcast.dates, "%m%d"))

  # Create dataframes to store the daily simulations (100 sets of simulations series are created here)
  daily.tmax.simulation <- data.frame(matrix(ncol = num.sets + 1, nrow = length(fcast.dates)))
  daily.tmin.simulation <- data.frame(matrix(ncol = num.sets + 1, nrow = length(fcast.dates)))
  daily.temp.simulation <- data.frame(matrix(ncol = num.sets + 1, nrow = length(fcast.dates)))
  daily.prcp.simulation <- data.frame(matrix(ncol = num.sets + 1, nrow = length(fcast.dates)))
  colnames(daily.tmax.simulation) <- c("Dates", as.character(c(1:num.sets)))
  colnames(daily.tmin.simulation) <- c("Dates", as.character(c(1:num.sets)))
  colnames(daily.temp.simulation) <- c("Dates", as.character(c(1:num.sets)))
  colnames(daily.prcp.simulation) <- c("Dates", as.character(c(1:num.sets)))
  daily.tmax.simulation$Dates <- fcast.dates
  daily.tmin.simulation$Dates <- fcast.dates
  daily.temp.simulation$Dates <- fcast.dates
  daily.prcp.simulation$Dates <- fcast.dates

  # Block bootstrap is accomplished in the folloiwng for loop
  for (num.simulations in 1:num.sets) {

    # Forecast date storage:
    fcast.tmax.str <- NULL
    fcast.tmin.str <- NULL
    fcast.prcp.str <- NULL

    tmax.monthly.xts.new <- tmax.monthly.xts
    tmin.monthly.xts.new <- tmin.monthly.xts
    prcp.monthly.xts.new <- prcp.monthly.log.xts

    # Simulation for the first yeat is produced first
    # The date for the first year is created with:
    fcast.dates.first <- seq(as.Date(paste0(fcast.starting.yr,"-01-01")), as.Date(paste0(fcast.starting.yr,"-12-31")), by="days")
    fcast.monthly.dates.first <- seq(as.Date(paste0(fcast.starting.yr,"-01-01")), as.Date(paste0(fcast.starting.yr,"-12-01")), by="month")

    # Baseline for the first  year is added
    fcast.daily.baseline.tmax <- data.frame(Date = fcast.dates.first,
                                            Value = stats::setNames((arima.fcast.tmax$mean),
                                                                    as.Date(fcast.monthly.dates.first))[format(fcast.dates.first,
                                                                                                               format = "%Y-%m-01")])
    fcast.daily.baseline.tmin <- data.frame(Date = fcast.dates.first,
                                            Value = stats::setNames((arima.fcast.tmin$mean),
                                                                    as.Date(fcast.monthly.dates.first))[format(fcast.dates.first,
                                                                                                               format = "%Y-%m-01")])
    fcast.daily.baseline.prcp <- data.frame(Date = fcast.dates.first,
                                            Value = stats::setNames(BoxCox.rev(arima.fcast.prcp$mean, lambda),
                                                                    as.Date(fcast.monthly.dates.first))[format(fcast.dates.first,
                                                                                                               format = "%Y-%m-01")])
    for (ii in 1:num.fcast.yr) {

      current.fcast.yr <- fcast.starting.yr + ii - 1



      current.year.dates <- data.frame("fcast.date" = seq(as.Date(paste0(current.fcast.yr, "-01-01")), as.Date(paste0(current.fcast.yr, "-12-31")), by="days"))

      # Creat block residuals
      bt.hist.years <- hist.years[-c(1:30)]
      # Excluding the December residuals from the last year
      bt.hist.yr <- bt.hist.years[c(sample(c(1:length(bt.hist.years)), 11), sample(c(1:(length(bt.hist.years) - 1)), 1))]

      select.date.year <- NULL

      for (jj in 1:12) {
        select.date <- sampling_date(current.fcast.yr, bt.hist.yr, jj, extra.days.for.bt)
        select.date.eva <- data.frame("Date" = select.date$select.date)
        tmax.daily.na.eva <- merge(select.date.eva, tmax.daily.na.df)
        tmin.daily.na.eva <- merge(select.date.eva, tmin.daily.na.df)
        prcp.daily.na.eva <- merge(select.date.eva, prcp.daily.na.df)

        # Checking if there is missing values in the historical data
        # Resampling if there is
        if (any(is.na(tmax.daily.na.eva) == TRUE) | any(is.na(tmin.daily.na.eva) == TRUE) | any(is.na(prcp.daily.na.eva) == TRUE)){
          select.date <- sampling_date(current.fcast.yr, bt.hist.yr, jj, extra.days.for.bt)
          select.date.eva <- data.frame("Date" = select.date$select.date)
          tmax.daily.na.eva <- merge(select.date.eva, tmax.daily.na.df)
          tmin.daily.na.eva <- merge(select.date.eva, tmin.daily.na.df)
          prcp.daily.na.eva <- merge(select.date.eva, prcp.daily.na.df)
        }

        if (any(is.na(tmax.daily.na.eva) == TRUE) | any(is.na(tmin.daily.na.eva) == TRUE) | any(is.na(prcp.daily.na.eva) == TRUE)){
          select.date <- sampling_date(current.fcast.yr, bt.hist.yr, jj, extra.days.for.bt)
          select.date.eva <- data.frame("Date" = select.date$select.date)
        }

        select.date.year <- rbind(select.date.year, select.date)
      }

      # select.date.df <- cbind(current.year.dates, select.date.year)

      select.date.hist <- data.frame("Date" = select.date.year$select.date)
      hist.select.daily.tmax <- merge(select.date.hist, tmax.daily.all.df, by = "Date", sort = FALSE, all.x = TRUE)
      hist.select.daily.tmin <- merge(select.date.hist, tmin.daily.all.df, by = "Date", sort = FALSE, all.x = TRUE)
      hist.select.daily.prcp <- merge(select.date.hist, prcp.daily.all.df, by = "Date", sort = FALSE, all.x = TRUE)
      hist.select.daily.prcp$Value[hist.select.daily.prcp$Value == 0] <- NA

      select.date.fitted <- data.frame("Date" = select.date.year$baseline.date)
      hist.fitted.select.daily.tmax <- merge(select.date.fitted, fitted.daily.tmax, by = "Date", sort = FALSE, all.x = TRUE)
      hist.fitted.select.daily.tmin <- merge(select.date.fitted, fitted.daily.tmin, by = "Date", sort = FALSE, all.x = TRUE)
      hist.fitted.select.daily.prcp <- merge(select.date.fitted, fitted.daily.prcp, by = "Date", sort = FALSE, all.x = TRUE)

      residuals.daily.tmax <- hist.select.daily.tmax$Value - hist.fitted.select.daily.tmax$Value
      residuals.daily.tmin <- hist.select.daily.tmin$Value - hist.fitted.select.daily.tmin$Value
      residuals.daily.prcp <- hist.select.daily.prcp$Value - hist.fitted.select.daily.prcp$Value

      # Store the new simulation data:
      new.fcast.tmax <- residuals.daily.tmax + fcast.daily.baseline.tmax$Value
      fcast.tmax.str <- c(fcast.tmax.str, new.fcast.tmax)
      new.fcast.tmin <- residuals.daily.tmin + fcast.daily.baseline.tmin$Value
      fcast.tmin.str <- c(fcast.tmin.str, new.fcast.tmin)
      new.fcast.prcp <- residuals.daily.prcp + fcast.daily.baseline.prcp$Value
      new.fcast.prcp[is.na(new.fcast.prcp)] <- 0
      new.fcast.prcp[new.fcast.prcp < 0] <- 0
      fcast.prcp.str <- c(fcast.prcp.str, new.fcast.prcp)

      if (ii != num.fcast.yr) {

        year.dates <- seq(as.Date(paste0(current.fcast.yr, "-01-01")), as.Date(paste0(current.fcast.yr, "-12-31")), by="days")

        # The block residuals are added into the baseline series
        new.fcast.daily.tmax <- xts::xts(new.fcast.tmax, as.Date(year.dates, "%Y-%m-%d"))
        new.fcast.monthly.tmax <- xts::apply.monthly(new.fcast.daily.tmax, FUN = mean, na.rm = TRUE)
        tmax.monthly.xts.new <- rbind(tmax.monthly.xts.new, new.fcast.monthly.tmax)
        # Again, becasue ARIMA function expects a ts object, needs to transfer the xts objects to ts objects
        tmax.monthly.ts.new <- stats::ts(tmax.monthly.xts.new, frequency = 12, start = c(hist.years[1], 1))
        arima.fit.tmax.new <- forecast::Arima(tmax.monthly.ts.new, model = arima.fit.tmax)
        arima.fcast.tmax.new <- forecast::forecast(arima.fit.tmax.new, h = 12)

        new.fcast.daily.tmin <- xts::xts(new.fcast.tmin, as.Date(year.dates, "%Y-%m-%d"))
        new.fcast.monthly.tmin <- xts::apply.monthly(new.fcast.daily.tmin, FUN = mean, na.rm = TRUE)
        tmin.monthly.xts.new <- rbind(tmin.monthly.xts.new, new.fcast.monthly.tmin)
        tmin.monthly.ts.new <- stats::ts(tmin.monthly.xts.new, frequency = 12, start = c(hist.years[1], 1))
        arima.fit.tmin.new <- forecast::Arima(tmin.monthly.ts.new, model = arima.fit.tmin)
        arima.fcast.tmin.new <- forecast::forecast(arima.fit.tmin.new, h = 12)

        # Precipitation series needs to transfer the NA values (i.e., representing dry days) back
        new.fcast.daily.prcp <- xts::xts(new.fcast.prcp, as.Date(year.dates, "%Y-%m-%d"))
        new.fcast.monthly.prcp <- xts::apply.monthly(new.fcast.daily.prcp, FUN = mean, na.rm = TRUE)
        prcp.monthly.xts.new <- rbind(prcp.monthly.xts.new, new.fcast.monthly.prcp)
        prcp.monthly.ts.new <- stats::ts(prcp.monthly.xts.new, frequency = 12, start = c(hist.years[1], 1))
        # Reverse Box-Cox transformation is applied here as well
        prcp.monthly.ts.trans.new <- forecast::BoxCox(prcp.monthly.ts.new, lambda)
        arima.fit.prcp.new <- forecast::Arima(prcp.monthly.ts.trans.new, model = arima.fit.prcp)
        arima.fcast.prcp.new <- forecast::forecast(arima.fit.prcp.new, h = 12)

        next.year.fcast <- current.fcast.yr + 1

        # Baseline for the next forecasting year is added
        fcast.dates.new <- seq(as.Date(paste0(next.year.fcast, "-01-01")), as.Date(paste0(next.year.fcast, "-12-31")), by = "days")
        fcast.month.new <- seq(as.Date(paste0(next.year.fcast, "-01-01")), as.Date(paste0(next.year.fcast, "-12-01")), by = "month")
        fcast.daily.baseline.tmax <- data.frame(Date = fcast.dates.new,
                                                Value = stats::setNames(arima.fcast.tmax.new$mean,
                                                                        as.Date(fcast.month.new))[format(fcast.dates.new,
                                                                                                         format = "%Y-%m-01")])
        fcast.daily.baseline.tmin <- data.frame(Date = fcast.dates.new,
                                                Value = stats::setNames((arima.fcast.tmin.new$mean),
                                                                        as.Date(fcast.month.new))[format(fcast.dates.new,
                                                                                                         format = "%Y-%m-01")])
        fcast.daily.baseline.prcp <- data.frame(Date = fcast.dates.new,
                                                Value = stats::setNames(BoxCox.rev(arima.fcast.prcp.new$mean, lambda),
                                                                        as.Date(fcast.month.new))[format(fcast.dates.new,
                                                                                                         format = "%Y-%m-01")])

      }
      # Daily simulation for the next year is continued in the loop
    }

    fcast.daily.data.tmax <- fcast.tmax.str
    fcast.daily.data.tmin <- fcast.tmin.str
    fcast.daily.data.temp <- (fcast.daily.data.tmax + fcast.daily.data.tmin) / 2
    fcast.daily.log.prcp <- fcast.prcp.str
    fcast.daily.log.prcp[is.na(fcast.daily.log.prcp)] <- 0
    fcast.daily.log.prcp[fcast.daily.log.prcp < 0] <- 0
    # Reverse log-transformation is finally applied here to transform the precipitation series into original scale
    fcast.daily.data.prcp <- exp(fcast.daily.log.prcp) - 1
    fcast.daily.data.prcp[fcast.daily.data.prcp < 0.01] <- 0

    fcast.daily.xts.tmax <- xts::xts(fcast.daily.data.tmax, as.Date(fcast.dates, "%Y-%m-%d"))
    fcast.daily.xts.tmin <- xts::xts(fcast.daily.data.tmin, as.Date(fcast.dates, "%Y-%m-%d"))
    fcast.daily.xts.temp <- xts::xts(fcast.daily.data.temp, as.Date(fcast.dates, "%Y-%m-%d"))
    fcast.daily.xts.prcp <- xts::xts(fcast.daily.data.prcp, as.Date(fcast.dates, "%Y-%m-%d"))

    daily.tmax.simulation[, as.character(num.simulations)] <- round(as.numeric(fcast.daily.xts.tmax), digits = 2)
    daily.tmin.simulation[, as.character(num.simulations)] <- round(as.numeric(fcast.daily.xts.tmin), digits = 2)
    daily.temp.simulation[, as.character(num.simulations)] <- round(as.numeric(fcast.daily.xts.temp), digits = 2)
    daily.prcp.simulation[, as.character(num.simulations)] <- round(as.numeric(fcast.daily.xts.prcp), digits = 3)

  }


  fcast.output <- list("hist.obs" = hist.obs.daily,
                       "fcast.tmax" = daily.tmax.simulation, "fcast.tmin" = daily.tmin.simulation,
                       "fcast.temp" = daily.temp.simulation, "fcast.prcp" = daily.prcp.simulation)

  return(fcast.output)

}

