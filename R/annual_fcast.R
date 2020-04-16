#' @export
fcast_annual <- function(hist.obs, num.fcast.yrs = 20, fcast.starting.yr = 2020, forecast.model = "ARIMA") {

  #' @import forecast
  #' @import dplyr
  #' @import Kendall
  #' @import extRemes

  #' @importFrom stats fitted
  #' @importFrom stats lm
  #' @importFrom stats median
  #' @importFrom stats na.omit
  #' @importFrom stats predict.lm
  #' @importFrom stats qnorm
  #' @importFrom stats quantile
  #' @importFrom stats residuals
  #' @importFrom stats shapiro.test


  if (forecast.model == "ARIMA") {
    fcast.output <- ARIMA_fcast_annual(hist.obs, num.fcast.yrs, fcast.starting.yr)
  }
  if (forecast.model == "linear") {
    fcast.output <- LR_fcast_annual(hist.obs, num.fcast.yrs, fcast.starting.yr)
  }
  if (forecast.model == "baselineall") {
    fcast.output <- Baselineall_fcast_annual(hist.obs, num.fcast.yrs, fcast.starting.yr)
  }
  if (forecast.model == "baseline30") {
    fcast.output <- Bl30_fcast_annual(hist.obs, num.fcast.yrs, fcast.starting.yr)
  }
  if (forecast.model == "gev") {
    fcast.output <- GEV_fcast_annual(hist.obs, num.fcast.yrs, fcast.starting.yr)
  }
  return(fcast.output)
}

BoxCox.rev <- function(x, lambda)  {
  if (lambda == 0) {
    exp(x)
  } else {
    (lambda*x +1)^(1/lambda)
  }
}

ARIMA_fcast_annual <- function(hist.obs, num.fcast.yrs, fcast.starting.yr) {


  # as some cities may have missing years in the recent record, the following steps are used to allow a fixable starting year
  # for forecasting.
  var.name <- colnames(hist.obs)[2]
  colnames(hist.obs) <- c("year", "obs")

  `%>%` <- dplyr::`%>%`
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., is.na(obs) == FALSE & year <= fcast.starting.yr - 1) ->  hist.fit.obs

  hist.yrs <- hist.fit.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., year >= fcast.st.yr) ->  hist.eva.obs
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1
  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1
  num.eva.yrs <- nrow(hist.eva.obs)
  hist.eva.data <- hist.eva.obs[, 2]
  hist.obs.data <- hist.fit.obs[, 2]
  fcast.yrs <- c(fcast.st.yr:fcast.ed.yr)
  fit.yrs <- hist.yrs

  trend.test <- Kendall::MannKendall(hist.obs.data)

  if (trend.test$sl < 0.2) {
    apply.drift <- TRUE
  } else {
    apply.drift <- FALSE
  }

  # use ARIMA to forecast time series
  arima.test.fit <- forecast::auto.arima(hist.obs.data, allowdrift = apply.drift)

  # the use of any(arimaorder(arima.test.fit) != 0) is considering if one variable has same numbers for all of the years
  # for example, a city in a warm region may have no days with daily minimum temperature below 32F
  if (any(forecast::arimaorder(arima.test.fit) != 0)) {
    normal.test <- shapiro.test(arima.test.fit$residuals)
    if (normal.test$p.value < 0.05) {
      apply.BoxCox <- TRUE
    } else {
      apply.BoxCox <- FALSE
    }
  } else{
    normal.test <- shapiro.test(hist.obs.data)
    if (normal.test$p.value < 0.05) {
      apply.BoxCox <- TRUE
    } else {
      apply.BoxCox <- FALSE
    }
  }


  if (apply.BoxCox == TRUE) {

    adjust.value <- 0

    if (min(hist.obs.data, na.rm = TRUE) == 0) {
      adjust.value <- min(hist.obs.data[hist.obs.data > 0], na.rm = TRUE) * 1/2
    }

    if (min(hist.obs.data, na.rm = TRUE) < 0) {
      if (max(hist.obs.data, na.rm = TRUE) < 0) {
        adjust.value <- abs(min(hist.obs.data, na.rm = TRUE) + max(hist.obs.data, na.rm = TRUE)* 1/2)
      } else {
        adjust.value <- min(hist.obs.data[hist.obs.data > 0], na.rm = TRUE) * 1/2 + abs(min(hist.obs.data, na.rm = TRUE))
      }}

    hist.obs.data <- hist.obs.data + adjust.value

    # Re-assess the trend with transformed data
    lambda <- forecast::BoxCox.lambda(hist.obs.data, method = "loglik")
    hist.obs.trans <- forecast::BoxCox(hist.obs.data, lambda)
    trend.test.2 <- Kendall::MannKendall(hist.obs.trans)
    if (trend.test.2$sl < 0.2) {
      apply.drift <- TRUE
    } else {
      apply.drift <- FALSE
    }

    auto.arima.fit <- forecast::auto.arima(hist.obs.data, allowdrift = apply.drift, lambda = lambda)
    arima.order <- forecast::arimaorder(auto.arima.fit)
    arima.fit <- forecast::Arima(hist.obs.data, order = arima.order, include.drift = apply.drift, lambda = auto.arima.fit$lambda)

    arima.fcast <- forecast::forecast(arima.fit, h = num.fcast.yrs.actual, level=c(80, 95))
    pred95.upp.arima <- arima.fcast$upper[,2] - adjust.value
    # Prediction level at 95% can return unrealistically large values when lambda = -1
    if (auto.arima.fit$lambda == -1 & any(pred95.upp.arima > max(hist.obs.data, na.rm = TRUE) * 1.5) & any(is.na(pred95.upp.arima) == FALSE)) {
      pred95.upp.arima <- rep(NA, length(arima.fcast$upper[,2]))
    }
    pred95.low.arima <- arima.fcast$lower[,2] - adjust.value
    pred80.upp.arima <- arima.fcast$upper[,1] - adjust.value
    pred80.low.arima <- arima.fcast$lower[,1] - adjust.value
    pred.mean <- arima.fcast$mean - adjust.value
    hist.fit <- fitted(arima.fit) - adjust.value
    residuals.fit <- hist.obs.data - hist.fit

    if (any(var.name == c("Avg.Temp", "Avg.Tmax", "Avg.Tmin", "TXx", "TNx", "TXn", "TNn")) == FALSE)  {
      pred95.upp.arima[as.numeric(pred95.upp.arima) < 0] <- 0
      pred95.low.arima[as.numeric(pred95.low.arima) < 0] <- 0
      pred80.upp.arima[as.numeric(pred80.upp.arima) < 0] <- 0
      pred80.low.arima[as.numeric(pred80.low.arima) < 0] <- 0
      pred.mean[as.numeric(pred.mean) < 0] <- 0
    }


  } else {

    arima.order <- forecast::arimaorder(forecast::auto.arima(hist.obs.data, allowdrift = apply.drift))
    # Again, here needs to consider if the variable has the same number for all of the years
    # In other words, the three orders for the ARIMA model will all be 0
    if (any(arima.order != 0)){
      arima.fit <- forecast::Arima(hist.obs.data, order = arima.order, include.drift = apply.drift)
    } else {
      arima.fit <- forecast::auto.arima(hist.obs.data)
    }
    arima.fcast <- forecast::forecast(arima.fit, h = num.fcast.yrs.actual, level=c(80, 95))
    pred95.upp.arima <- arima.fcast$upper[,2]
    pred95.low.arima <- arima.fcast$lower[,2]
    pred80.upp.arima <- arima.fcast$upper[,1]
    pred80.low.arima <- arima.fcast$lower[,1]
    pred.mean <- arima.fcast$mean
    hist.fit <- fitted(arima.fit)
    residuals.fit <- residuals(arima.fit)

    if (any(var.name == c("Avg.Temp", "Avg.Tmax", "Avg.Tmin", "TXx", "TNx", "TXn", "TNn")) == FALSE)  {
      pred95.upp.arima[as.numeric(pred95.upp.arima) < 0] <- 0
      pred95.low.arima[as.numeric(pred95.low.arima) < 0] <- 0
      pred80.upp.arima[as.numeric(pred80.upp.arima) < 0] <- 0
      pred80.low.arima[as.numeric(pred80.low.arima) < 0] <- 0
      pred.mean[as.numeric(pred.mean) < 0] <- 0
    }

  }

  fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred.mean <-c(pred.mean, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred95.low <- c(pred95.low.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred95.upp <- c(pred95.upp.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred80.low <- c(pred80.low.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred80.upp <- c(pred80.upp.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))

  hist.eva.yrs.df <- c(hist.eva.obs[, 1], rep(NA, (length(fit.yrs) - num.eva.yrs)))
  hist.eva.obs.df <- c(hist.eva.data, rep(NA, (length(fit.yrs) - num.eva.yrs)))

  fcast.output <- data.frame("fcast.yrs" = fcast.yrs.df, "median" = pred.mean, "low95" = pred95.low,
                             "upp95" = pred95.upp, "low80" = pred80.low, "upp80" = pred80.upp,
                             "hist.yrs" = fit.yrs, "hist.fit" = hist.fit, "obs" = hist.obs.data, "res" = residuals.fit,
                             "hist.eva.yrs" = hist.eva.yrs.df, "hist.eva.obs" = hist.eva.obs.df)
  return(fcast.output)
}

LR_fcast_annual <- function(hist.obs, num.fcast.yrs, fcast.starting.yr) {

  colnames(hist.obs) <- c("year", "obs")

  `%>%` <- dplyr::`%>%`
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., is.na(obs) == FALSE & year <= fcast.starting.yr - 1) ->  hist.fit.obs
  hist.yrs <- hist.fit.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., year >= fcast.st.yr) ->  hist.eva.obs
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1

  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1

  num.eva.yrs <- nrow(hist.eva.obs)
  hist.eva.data <- hist.eva.obs[, 2]
  hist.obs.data <- hist.fit.obs[, 2]

  fcast.yrs <- c(fcast.st.yr:fcast.ed.yr)
  fit.yrs <- hist.yrs


  hist.obs.df <- data.frame("year" = fit.yrs, "obs" = hist.obs.data)
  LR.fit <- lm(obs ~ year, hist.obs.df)
  hist.fit <- LR.fit$fitted.values
  residuals.fit <- LR.fit$residuals
  fcast.yrs.df <- data.frame("year" = fcast.yrs)
  LR.fcast_80 <- predict.lm(LR.fit, fcast.yrs.df, interval = "prediction", level = 0.8)[c(1:num.fcast.yrs.actual),]
  LR.fcast_95 <- predict.lm(LR.fit, fcast.yrs.df, interval = "prediction", level = 0.95)[c(1:num.fcast.yrs.actual),]
  fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred.mean <-c(LR.fcast_80[, 1], rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred95.low <- c(LR.fcast_95[, 2], rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred95.upp <- c(LR.fcast_95[, 3], rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred80.low <- c(LR.fcast_80[, 2], rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred80.upp <- c(LR.fcast_80[, 3], rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))

  hist.eva.yrs.df <- c(hist.eva.obs[, 1], rep(NA, (length(fit.yrs) - num.eva.yrs)))
  hist.eva.obs.df <- c(hist.eva.data, rep(NA, (length(fit.yrs) - num.eva.yrs)))

  fcast.output <- data.frame("fcast.yrs" = fcast.yrs.df, "median" = pred.mean, "low95" = pred95.low,
                             "upp95" = pred95.upp, "low80" = pred80.low, "upp80" = pred80.upp,
                             "hist.yrs" = fit.yrs, "hist.fit" = hist.fit, "obs" = hist.obs.data, "res" = residuals.fit,
                             "hist.eva.yrs" = hist.eva.yrs.df, "hist.eva.obs" = hist.eva.obs.df)
  return(fcast.output)
}

Baselineall_fcast_annual <- function(hist.obs, num.fcast.yrs, fcast.starting.yr) {

  colnames(hist.obs) <- c("year", "obs")

  `%>%` <- dplyr::`%>%`
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., is.na(obs) == FALSE & year <= fcast.starting.yr - 1) ->  hist.fit.obs
  hist.yrs <- hist.fit.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., year >= fcast.st.yr) ->  hist.eva.obs
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1

  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1

  num.eva.yrs <- nrow(hist.eva.obs)
  hist.eva.data <- hist.eva.obs[, 2]
  hist.obs.data <- hist.fit.obs[, 2]

  fcast.yrs <- c(fcast.st.yr:fcast.ed.yr)
  fit.yrs <- hist.yrs

  fcast.data <- quantile(hist.obs.data, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), na.rm = TRUE)

  hist.fit <- rep(fcast.data[3], length(hist.yrs))
  residuals.fit <- hist.obs.data - hist.fit
  fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred.mean <-c(rep(fcast.data[3], num.fcast.yrs), rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred95.low <- c(rep(fcast.data[1], num.fcast.yrs), rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred95.upp <- c(rep(fcast.data[5], num.fcast.yrs), rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred80.low <- c(rep(fcast.data[2], num.fcast.yrs), rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred80.upp <- c(rep(fcast.data[4], num.fcast.yrs), rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))

  hist.eva.yrs.df <- c(hist.eva.obs[, 1], rep(NA, (length(fit.yrs) - num.eva.yrs)))
  hist.eva.obs.df <- c(hist.eva.data, rep(NA, (length(fit.yrs) - num.eva.yrs)))

  fcast.output <- data.frame("fcast.yrs" = fcast.yrs.df, "median" = pred.mean, "low95" = pred95.low,
                             "upp95" = pred95.upp, "low80" = pred80.low, "upp80" = pred80.upp,
                             "hist.yrs" = fit.yrs, "hist.fit" = hist.fit, "obs" = hist.obs.data, "res" = residuals.fit,
                             "hist.eva.yrs" = hist.eva.yrs.df, "hist.eva.obs" = hist.eva.obs.df)
  return(fcast.output)
}

Bl30_fcast_annual <- function(hist.obs, num.fcast.yrs, fcast.starting.yr) {

  colnames(hist.obs) <- c("year", "obs")

  `%>%` <- dplyr::`%>%`
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., is.na(obs) == FALSE & year <= fcast.starting.yr - 1) ->  hist.fit.obs
  hist.yrs <- hist.fit.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., year >= fcast.st.yr) ->  hist.eva.obs
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1

  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1

  num.eva.yrs <- nrow(hist.eva.obs)
  hist.eva.data <- hist.eva.obs[, 2]
  hist.obs.data <- hist.fit.obs[, 2]

  fcast.yrs <- c(fcast.st.yr:fcast.ed.yr)
  fit.yrs <- hist.yrs

  hist.obs.df <- data.frame("year" = hist.yrs, "obs" = hist.obs.data)
  row.st.num <- nrow(hist.obs.df) - 29
  hist.obs.df %>% dplyr::slice(., row.st.num:dplyr::n()) -> hist.level.df
  fit.yrs <- hist.level.df$year


  fcast.data <- quantile(hist.level.df$obs, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), na.rm = TRUE)

  hist.fit <- rep(fcast.data[3], 30)
  residuals.fit <- hist.level.df$obs - hist.fit
  fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(hist.yrs) - num.fcast.yrs.actual)))
  pred.mean <-c(rep(fcast.data[3], num.fcast.yrs), rep(NA, (length(hist.yrs) - num.fcast.yrs.actual)))
  pred95.low <- c(rep(fcast.data[1], num.fcast.yrs), rep(NA, (length(hist.yrs) - num.fcast.yrs.actual)))
  pred95.upp <- c(rep(fcast.data[5], num.fcast.yrs), rep(NA, (length(hist.yrs) - num.fcast.yrs.actual)))
  pred80.low <- c(rep(fcast.data[2], num.fcast.yrs), rep(NA, (length(hist.yrs) - num.fcast.yrs.actual)))
  pred80.upp <- c(rep(fcast.data[4], num.fcast.yrs), rep(NA, (length(hist.yrs) - num.fcast.yrs.actual)))


  hist.fit <- c(rep(NA, (nrow(hist.obs.df) - 30)), hist.fit)
  residuals.fit <- c(rep(NA, (nrow(hist.obs.df) - 30)), residuals.fit)

  hist.eva.yrs.df <- c(hist.eva.obs[, 1], rep(NA, (length(hist.yrs) - num.eva.yrs)))
  hist.eva.obs.df <- c(hist.eva.data, rep(NA, (length(hist.yrs) - num.eva.yrs)))

  fcast.output <- data.frame("fcast.yrs" = fcast.yrs.df, "median" = pred.mean, "low95" = pred95.low,
                             "upp95" = pred95.upp, "low80" = pred80.low, "upp80" = pred80.upp,
                             "hist.yrs" = hist.yrs, "hist.fit" = hist.fit, "obs" = hist.obs.data, "res" = residuals.fit,
                             "hist.eva.yrs" = hist.eva.yrs.df, "hist.eva.obs" = hist.eva.obs.df)
  return(fcast.output)
}

GEV_fcast_annual <- function(hist.obs, num.fcast.yrs, fcast.starting.yr) {

  colnames(hist.obs) <- c("year", "obs")

  `%>%` <- dplyr::`%>%`
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., is.na(obs) == FALSE & year <= fcast.starting.yr - 1) ->  hist.fit.obs
  hist.yrs <- hist.fit.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  hist.obs %>% dplyr::select(., year, obs) %>% dplyr::filter(., year >= fcast.st.yr) ->  hist.eva.obs
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1

  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1

  num.eva.yrs <- nrow(hist.eva.obs)
  hist.eva.data <- hist.eva.obs[, 2]
  hist.obs.data <- hist.fit.obs[, 2]

  fcast.yrs <- c(fcast.st.yr:fcast.ed.yr)
  fit.yrs <- hist.yrs

  trend.test <- Kendall::MannKendall(hist.obs.data)
  if (trend.test$sl < 0.2) {
    apply.trend <- TRUE
  } else {
    apply.trend <- FALSE
  }

  years.df <- data.frame("years" = fit.yrs)

  gev.test <- tryCatch({
    gev.fit <- extRemes::fevd(x = as.numeric(hist.obs.data), data = years.df, na.action = na.omit,
                    location.fun = ~ years + 1, type = "GEV", method = "MLE")
    apply.nonsGEV <- TRUE
  }, error = function(error) {
    apply.nonsGEV <- FALSE
    return(apply.nonsGEV)
  })

  if (apply.trend == TRUE & gev.test == TRUE) {
    gev.fit <- extRemes::fevd(x = as.numeric(hist.obs.data), data = years.df, na.action = na.omit,
                    location.fun = ~ years + 1, type = "GEV", method = "MLE")
    apply.nonsGEV <- TRUE
  } else {
    apply.nonsGEV <- FALSE
    gev.fit <- extRemes::fevd(x = as.numeric(hist.obs.data), type = "GEV", method = "Lmoments")
    gev.mean.fit <- rep(extRemes::return.level(gev.fit, return.period = 2), length(fit.yrs))
    fcast.mean.gev <- rep(extRemes::return.level(gev.fit, return.period = 2), num.fcast.yrs.actual)
    pred95.upp.gev <- rep(extRemes::return.level(gev.fit, return.period = 40), num.fcast.yrs.actual)
    pred95.low.gev <- rep(extRemes::return.level(gev.fit, return.period = 100/97.5), num.fcast.yrs.actual)
    pred80.upp.gev <- rep(extRemes::return.level(gev.fit, return.period = 10), num.fcast.yrs.actual)
    pred80.low.gev <- rep(extRemes::return.level(gev.fit, return.period = 10/9), num.fcast.yrs.actual)

    hist.fit <- gev.mean.fit
    residuals.fit <- hist.obs.data - gev.mean.fit
    fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
    pred.mean <-c(fcast.mean.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
    pred95.low <- c(pred95.low.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
    pred95.upp <- c(pred95.upp.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
    pred80.low <- c(pred80.low.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
    pred80.upp <- c(pred80.upp.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  }

  # GEV fitting with MLE can lead to unrealistic fitting, therefore it's crucial to check if the results of fitting
  # is realistic first.
  if (apply.nonsGEV == TRUE) {
      if (extRemes::erlevd(gev.fit, period = 2)[1] <= max(hist.obs.data, na.rm = TRUE)) {

      num.yrs <- length(extRemes::erlevd(gev.fit, period = 2))
      diff.yrs <- as.numeric(names(extRemes::erlevd(gev.fit, period = 2))[num.yrs]) - as.numeric(names(extRemes::erlevd(gev.fit, period = 2))[num.yrs - 1])

      gev.mean.fit <- extRemes::erlevd(gev.fit, period = 2)
      gev.mean.fcast <- seq(from = extRemes::erlevd(gev.fit, period = 2)[num.yrs] + (extRemes::erlevd(gev.fit, period = 2)[num.yrs] - extRemes::erlevd(gev.fit, period = 2)[num.yrs - 1]) / diff.yrs,
                      by = (extRemes::erlevd(gev.fit, period = 2)[num.yrs] -  extRemes::erlevd(gev.fit, period = 2)[num.yrs - 1]) / diff.yrs,
                      length.out = num.fcast.yrs.actual)

      pred95.upp.gev <- seq(from = extRemes::erlevd(gev.fit, period = 40)[num.yrs] + (extRemes::erlevd(gev.fit, period = 40)[num.yrs] - extRemes::erlevd(gev.fit, period = 40)[num.yrs - 1]) / diff.yrs,
                               by = (extRemes::erlevd(gev.fit, period = 40)[num.yrs] -  extRemes::erlevd(gev.fit, period = 40)[num.yrs - 1]) / diff.yrs,
                               length.out = num.fcast.yrs.actual)
      pred95.low.gev <- seq(from = extRemes::erlevd(gev.fit, period = 100/97.5)[num.yrs] + (extRemes::erlevd(gev.fit, period = 100/97.5)[num.yrs] - extRemes::erlevd(gev.fit, period = 100/97.5)[num.yrs - 1]) / diff.yrs,
                               by = (extRemes::erlevd(gev.fit, period = 100/97.5)[num.yrs] -  extRemes::erlevd(gev.fit, period = 100/97.5)[num.yrs - 1]) / diff.yrs,
                               length.out = num.fcast.yrs.actual)
      pred80.upp.gev <- seq(from = extRemes::erlevd(gev.fit, period = 10)[num.yrs] + (extRemes::erlevd(gev.fit, period = 10)[num.yrs] - extRemes::erlevd(gev.fit, period = 10)[num.yrs - 1]) / diff.yrs,
                               by = (extRemes::erlevd(gev.fit, period = 10)[num.yrs] -  extRemes::erlevd(gev.fit, period = 10)[num.yrs - 1]) / diff.yrs,
                               length.out = num.fcast.yrs.actual)
      pred80.low.gev <- seq(from = extRemes::erlevd(gev.fit, period = 10/9)[num.yrs] + (extRemes::erlevd(gev.fit, period = 10/9)[num.yrs] - extRemes::erlevd(gev.fit, period = 10/9)[num.yrs - 1]) / diff.yrs,
                               by = (extRemes::erlevd(gev.fit, period = 10/9)[num.yrs] -  extRemes::erlevd(gev.fit, period = 10/9)[num.yrs - 1]) / diff.yrs,
                               length.out = num.fcast.yrs.actual)
      fcast.mean.gev <- seq(from = extRemes::erlevd(gev.fit, period = 2)[num.yrs] + (extRemes::erlevd(gev.fit, period = 2)[num.yrs] -  extRemes::erlevd(gev.fit, period = 2)[num.yrs - 1]) / diff.yrs,
                               by = (extRemes::erlevd(gev.fit, period = 2)[num.yrs] -  extRemes::erlevd(gev.fit, period = 2)[num.yrs - 1]) / diff.yrs,
                               length.out = num.fcast.yrs.actual)
      hist.fit <- gev.mean.fit
      residuals.fit <- hist.obs.data - gev.mean.fit
      fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
      pred.mean <-c(fcast.mean.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
      pred95.low <- c(pred95.low.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
      pred95.upp <- c(pred95.upp.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
      pred80.low <- c(pred80.low.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
      pred80.upp <- c(pred80.upp.gev, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))

    } else {
      hist.fit <- rep(NA, length(fit.yrs))
      residuals.fit <- rep(NA, length(fit.yrs))
      fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
      pred.mean <- rep(NA, length(fit.yrs))
      pred95.low <- rep(NA, length(fit.yrs))
      pred95.upp <- rep(NA, length(fit.yrs))
      pred80.low <- rep(NA, length(fit.yrs))
      pred80.upp <- rep(NA, length(fit.yrs))
    }
  }

  hist.eva.yrs.df <- c(hist.eva.obs[, 1], rep(NA, (length(fit.yrs) - num.eva.yrs)))
  hist.eva.obs.df <- c(hist.eva.data, rep(NA, (length(fit.yrs) - num.eva.yrs)))

  fcast.output <- data.frame("fcast.yrs" = fcast.yrs.df, "median" = pred.mean, "low95" = pred95.low,
                             "upp95" = pred95.upp, "low80" = pred80.low, "upp80" = pred80.upp,
                             "hist.yrs" = fit.yrs, "hist.fit" = hist.fit, "obs" = hist.obs.data, "res" = residuals.fit,
                             "hist.eva.yrs" = hist.eva.yrs.df, "hist.eva.obs" = hist.eva.obs.df)
  return(fcast.output)
}
