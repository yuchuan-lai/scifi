# Functions used for ARIMA_conf_level
ARIMA_init_fcast_func <- function(hist.obs, num.fcast.yrs, fcast.starting.yr, return.period) {

  # as some cities may have missing years in the recent record, the following steps are used to allow a fixable starting year
  # for forecasting.

  colnames(hist.obs) <- c("year", "obs")

  hist.yrs <- hist.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1
  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1
  hist.obs.data <- hist.obs[, 2]
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

    pred.level <- 100 - 2 * 100 / return.period

    arima.fcast <- forecast::forecast(arima.fit, h = num.fcast.yrs.actual, level = pred.level)
    pred.upp.arima <- arima.fcast$upper - adjust.value
    pred.low.arima <- arima.fcast$lower - adjust.value
    pred.mean <- arima.fcast$mean - adjust.value
    hist.fit <- fitted(arima.fit) - adjust.value
    residuals.fit <- hist.obs.data - hist.fit

    arima.fit.mean.trans <- forecast::BoxCox(fitted(arima.fit), auto.arima.fit$lambda)
    arima.fit.std <- rep(sqrt(arima.fit$sigma2), length(hist.yrs))

    percentile.level <- 1 - 1 / return.period
    hist.return.level <- BoxCox.rev(qnorm(percentile.level, arima.fit.mean.trans, arima.fit.std), auto.arima.fit$lambda) - adjust.value

  } else {

    arima.order <- forecast::arimaorder(forecast::auto.arima(hist.obs.data, allowdrift = apply.drift))
    # Again, here needs to consider if the variable has the same number for all of the years
    # In other words, the three orders for the ARIMA model will all be 0
    if (any(arima.order != 0)){
      arima.fit <- forecast::Arima(hist.obs.data, order = arima.order, include.drift = apply.drift)
    } else {
      arima.fit <- forecast::auto.arima(hist.obs.data)
    }

    pred.level <- 100 - 2 * 100 / return.period

    arima.fcast <- forecast::forecast(arima.fit, h = num.fcast.yrs.actual, level = pred.level)
    pred.upp.arima <- arima.fcast$upper
    pred.low.arima <- arima.fcast$lower
    pred.mean <- arima.fcast$mean
    hist.fit <- fitted(arima.fit)
    residuals.fit <- residuals(arima.fit)

    arima.fit.mean <- fitted(arima.fit)
    arima.fit.std <- rep(sqrt(arima.fit$sigma2), length(hist.yrs))

    percentile.level <- 1 - 1 / return.period
    hist.return.level <- qnorm(percentile.level, arima.fit.mean, arima.fit.std)

  }

  fcast.yrs.df <- c(fcast.yrs, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred.mean <-c(pred.mean, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred.low <- c(pred.low.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  pred.upp <- c(pred.upp.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))
  return.level <- c(pred.upp.arima, rep(NA, (length(fit.yrs) - num.fcast.yrs.actual)))

  fcast.output <- data.frame("fcast.yrs" = fcast.yrs.df, "median" = pred.mean, "pred.low" = pred.low,
                             "pred.upp" = pred.upp, "fcast.return.level" = return.level,
                             "hist.yrs" = fit.yrs, "hist.fit" = hist.fit, "hist.obs" = hist.obs.data, "res" = residuals.fit, "hist.return.level" = hist.return.level)
  return(fcast.output)
}

ARIMA_bootstrap_fcast_func <- function(hist.obs.data, ARIMA.init.order, num.fcast.yrs, ARIMA.lambda, return.period) {

  if (ARIMA.lambda != FALSE) {

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

    hist.obs.data <- unlist(hist.obs.data + adjust.value, use.names = FALSE)

    arima.fit <- forecast::Arima(hist.obs.data, order = ARIMA.init.order, lambda = ARIMA.lambda)

    pred.level <- 100 - 2 * 100 / return.period

    arima.fcast <- forecast::forecast(arima.fit, h = num.fcast.yrs, level = pred.level)
    pred.upp.arima <- arima.fcast$upper - adjust.value

    arima.fit.mean.trans <- forecast::BoxCox(fitted(arima.fit), ARIMA.lambda)
    arima.fit.std <- rep(sqrt(arima.fit$sigma2), length(hist.obs.data))

    percentile.level <- 1 - 1 / return.period
    hist.return.level <- BoxCox.rev(qnorm(percentile.level, arima.fit.mean.trans, arima.fit.std), ARIMA.lambda)  - adjust.value

  }
  if (ARIMA.lambda == FALSE)  {

    hist.obs.data <- unlist(hist.obs.data, use.names = FALSE)

    if (any(ARIMA.init.order != 0)){
      arima.fit <- forecast::Arima(hist.obs.data, order = ARIMA.init.order)
    } else {
      arima.fit <- forecast::auto.arima(hist.obs.data)
    }

    pred.level <- 100 - 2 * 100 / return.period

    arima.fcast <- forecast::forecast(arima.fit, h = num.fcast.yrs, level = pred.level)
    pred.upp.arima <- arima.fcast$upper

    arima.fit.mean <- fitted(arima.fit)
    arima.fit.std <- rep(sqrt(arima.fit$sigma2), length(hist.obs.data))

    percentile.level <- 1 - 1 / return.period

    hist.return.level <- qnorm(percentile.level, arima.fit.mean, arima.fit.std)
  }

  return.level <- c(hist.return.level, pred.upp.arima)

  return(return.level)
}

ARIMA_init_order <- function(hist.obs, num.fcast.yrs, fcast.starting.yr) {

  hist.yrs <- hist.obs$year
  fcast.st.yr <- hist.yrs[length(hist.yrs)] + 1
  fcast.ed.yr <- fcast.starting.yr + num.fcast.yrs - 1
  num.fcast.yrs.actual <- fcast.ed.yr - fcast.st.yr + 1
  hist.obs.data <- hist.obs[, 2]
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
    arima.fit.stats <- data.frame("p" = arima.order[1], "d" = arima.order[2], "q" = arima.order[3],
                                  "apply.BoxCox" = apply.BoxCox, "lambda" = auto.arima.fit$lambda)

  } else {

    arima.order <- forecast::arimaorder(forecast::auto.arima(hist.obs.data, allowdrift = apply.drift))
    # Again, here needs to consider if the variable has the same number for all of the years
    # In other words, the three orders for the ARIMA model will all be 0
    if (any(arima.order != 0)){
      arima.fit.stats <- data.frame("p" = arima.order[1], "d" = arima.order[2], "q" = arima.order[3],
                                    "apply.BoxCox" = apply.BoxCox)
          } else {
            arima.fit.stats <- data.frame("p" = 0, "d" = 0, "q" = 0, "apply.drift" = apply.drift,
                                          "apply.BoxCox" = apply.BoxCox)
      }


  }

  return(arima.fit.stats)
}

# The main function for the estimation of confidence intervals for ARIMA forecasting of annual extremes
#' @export
fcast_annual_cl <- function(hist.obs, num.fcast.yrs = 20, fcast.starting.yr = 2020, return.period = 5, num.bt = 100, confidence.level = 0.9) {

  #' @importFrom matrixStats rowQuantiles

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

  ARIMA.init.fcast.output <- ARIMA_init_fcast_func(hist.fit.obs, num.fcast.yrs, fcast.st.yr, return.period)

  ARIMA.init.stats <- ARIMA_init_order(hist.fit.obs, num.fcast.yrs, fcast.st.yr)

  apply.BoxCox <- ARIMA.init.stats$apply.BoxCox
  if (apply.BoxCox) {
    ARIMA.init.order <- c(ARIMA.init.stats$p, ARIMA.init.stats$d, ARIMA.init.stats$q)
    ARIMA.lambda <- ARIMA.init.stats$lambda
  } else {
    ARIMA.init.order <- c(ARIMA.init.stats$p, ARIMA.init.stats$d, ARIMA.init.stats$q)
    ARIMA.lambda <- FALSE
  }

  ARIMA.init.fcast.output %>% dplyr::select(., hist.fit) -> init.fitting
  bt.include.num <- nrow(init.fitting) - ARIMA.init.stats$p - ARIMA.init.stats$d
  ARIMA.init.fcast.output %>% dplyr::select(., res) %>% utils::tail(., bt.include.num) -> init.residuals
  num.hist.yrs <- nrow(init.fitting)

  bt.return.level.str <- NULL

  for (ii in 1:num.bt) {
    bt.residuals <- sample(t(init.residuals), size = num.hist.yrs, replace = TRUE)

    bt.hist.obs <- init.fitting + bt.residuals
    bt.return.level <- ARIMA_bootstrap_fcast_func(bt.hist.obs, ARIMA.init.order, num.fcast.yrs, ARIMA.lambda, return.period)
    bt.return.level.str <- cbind(bt.return.level.str, bt.return.level)

  }

  conf.level.low <- (1 - confidence.level) / 2
  conf.level.upp <- (1 - confidence.level) / 2 + confidence.level

  conf.level.return.level <- matrixStats::rowQuantiles(bt.return.level.str, probs = c(conf.level.low, conf.level.upp))

  ARIMA.init.fcast.output %>% dplyr::select(., fcast.yrs, fcast.return.level, hist.yrs, hist.obs, hist.return.level) -> return.level.df

  year.all <- c(hist.fit.obs$year, fcast.st.yr:(fcast.st.yr + num.fcast.yrs - 1))
  bt.return.level <- data.frame("year" = year.all, "lowerB" = conf.level.return.level[, 1], "upperB" = conf.level.return.level[, 2])

  bt.return.level %>% dplyr::filter(., year <= fcast.st.yr - 1) %>% dplyr::select(., - year) -> hist.return.level.conf
  bt.return.level %>% dplyr::filter(., year >= fcast.st.yr) %>% dplyr::select(., - year)  -> fcast.return.level.conf

  colnames(hist.return.level.conf) <- c("hist.lowerB", "hist.upperB")
  colnames(fcast.return.level.conf) <- c("fcast.lowerB", "fcast.upperB")

  num.fit.yrs <- nrow(hist.fit.obs)
  fcast.return.level.fill <- data.frame(matrix(rep(NA, (num.fit.yrs - num.fcast.yrs.actual) * 2), ncol = 2))
  colnames(fcast.return.level.fill) <- c("fcast.lowerB", "fcast.upperB")

  hist.eva.yrs.df <- c(hist.eva.obs[, 1], rep(NA, (num.fit.yrs - num.eva.yrs)))
  hist.eva.obs.df <- c(hist.eva.data, rep(NA, (num.fit.yrs - num.eva.yrs)))

  output.df <- data.frame(return.level.df, hist.return.level.conf, return.level.df, rbind(fcast.return.level.conf, fcast.return.level.fill),
                     "hist.eva.yrs" = hist.eva.yrs.df, "hist.eva.obs" = hist.eva.obs.df)

  col.order <- c("fcast.yrs", "fcast.return.level", "fcast.lowerB", "fcast.upperB",
                 "hist.yrs", "hist.obs", "hist.return.level", "hist.lowerB", "hist.upperB",
                 "hist.eva.yrs", "hist.eva.obs")
  output.df <- output.df[, col.order]


  return(output.df)
}

