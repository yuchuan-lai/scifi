#' @export
plot_annual <- function (fcast.output) {

  #' @importFrom ggplot2 geom_point
  #' @importFrom ggplot2 geom_line
  #' @importFrom ggplot2 geom_ribbon
  #' @importFrom ggplot2 guide_legend
  #' @importFrom ggplot2 scale_fill_manual
  #' @importFrom ggplot2 scale_color_manual
  #' @importFrom ggplot2 theme
  #' @importFrom ggplot2 element_blank
  #'
  `%>%` <- dplyr::`%>%`
  fcast.output %>% dplyr::select(., c(fcast.yrs, median, low95, upp95, low80, upp80)) %>% dplyr::filter(., is.na(median) == FALSE) -> fcast.df
  fcast.output %>% dplyr::select(., c(hist.yrs, hist.fit, low95, obs, res)) %>% dplyr::filter(., is.na(hist.fit) == FALSE) -> hist.obs.df
  fcast.output %>% dplyr::select(., c(hist.eva.yrs, hist.eva.obs)) %>% dplyr::filter(., is.na(hist.eva.obs) == FALSE) -> hist.eva.df

  if (nrow(hist.eva.df) != 0) {
    fcast.p <- ggplot2::ggplot() + ggplot2::geom_point(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = obs, col = "#FF9900")) +
      ggplot2::geom_line(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = hist.fit, col = "red"), size = 1) +
      ggplot2::geom_ribbon(data = fcast.df, ggplot2::aes(x = fcast.yrs, ymin = low95, ymax= upp95, fill = "#FF660040")) +
      ggplot2::geom_ribbon(data = fcast.df, ggplot2::aes(x = fcast.yrs, ymin = low80, ymax= upp80, fill = "#FF660060")) +
      ggplot2::geom_line(data = fcast.df, ggplot2::aes(x = fcast.yrs, y = median, col = "#660000"), size = 1) +
      ggplot2::geom_point(data = hist.eva.df, ggplot2::aes(x = hist.eva.yrs, y = hist.eva.obs, col = "#CC0033")) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "right") +
      ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "blank"), shape = c(NA, NA))),
                                 labels = c("80%", "95%"), values = c("#FF660040", "#FF660070") , name = "Prediction intervals") +
      ggplot2::scale_color_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "solid", "blank", "solid"), shape = c(16, NA, 16, NA),
                                                                                    color = c("#FF9900", "red","#CC0033", "#660000"))),
                                  labels = c("historical observations", "model fitting", "new observations", "point forecast"),
                                  values = c("#660000", "#CC0033", "#FF9900", "red"), name = " ") +
      ggplot2::xlab("year")
  } else {
    fcast.p <- ggplot2::ggplot() + ggplot2::geom_point(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = obs, col = "#FF9900")) +
      ggplot2::geom_line(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = hist.fit, col = "red"), size = 1) +
      ggplot2::geom_ribbon(data = fcast.df, ggplot2::aes(x = fcast.yrs, ymin = low95, ymax= upp95, fill = "#FF660040")) +
      ggplot2::geom_ribbon(data = fcast.df, ggplot2::aes(x = fcast.yrs, ymin = low80, ymax= upp80, fill = "#FF660060")) +
      ggplot2::geom_line(data = fcast.df, ggplot2::aes(x = fcast.yrs, y = median, col = "#660000"), size = 1) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "right") +
      ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "blank"), shape = c(NA, NA))),
                                 labels = c("80%", "95%"), values = c("#FF660040", "#FF660070") , name = "Prediction intervals") +
      ggplot2::scale_color_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "solid", "solid"), shape = c(16, NA, NA),
                                                                                    color = c("#FF9900", "red", "#660000"))),
                                  labels = c("historical observations", "model fitting", "point forecast"),
                                  values = c("#660000", "#FF9900", "red"), name = " ") +
      ggplot2::xlab("year")
  }

  return(fcast.p)
}

#' @export
plot_annual_cl <- function (fcast.output) {

  `%>%` <- dplyr::`%>%`
  fcast.output %>% dplyr::select(., c(fcast.yrs, fcast.return.level, fcast.lowerB, fcast.upperB)) %>% dplyr::filter(., is.na(fcast.return.level) == FALSE) -> fcast.df
  fcast.output %>% dplyr::select(., c(hist.yrs, hist.obs, hist.return.level, hist.lowerB, hist.upperB)) %>% dplyr::filter(., is.na(hist.obs) == FALSE) -> hist.obs.df
  fcast.output %>% dplyr::select(., c(hist.eva.yrs, hist.eva.obs)) %>% dplyr::filter(., is.na(hist.eva.obs) == FALSE) -> hist.eva.df

  if (nrow(hist.eva.df) != 0) {
    fcast.p <- ggplot2::ggplot() + ggplot2::geom_point(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = hist.obs, col = "#FF9900")) +
      ggplot2::geom_line(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = hist.return.level, col = "red"), size = 1) +
      ggplot2::geom_ribbon(data = hist.obs.df, ggplot2::aes(x = hist.yrs, ymin = hist.lowerB, ymax = hist.upperB, fill = "#FF660040")) +
      ggplot2::geom_ribbon(data = fcast.df, ggplot2::aes(x = fcast.yrs, ymin = fcast.lowerB, ymax = fcast.upperB, fill = "#FF660060")) +
      ggplot2::geom_line(data = fcast.df, ggplot2::aes(x = fcast.yrs, y = fcast.return.level, col = "#660000"), size = 1) +
      ggplot2::geom_point(data = hist.eva.df, ggplot2::aes(x = hist.eva.yrs, y = hist.eva.obs, col = "#CC0033")) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "right") +
      ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "blank"), shape = c(NA, NA))),
                                 labels = c("historical", "forecasted"), values = c("#FF660040", "#FF660070") , name = "confidence levels") +
      ggplot2::scale_color_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "solid", "blank", "solid"), shape = c(16, NA, 16, NA),
                                                                                    color = c("#FF9900", "red","#CC0033", "#660000"))),
                                  labels = c("historical observations", "historical initial return level", "new observations", "forecasted initial return level"),
                                  values = c("#660000", "#CC0033", "#FF9900", "red"), name = " ") +
      ggplot2::xlab("year")

  } else {
    fcast.p <- ggplot2::ggplot() + ggplot2::geom_point(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = hist.obs, col = "#FF9900")) +
      ggplot2::geom_line(data = hist.obs.df, ggplot2::aes(x = hist.yrs, y = hist.return.level, col = "red"), size = 1) +
      ggplot2::geom_ribbon(data = hist.obs.df, ggplot2::aes(x = hist.yrs, ymin = hist.lowerB, ymax = hist.upperB, fill = "#FF660040")) +
      ggplot2::geom_ribbon(data = fcast.df, ggplot2::aes(x = fcast.yrs, ymin = fcast.lowerB, ymax = fcast.upperB, fill = "#FF660060")) +
      ggplot2::geom_line(data = fcast.df, ggplot2::aes(x = fcast.yrs, y = fcast.return.level, col = "#660000"), size = 1) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "right") +
      ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "blank"), shape = c(NA, NA))),
                                 labels = c("historical", "forecasted"), values = c("#FF660040", "#FF660070") , name = "confidence levels") +
      ggplot2::scale_color_manual(guide = ggplot2::guide_legend(override.aes = list(linetype = c("blank", "solid", "solid"), shape = c(16, NA, NA),
                                                                                    color = c("#FF9900", "red", "#660000"))),
                                  labels = c("historical observations", "historical initial return level", "new observations", "forecasted initial return level"),
                                  values = c("#660000", "#FF9900", "red"), name = " ") +
      ggplot2::xlab("year")
  }

  return(fcast.p)
}

#' @export
plot_daily <- function (fcast.output, plot.var.name = "temp") {

  #' @importFrom  reshape2 melt
  #' @importFrom  gridExtra grid.arrange

  if ((plot.var.name %in% c("temp", "tmax", "tmin", "prcp")) == FALSE){
    stop("Please check the input var name: `temp`, `tmax`, `tmin`, `prcp` are acceptable ones")
  }

  `%>%` <- dplyr::`%>%`
  fcast.output$hist.obs -> hist.obs.daily

  if (plot.var.name == "temp") {

    # Plotting the probability density functions of daily data
    fcast.output$fcast.temp -> fcast.temp.daily
    fcast.st.yr <- lubridate::year(fcast.temp.daily$Dates[1])
    fcast.ed.yr <- lubridate::year(fcast.temp.daily$Dates[nrow(fcast.temp.daily)])

    hist.daily.tmax.xts <- xts::xts(hist.obs.daily$tmax, as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
    hist.daily.tmin.xts <- xts::xts(hist.obs.daily$tmin, as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
    hist.daily.temp.xts <- (hist.daily.tmax.xts + hist.daily.tmin.xts) / 2
    hist.daily.temp.df <- data.frame("date" = zoo::index(hist.daily.temp.xts), "variable" = rep("obs", length(hist.daily.temp.xts)), "value" = hist.daily.temp.xts)

    num.simu <- ncol(fcast.temp.daily) - 1
    colnames(fcast.temp.daily) <- c("date", paste0("simu.", seq(1:num.simu)))
    fcast.temp.daily.df <- reshape2::melt(fcast.temp.daily, id = "date")

    temp.daily.df <- rbind(fcast.temp.daily.df, hist.daily.temp.df)
    temp.daily.df %>% dplyr::filter(., is.na(value) == FALSE) -> temp.daily.df
    plot.color <- c(rep("#0000FF50", num.simu), "#333333")

    daily.PDF <- ggplot2::ggplot(data = temp.daily.df, ggplot2::aes(x = value, color = variable)) + ggplot2::geom_density() +
      ggplot2::scale_color_manual(values = plot.color) +
      ggplot2::xlab("Daily temperature (F)") +
      ggplot2::ylab("Probability density") +
      ggplot2::theme(legend.position = "none")

    # Calculate the annual average with historical observations and simulated sereis
    # Historical observations need to be filtered out with years that have substantial missing data
    hist.annual.temp <- xts::apply.yearly(hist.daily.temp.xts, FUN = mean, na.rm = TRUE)
    temp.na.count <- xts::xts(is.na(hist.daily.temp.xts), as.Date(zoo::index(hist.daily.temp.xts), format = "%Y-%m-%d"))
    temp.na.count.annual <- xts::apply.yearly(temp.na.count, FUN = sum)
    city.hist.years.all <- lubridate::year(temp.na.count.annual)[temp.na.count.annual <= 10]

    hist.annual.temp.df <- data.frame("year" = lubridate::year(zoo::index(hist.annual.temp)),
                                      "variable" = rep("obs", length(hist.annual.temp)),
                                      "value" = hist.annual.temp)
    hist.annual.temp.df %>% dplyr::filter(., year %in% city.hist.years.all) -> hist.annual.temp.df

    annual_avg_cal <- function(simu.daily.data) {
      simu.daily.xts <- xts::xts(simu.daily.data, as.Date(fcast.temp.daily$date, format = "%Y-%m-%d"))
      simu.annual.avg <- xts::apply.yearly(simu.daily.xts, FUN = mean, na.rm = TRUE)
      return(simu.annual.avg)
    }

    fcast.temp.annual <- apply(fcast.temp.daily[-c(1)], 2, annual_avg_cal)
    fcast.temp.annual <- data.frame("year" = fcast.st.yr:fcast.ed.yr, fcast.temp.annual)
    fcast.temp.annual.df <- reshape2::melt(fcast.temp.annual, id = "year")


    annual.trend.p <- ggplot2::ggplot() + ggplot2::geom_line(data = fcast.temp.annual.df, ggplot2::aes(x = year, y = value, color = variable)) +
      ggplot2::scale_color_manual(values = rep("#0000FF50", num.simu)) +
      ggplot2::geom_point(data = hist.annual.temp.df, ggplot2::aes(x = year, y = value), color = "#333333") +
      ggplot2::xlab("Year") +
      ggplot2::ylab("Annual average daily mean temperature (F)") +
      ggplot2::theme(legend.position = "none")

  }

  if (plot.var.name == "tmax") {

    fcast.output$fcast.tmax -> fcast.tmax.daily
    fcast.st.yr <- lubridate::year(fcast.tmax.daily$Dates[1])
    fcast.ed.yr <- lubridate::year(fcast.tmax.daily$Dates[nrow(fcast.tmax.daily)])

    hist.daily.tmax.xts <- xts::xts(hist.obs.daily$tmax, as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
    hist.daily.tmax.df <- data.frame("date" = zoo::index(hist.daily.tmax.xts), "variable" = rep("obs", length(hist.daily.tmax.xts)), "value" = hist.daily.tmax.xts)

    num.simu <- ncol(fcast.tmax.daily) - 1
    colnames(fcast.tmax.daily) <- c("date", paste0("simu.", seq(1:num.simu)))
    fcast.tmax.daily.df <- reshape2::melt(fcast.tmax.daily, id = "date")

    tmax.daily.df <- rbind(fcast.tmax.daily.df, hist.daily.tmax.df)
    tmax.daily.df %>% dplyr::filter(., is.na(value) == FALSE) -> tmax.daily.df
    plot.color <- c(rep("#0000FF50", num.simu), "#333333")

    daily.PDF <- ggplot2::ggplot(data = tmax.daily.df, ggplot2::aes(x = value, color = variable)) + ggplot2::geom_density() +
      ggplot2::scale_color_manual(values = plot.color) +
      ggplot2::xlab("Daily maximum temperature (F)") +
      ggplot2::ylab("Probability density") +
      ggplot2::theme(legend.position = "none")

    # Calculate the annual average with historical observations and simulated sereis
    # Historical observations need to be filtered out with years that have substantial missing data
    hist.annual.tmax <- xts::apply.yearly(hist.daily.tmax.xts, FUN = mean, na.rm = TRUE)
    tmax.na.count <- xts::xts(is.na(hist.daily.tmax.xts), as.Date(zoo::index(hist.daily.tmax.xts), format = "%Y-%m-%d"))
    tmax.na.count.annual <- xts::apply.yearly(tmax.na.count, FUN = sum)
    city.hist.years.all <- lubridate::year(tmax.na.count.annual)[tmax.na.count.annual <= 10]

    hist.annual.tmax.df <- data.frame("year" = lubridate::year(zoo::index(hist.annual.tmax)),
                                      "variable" = rep("obs", length(hist.annual.tmax)),
                                      "value" = hist.annual.tmax)
    hist.annual.tmax.df %>% dplyr::filter(., year %in% city.hist.years.all) -> hist.annual.tmax.df

    annual_avg_cal <- function(simu.daily.data) {
      simu.daily.xts <- xts::xts(simu.daily.data, as.Date(fcast.tmax.daily$date, format = "%Y-%m-%d"))
      simu.annual.avg <- xts::apply.yearly(simu.daily.xts, FUN = mean, na.rm = TRUE)
      return(simu.annual.avg)
    }

    fcast.tmax.annual <- apply(fcast.tmax.daily[-c(1)], 2, annual_avg_cal)
    fcast.tmax.annual <- data.frame("year" = fcast.st.yr:fcast.ed.yr, fcast.tmax.annual)
    fcast.tmax.annual.df <- reshape2::melt(fcast.tmax.annual, id = "year")

    annual.trend.p <- ggplot2::ggplot() + ggplot2::geom_line(data = fcast.tmax.annual.df, ggplot2::aes(x = year, y = value, color = variable)) +
      ggplot2::scale_color_manual(values = rep("#0000FF50", num.simu)) +
      ggplot2::geom_point(data = hist.annual.tmax.df, ggplot2::aes(x = year, y = value), color = "#333333") +
      ggplot2::xlab("Year") +
      ggplot2::ylab("Annual average daily maximum temperature (F)") +
      ggplot2::theme(legend.position = "none")

  }


  if (plot.var.name == "tmin") {

    fcast.output$fcast.tmin -> fcast.tmin.daily
    fcast.st.yr <- lubridate::year(fcast.tmin.daily$Dates[1])
    fcast.ed.yr <- lubridate::year(fcast.tmin.daily$Dates[nrow(fcast.tmin.daily)])

    hist.daily.tmin.xts <- xts::xts(hist.obs.daily$tmin, as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
    hist.daily.tmin.df <- data.frame("date" = zoo::index(hist.daily.tmin.xts), "variable" = rep("obs", length(hist.daily.tmin.xts)), "value" = hist.daily.tmin.xts)

    num.simu <- ncol(fcast.tmin.daily) - 1
    colnames(fcast.tmin.daily) <- c("date", paste0("simu.", seq(1:num.simu)))
    fcast.tmin.daily.df <- reshape2::melt(fcast.tmin.daily, id = "date")

    tmin.daily.df <- rbind(fcast.tmin.daily.df, hist.daily.tmin.df)
    tmin.daily.df %>% dplyr::filter(., is.na(value) == FALSE) -> tmin.daily.df
    plot.color <- c(rep("#0000FF50", num.simu), "#333333")

    daily.PDF <- ggplot2::ggplot(data = tmin.daily.df, ggplot2::aes(x = value, color = variable)) + ggplot2::geom_density() +
      ggplot2::scale_color_manual(values = plot.color) +
      ggplot2::xlab("Daily minimum temperature (F)") +
      ggplot2::ylab("Probability density") +
      ggplot2::theme(legend.position = "none")

    # Calculate the annual average with historical observations and simulated sereis
    # Historical observations need to be filtered out with years that have substantial missing data
    hist.annual.tmin <- xts::apply.yearly(hist.daily.tmin.xts, FUN = mean, na.rm = TRUE)
    tmin.na.count <- xts::xts(is.na(hist.daily.tmin.xts), as.Date(zoo::index(hist.daily.tmin.xts), format = "%Y-%m-%d"))
    tmin.na.count.annual <- xts::apply.yearly(tmin.na.count, FUN = sum)
    city.hist.years.all <- lubridate::year(tmin.na.count.annual)[tmin.na.count.annual <= 10]

    hist.annual.tmin.df <- data.frame("year" = lubridate::year(zoo::index(hist.annual.tmin)),
                                      "variable" = rep("obs", length(hist.annual.tmin)),
                                      "value" = hist.annual.tmin)
    hist.annual.tmin.df %>% dplyr::filter(., year %in% city.hist.years.all) -> hist.annual.tmin.df

    annual_avg_cal <- function(simu.daily.data) {
      simu.daily.xts <- xts::xts(simu.daily.data, as.Date(fcast.tmin.daily$date, format = "%Y-%m-%d"))
      simu.annual.avg <- xts::apply.yearly(simu.daily.xts, FUN = mean, na.rm = TRUE)
      return(simu.annual.avg)
    }

    fcast.tmin.annual <- apply(fcast.tmin.daily[-c(1)], 2, annual_avg_cal)
    fcast.tmin.annual <- data.frame("year" = fcast.st.yr:fcast.ed.yr, fcast.tmin.annual)
    fcast.tmin.annual.df <- reshape2::melt(fcast.tmin.annual, id = "year")

    annual.trend.p <- ggplot2::ggplot() + ggplot2::geom_line(data = fcast.tmin.annual.df, ggplot2::aes(x = year, y = value, color = variable)) +
      ggplot2::scale_color_manual(values = rep("#0000FF50", num.simu)) +
      ggplot2::geom_point(data = hist.annual.tmin.df, ggplot2::aes(x = year, y = value), color = "#333333") +
      ggplot2::xlab("Year") +
      ggplot2::ylab("Annual average daily minimum temperature (F)") +
      ggplot2::theme(legend.position = "none")
  }


  if (plot.var.name == "prcp") {

    fcast.output$fcast.prcp -> fcast.prcp.daily
    fcast.st.yr <- lubridate::year(fcast.prcp.daily$Dates[1])
    fcast.ed.yr <- lubridate::year(fcast.prcp.daily$Dates[nrow(fcast.prcp.daily)])

    hist.daily.prcp.xts <- xts::xts(hist.obs.daily$prcp, as.Date(hist.obs.daily$Date, format = "%Y-%m-%d"))
    hist.daily.prcp.df <- data.frame("date" = zoo::index(hist.daily.prcp.xts), "variable" = rep("obs", length(hist.daily.prcp.xts)), "value" = hist.daily.prcp.xts)

    num.simu <- ncol(fcast.prcp.daily) - 1
    colnames(fcast.prcp.daily) <- c("date", paste0("simu.", seq(1:num.simu)))
    fcast.prcp.daily.df <- reshape2::melt(fcast.prcp.daily, id = "date")

    prcp.daily.df <- rbind(fcast.prcp.daily.df, hist.daily.prcp.df)
    prcp.daily.df %>% dplyr::filter(., is.na(value) == FALSE & value >= 1/ 25.4) -> prcp.daily.df
    plot.color <- c(rep("#0000FF50", num.simu), "#333333")

    daily.PDF <- ggplot2::ggplot(data = prcp.daily.df, ggplot2::aes(x = value, color = variable)) + ggplot2::geom_density() +
      ggplot2::scale_color_manual(values = plot.color) +
      ggplot2::xlab("Daily precipitation from wet days (in.)") +
      ggplot2::ylab("Probability density") +
      ggplot2::coord_cartesian(xlim = c(0, stats::quantile(prcp.daily.df$value, 0.95, na.rm = TRUE))) +
      ggplot2::theme(legend.position = "none")

    # Calculate the annual average with historical observations and simulated sereis
    # Historical observations need to be filtered out with years that have substantial missing data
    hist.annual.prcp <- xts::apply.yearly(hist.daily.prcp.xts, FUN = sum, na.rm = TRUE)
    prcp.na.count <- xts::xts(is.na(hist.daily.prcp.xts), as.Date(zoo::index(hist.daily.prcp.xts), format = "%Y-%m-%d"))
    prcp.na.count.annual <- xts::apply.yearly(prcp.na.count, FUN = sum)
    city.hist.years.all <- lubridate::year(prcp.na.count.annual)[prcp.na.count.annual <= 10]

    hist.annual.prcp.df <- data.frame("year" = lubridate::year(zoo::index(hist.annual.prcp)),
                                      "variable" = rep("obs", length(hist.annual.prcp)),
                                      "value" = hist.annual.prcp)
    hist.annual.prcp.df %>% dplyr::filter(., year %in% city.hist.years.all) -> hist.annual.prcp.df

    annual_tot_cal <- function(simu.daily.data) {
      simu.daily.xts <- xts::xts(simu.daily.data, as.Date(fcast.prcp.daily$date, format = "%Y-%m-%d"))
      simu.annual.avg <- xts::apply.yearly(simu.daily.xts, FUN = sum, na.rm = TRUE)
      return(simu.annual.avg)
    }

    fcast.prcp.annual <- apply(fcast.prcp.daily[-c(1)], 2, annual_tot_cal)
    fcast.prcp.annual <- data.frame("year" = fcast.st.yr:fcast.ed.yr, fcast.prcp.annual)
    fcast.prcp.annual.df <- reshape2::melt(fcast.prcp.annual, id = "year")

    annual.trend.p <- ggplot2::ggplot() + ggplot2::geom_line(data = fcast.prcp.annual.df, ggplot2::aes(x = year, y = value, color = variable)) +
      ggplot2::scale_color_manual(values = rep("#0000FF50", num.simu)) +
      ggplot2::geom_point(data = hist.annual.prcp.df, ggplot2::aes(x = year, y = value), color = "#333333") +
      ggplot2::xlab("Year") +
      ggplot2::ylab("Annual total precipitation (in.)") +
      ggplot2::theme(legend.position = "none")
  }

  fcast.p <- gridExtra::grid.arrange(daily.PDF, annual.trend.p, widths = c(1, 1))

  return(fcast.p)
}
