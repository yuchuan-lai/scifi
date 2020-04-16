#' @export
citylist <- function() {

  # Access the city list for the available compiled temperature and precipitation records
  # Records can be downloaded directly from:
  # https://kilthub.cmu.edu/projects/Use_of_historical_data_to_assess_regional_climate_change/61538
  # Descriptions of the compiled historical daily data can be found at:
  # https://doi.org/10.1175/JCLI-D-18-0630.1

  # A github repository was used to store the links to the compiled records
  # (March 2020 updated - with records for 2019)
  #' @importFrom utils read.csv
  annual.urls <- read.csv("https://github.com/yuchuan-lai/Historical-City-ClimData/raw/master/annual.url.csv")[-c(1)]
  url.city_info <- annual.urls$url[annual.urls$name == "city_info.csv"]
  city.list <- read.csv(as.character(url.city_info))[-c(1)]

  return(city.list)
}

#' @export
download <- function(city.name, record.type = "annual") {

  # Return error message if the record type is not available
  if (record.type != "daily" & record.type != "annual") {
    stop("Please re-check the record type. Available record types: `annual`, `daily` are acceptable ones")
  }

  # Downloading the historical annual data for the selected city
  if (record.type == "annual") {
    annual.urls <- read.csv("https://github.com/yuchuan-lai/Historical-City-ClimData/raw/master/annual.url.csv")[-c(1)]
    url.city_info <- annual.urls$url[annual.urls$name == "city_info.csv"]
    city.info <- read.csv(as.character(url.city_info))[-c(1)]
    city.name <- city.info$Name[city.info$Name == city.name][1]
    url.city.annual <- annual.urls$url[annual.urls$name == paste0(city.name,".csv")]
    city.record <- read.csv(as.character(url.city.annual))[-c(1)]
  }

  # Donwloading the historical daily data for the selected city
  if (record.type == "daily") {
    daily.urls <- read.csv("https://github.com/yuchuan-lai/Historical-City-ClimData/raw/master/daily.url.csv")[-c(1)]
    url.city_info <- daily.urls$url[daily.urls$name == "city_info.csv"]
    city.info <- read.csv(as.character(url.city_info))[-c(1)]
    city.ID <- city.info$ID[city.info$Name == city.name][1]
    url.city.daily <- daily.urls$url[daily.urls$name == paste0(city.ID,".csv")]
    city.record <- read.csv(as.character(url.city.daily))[-c(1)]
  }


  return(city.record)
}

#' @export
select <- function(city.hist.data, var.name, year.col.name = "Year", thresholds = 10){

  `%>%` <- dplyr::`%>%`
  # find the index based on the acceptable thresolds

  if (var.name == "Avg.Temp") {
    exclude.index <- city.hist.data$Miss.Tmax > thresholds & city.hist.data$Miss.Tmin > thresholds
  }
  if (any(var.name == c("Avg.Tmax", "Max.Tmax", "Min.Tmax", "Warm.Days", "Cold.Days"))){
    exclude.index <- city.hist.data$Miss.Tmax > thresholds
  }
  if (any(var.name == c("Avg.Tmin", "Max.Tmin", "Min.Tmin", "Warm.Nights", "Cold.Nights"))){
    exclude.index <- city.hist.data$Miss.Tmin > thresholds
  }
  if (any(var.name == c("ToT.Prcp", "Max.1.day.P", "Max.5.day.P", "Wet.Days"))){
    exclude.index <- city.hist.data$Miss.Prcp > thresholds
  }

  city.hist.data %>% dplyr::filter(., is.na(!!as.symbol(var.name)) == FALSE) %>% dplyr::select(., c(!!as.symbol(year.col.name), !!as.symbol(var.name))) -> hist.var.obs

  hist.var.obs %>% dplyr::select(., !!as.symbol(year.col.name)) -> year.col
  hist.var.obs %>% dplyr::select(., !!as.symbol(var.name)) -> var.col

  hist.var.obs <- data.frame(year.col, var.col)
  colnames(hist.var.obs) <- c("year", var.name)
  hist.var.obs[exclude.index, 2] <- NA

  return(hist.var.obs)
}
