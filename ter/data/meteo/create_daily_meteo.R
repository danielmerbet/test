library(dplyr)

dir <- "~/Documents/NEREIDA/"
predam <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_data_predam.csv"))
postdam <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_data_postdam.csv"))
coast <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_data_coast.csv"))
outlet <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_data_outlet.csv"))

coast$date <- as.Date(coast$date)
predam$date <- as.Date(predam$date)
postdam$date <- as.Date(postdam$date)
outlet$date <- as.Date(outlet$date)

###tmin and tmax
coast_daily <- coast %>%
  group_by(date) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
tmin <- coast %>%
  group_by(date) %>%
  summarise(tmin = min(temperature_2m, na.rm = TRUE))
coast_daily$tmin <- tmin$tmin
tmax <- coast %>%
  group_by(date) %>%
  summarise(tmax = max(temperature_2m, na.rm = TRUE))
coast_daily$tmax <- tmax$tmax

predam_daily <- predam %>%
  group_by(date) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
tmin <- predam %>%
  group_by(date) %>%
  summarise(tmin = min(temperature_2m, na.rm = TRUE))
predam_daily$tmin <- tmin$tmin
tmax <- predam %>%
  group_by(date) %>%
  summarise(tmax = max(temperature_2m, na.rm = TRUE))
predam_daily$tmax <- tmax$tmax

postdam_daily <- postdam %>%
  group_by(date) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
tmin <- postdam %>%
  group_by(date) %>%
  summarise(tmin = min(temperature_2m, na.rm = TRUE))
postdam_daily$tmin <- tmin$tmin
tmax <- postdam %>%
  group_by(date) %>%
  summarise(tmax = max(temperature_2m, na.rm = TRUE))
postdam_daily$tmax <- tmax$tmax

outlet_daily <- outlet %>%
  group_by(date) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
tmin <- outlet %>%
  group_by(date) %>%
  summarise(tmin = min(temperature_2m, na.rm = TRUE))
outlet_daily$tmin <- tmin$tmin
tmax <- outlet %>%
  group_by(date) %>%
  summarise(tmax = max(temperature_2m, na.rm = TRUE))
outlet_daily$tmax <- tmax$tmax

#precipitation is acucmulated
pr <- coast %>%
  group_by(date) %>%
  summarise(pr = sum(precipitation, na.rm = TRUE))
coast_daily$pr <- pr$pr
coast_daily$precipitation <- NULL

pr <- predam %>%
  group_by(date) %>%
  summarise(pr = sum(precipitation, na.rm = TRUE))
predam_daily$pr <- pr$pr
predam_daily$precipitation <- NULL

pr <- postdam %>%
  group_by(date) %>%
  summarise(pr = sum(precipitation, na.rm = TRUE))
postdam_daily$pr <- pr$pr
postdam_daily$precipitation <- NULL

pr <- outlet %>%
  group_by(date) %>%
  summarise(pr = sum(precipitation, na.rm = TRUE))
outlet_daily$pr <- pr$pr
outlet_daily$precipitation <- NULL

###potential evaporation
dates <- seq(from = coast$date[1], by = 1, len = nrow(coast_daily))
compute_Ra <- function(date, phi) {
  doy <- yday(date)  # Day of year
  
  # Constants
  Gsc <- 0.0820  # MJ/m²/min (solar constant)
  
  # Compute dr (inverse relative distance Earth-Sun)
  dr <- 1 + 0.033 * cos((2 * pi / 365) * doy)
  
  # Compute solar declination (delta)
  delta <- 0.409 * sin((2 * pi / 365) * doy - 1.39)
  
  # Compute sunset hour angle (omega_s)
  omega_s <- acos(-tan(phi) * tan(delta))
  
  # Compute Ra (MJ/m²/day)
  Ra <- (24 * 60 / pi) * Gsc * dr * 
    (omega_s * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(omega_s))
  
  return(Ra)
}
hargreaves_samani <- function(tmax, tmin, Ra) {
  #tmean <- (tmax + tmin) / 2
  #PET <- 0.0023 * (tmean + 17.8) * sqrt(tmax - tmin) * Ra
  trng <- tmax - tmin
  trng[trng < 0] <- 0  # Avoid negative values
  PET <- 0.0023 * (Ra / 2.45) * ((tmax + tmin) / 2 + 17.8) * sqrt(trng)
  return(PET)
}

#coast
target_lat <- 41.97
phi <- target_lat * pi / 180  # Convert to radians
Ra_values <- sapply(dates, function(date) compute_Ra(date, phi))
PET <- hargreaves_samani(coast_daily$tmax, coast_daily$tmin, Ra_values)
coast_daily$pet <- PET
range(coast_daily$pet)

#predam
target_lat <- 42.3
phi <- target_lat * pi / 180  # Convert to radians
Ra_values <- sapply(dates, function(date) compute_Ra(date, phi))
PET <- hargreaves_samani(predam_daily$tmax, predam_daily$tmin, Ra_values)
predam_daily$pet <- PET
range(predam_daily$pet)

#postdam
target_lat <- 42.06
phi <- target_lat * pi / 180  # Convert to radians
Ra_values <- sapply(dates, function(date) compute_Ra(date, phi))
PET <- hargreaves_samani(postdam_daily$tmax, postdam_daily$tmin, Ra_values)
postdam_daily$pet <- PET
range(postdam_daily$pet)

#outlet
target_lat <- 42.06
phi <- target_lat * pi / 180  # Convert to radians
Ra_values <- sapply(dates, function(date) compute_Ra(date, phi))
PET <- hargreaves_samani(outlet_daily$tmax, outlet_daily$tmin, Ra_values)
outlet_daily$pet <- PET
range(outlet_daily$pet)

#save data
write.csv(coast_daily, file=paste0(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_coast.csv")),
          quote = F, row.names = F)
write.csv(predam_daily, file=paste0(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_predam.csv")),
          quote = F, row.names = F)
write.csv(postdam_daily, file=paste0(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_postdam.csv")),
          quote = F, row.names = F)
write.csv(outlet_daily, file=paste0(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_outlet.csv")),
          quote = F, row.names = F)
