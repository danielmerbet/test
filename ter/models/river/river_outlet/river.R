#This code has two main steps: the (1) river forecasting modeling and (2) probabilistic plots:
#Created by D. Mercado-Bettín

library(lubridate); library(sp); library(airGR);

dir <- "~/Documents/NEREIDA/"
area <- 1713500000 #all basin: 3171660000.000 #m², baix: 1713500000

spinup <- 2

#load data
reanalysis_data <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_down.csv"))

forecast_data <- list(tp=forecast_list$tp, pev=forecast_list$pev)
reanalysis_data <- list(tp=reanalysis_int$tp, pev=reanalysis_int$pev)

forecast_dates <- seq(from=as.Date("2024-10-01"), by=1, length.out=215)

pos_forecast <- forecast_dates %in% dates_total_forecast
pos_reanalysis <- as.Date(reanalysis_data$tp$Dates$start) %in% dates_total

discharge_data <- matrix(, 51, length(dates_total_forecast))
for (member in 1:51){
  #discharge_out <- c()
  #for (y in 1:length(years)){
    #if (months[1]==1){ #if forecast starts in january we only need data from previous years
    #  reanalysis_data_var <- subsetGrid(reanalysis_data, years = (years_rean[y]-spinup):(year-1))
    #}else{ #otherwise we need some month from the same year of the forecast
    #  reanalysis_data_var1 <- lapply(reanalysis_data, function(x) subsetGrid(x, years = (years_rean[y]-spinup):(years_rean[y]-1)))
    #  reanalysis_data_var2 <- lapply(reanalysis_data, function(x) subsetGrid(x, years = years_rean[y], season = 1:(months[1]-1)))       
    #  reanalysis_data_var <- lapply(names(reanalysis_data), function(x) bindGrid(reanalysis_data_var1[[x]], reanalysis_data_var2[[x]], dimension = c("time")))
    #  names(reanalysis_data_var) <- names(reanalysis_data)
    #  reanalysis_data_var1 <- NULL
    #  reanalysis_data_var2 <- NULL
    #}
    #forecast_data_var <- lapply(forecast_data, function(x) subsetGrid(x, years = years_fore[y]+1, members = member))
    #forecast_data_var <- lapply(forecast_data_var, function(x){
    #  x$Dates$start <- paste(as.Date(x$Dates$start), "00:00:00 GMT")
    #  x$Dates$end <- paste(as.Date(x$Dates$start), "00:00:00 GMT")
    #  return(x)
    #})
    #climate_var <- lapply(names(forecast_data), function(x) bindGrid(reanalysis_data_var[[x]], forecast_data_var[[x]], dimension = c("time"))) 
    #names(climate_var) <- names(forecast_data_var)
    
    #climate_var <- lapply(names(forecast_data), function(x) bindGrid(reanalysis_data_var[[x]], forecast_data_var[[x]], dimension = c("time"))) 
    #names(climate_var) <- names(forecast_data_var)
    
    meteo <- data.frame(dates=dates_total,
                        P=c(reanalysis_data$tp$Data[pos_reanalysis], forecast_list$tp[pos_forecast,member]),
                        E=c(reanalysis_data$pev$Data[pos_reanalysis], forecast_list$pev[pos_forecast,member]))
    meteo$P[meteo$P<0] <- 0
    meteo$E[meteo$E<0] <- 0
    
    #run hydrologic model
    InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = as.POSIXct(meteo$dates),
                                     Precip = meteo$P, PotEvap = meteo$E)
    WarmUp_ini <- as.Date(InputsModel$DatesR[1])
    WarmUp_end <- date_ini-1
    Run_ini <- WarmUp_end+1
    Run_end <- meteo$dates[length(meteo$dates)]
    IndPeriod_WarmUp <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_ini), 
                            which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_end))
    Ind_Run <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_ini), 
                   which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_end))
    RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                   InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                   IndPeriod_WarmUp = IndPeriod_WarmUp)
    #Simulation run with reanalysis time range
    OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
    #pos_ini <- which(as.Date(OutputsModel$DatesR)==initial_reanalysis_date)
    #pos_end <- length(OutputsModel$DatesR)
    daily_discharge <- OutputsModel$Qsim*area/(1000*86400)
    
    #discharge_out <- c(discharge_out, daily_discharge)
    
  #}
  discharge_data[member,] <- daily_discharge
  print(paste("finished member:",member))
}

attr(discharge_data, "dimensions") <- c("member", "time")

#save data just in case, before converting to netcdf
save(discharge_data, file=paste0(getwd(), "/discharge_forecast.RData"))
