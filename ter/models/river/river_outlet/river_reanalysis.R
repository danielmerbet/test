library(airGR);
library(imputeTS)
library(lubridate); library(hydroGOF)
#para tener hydroGOF: 
#paso 1: install.packages("raster"):
#paso2:install.packages("hydroGOF")
#library(hydroGOF)

dir <- "~/Documents/NEREIDA/"
area <- 1713500000 #all basin: 3171660000.000 #m², baix: 1713500000

#load climate
meteo <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_down.csv"))
meteo$date <- as.Date(meteo$date)
initial_reanalysis_date <- meteo$date[1] #initial date after calibration, 2 years for spin-up, since seasonal forecast stars in 1993, there is no problem

#load observed discharge (it must be changed every new forecast)
Q <- read.csv(paste0(dir,"ter/data/water/outlet/q_aca_outlet_ter.csv"))
colnames(Q) <- c("date", "Q")
Q$date <- as.Date(Q$date, format="%Y/%m/%d")
Q$Qls <- Q$Q*1000
Q$date <- as.Date(Q$date)

#Merge hydrologic and meteorological data
meteo_Q <- merge(meteo, Q, by="date")

#data ready to input in hydrologic model
data_airGR <- data.frame(DatesR=meteo_Q$date, 
                         P=meteo_Q$pr, 
                         T=meteo_Q$temperature_2m, 
                         E=meteo_Q$pet, 
                         Qls=meteo_Q$Qls, 
                         Qmm=meteo_Q$Qls*86400/area) #m3/s   m/aream2    mm   86400s/d
#l 1m3        =  m  86400s 1000mm = mm      
#s 1000l Am2  =  s    d      1m   = d
data_airGR$Qls[data_airGR$Qls<0] <- NA
data_airGR$Qmm[data_airGR$Qmm<0] <- NA
#data_airGR$Qls <- na_ma(data_airGR$Qls)
#data_airGR$Qmm <- na_ma(data_airGR$Qmm)

summary(data_airGR)

#Creating inputs for model
data_airGR <- data_airGR[!duplicated(data_airGR$DatesR), ]

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = as.POSIXct(data_airGR$DatesR),
                                 Precip = data_airGR$P, PotEvap = data_airGR$E)
#,Qupstream = Qupstream, LengthHydro = LengthHydro,BasinAreas = BasinAreas

#str(InputsModel)

#Select time period to run and warm-up 
WarmUp_ini <- as.Date(InputsModel$DatesR[1])
WarmUp_end <- WarmUp_ini+500
Run_ini <- WarmUp_end+1
Run_end <- data_airGR$DatesR[length(data_airGR$DatesR)]
IndPeriod_WarmUp <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_ini), 
                        which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_end))
Ind_Run <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_ini), 
               which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_end))

#Run options
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IndPeriod_WarmUp = IndPeriod_WarmUp)
#str(RunOptions)

#Other statistical criteria for simulation (NSE)
#InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
#                               RunOptions = RunOptions, VarObs = "Q", Obs = data_airGR$Qmm[Ind_Run])
#str(InputsCrit)

#selected statistical criteria for simulation (KGE)
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, VarObs = "Q", 
                               Obs = data_airGR$Qmm[Ind_Run])
#str(InputsCrit)

#Calibration options (the ranges were taken according to literature)
cal_ranges <- as.matrix(data.frame(X1=c(50,2000), X2=c(-7.5,6), X3=c(10,450), X4=c(0.25,4.35)))
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, 
                                   FUN_CALIB = Calibration_Michel,
                                   SearchRanges = cal_ranges)
#str(CalibOptions)

#calculate best parameter fitting discharge observation
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, 
                                   RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, 
                                   CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J,
                                   verbose = TRUE)
Param <- OutputsCalib$ParamFinalR

#Simulation run
OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)

pdf(paste0(dir, "ter/river_down/Ter_GR4J.pdf"))
plot(OutputsModel, Qobs = data_airGR$Qmm[Ind_Run])
dev.off()

#write.csv(data.frame(obs=data_airGR$Qmm[Ind_Run], sim=OutputsModel$Qsim), 
#          paste0("/home/ry4902/Documents/Workflow/Complete/Ter_GR4J_dat.txt"),
#          quote = F, row.names = F)

#pdf(paste0(dir_out, "GR4J/KGE_cal.pdf"))
pdf(paste0(dir, "ter/river_down/KGE_cal.pdf"))
ggof(OutputsModel$Qsim, data_airGR$Qmm[Ind_Run], dates=data_airGR$DatesR[Ind_Run], 
     ftype = "dm", FUN = mean)
dev.off()

#save parmeters calibration for future simualtion (forecast)
save(Param, file=paste0(dir, "ter/river_down/parameter_calibration.RData"))

###add 
#upstream
## the reservoir withdraws 1 m3/s when it's possible considering the flow observed in the basin
Qupstream <- matrix(-sapply(data_airGR$Qls / 1000 - 1, function(x) {
  min(1, max(0, x, na.rm = TRUE))
}), ncol = 1)

## except between July and September when the reservoir releases 3 m3/s for low-flow mitigation
month <- as.numeric(format(data_airGR$DatesR, "%m"))
Qupstream[month >= 7 & month <= 9] <- 3
Qupstream <- Qupstream * 86400 ## Conversion in m3/day

## the reservoir is not an upstream subcachment: its areas is NA
BasinAreas <- c(NA, area)

## delay time between the reservoir and the catchment outlet is 2 days and the distance is 150 km
LengthHydro <- 10

## with a delay of 2 days for 150 km, the flow velocity is 75 km per day
Velocity <- (LengthHydro * 1e3 / 2) / (24 * 60 * 60) ## Conversion km/day -> m/s

#Creating inputs for model with reanalysis time range
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = as.POSIXct(meteo_Q$date),
                                 Precip = meteo_Q$pr, PotEvap = meteo_Q$pet)

WarmUp_ini <- as.Date(InputsModel$DatesR[1])
WarmUp_end <- WarmUp_ini+500
Run_ini <- WarmUp_end+1
Run_end <- meteo$date[length(meteo$date)]
IndPeriod_WarmUp <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_ini), 
                        which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_end))
Ind_Run <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_ini), 
               which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_end))
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IndPeriod_WarmUp = IndPeriod_WarmUp)
#Simulation run with reanalysis time range
OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, 
                              RunOptions = RunOptions, 
                              Param = Param)

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = as.POSIXct(meteo_Q$date),
                                 Precip = meteo_Q$pr, PotEvap = meteo_Q$pet,
                                 Qupstream = Qupstream, 
                                 LengthHydro = LengthHydro,
                                 BasinAreas = BasinAreas)

OutputsModel_lag <- RunModel_Lag(InputsModel = InputsModel,
                                 RunOptions = RunOptions,
                                 Param = Velocity,
                                 QcontribDown = OutputsModel)

#######
pdf(paste0(dir, "ter/river_down/Ter_GR4J_lag.pdf"))
plot(OutputsModel_lag, Qobs = data_airGR$Qmm[Ind_Run])
dev.off()

#write.csv(data.frame(obs=data_airGR$Qmm[Ind_Run], sim=OutputsModel$Qsim), 
#          paste0("/home/ry4902/Documents/Workflow/Complete/Ter_GR4J_dat.txt"),
#          quote = F, row.names = F)

#pdf(paste0(dir_out, "GR4J/KGE_cal.pdf"))
pdf(paste0(dir, "ter/river_down/KGE_cal_lag.pdf"))
ggof(OutputsModel$Qsim, data_airGR$Qmm[Ind_Run], dates=data_airGR$DatesR[Ind_Run], 
     ftype = "dm", FUN = mean)
dev.off()

#########3
pos_ini <- which(as.Date(OutputsModel$DatesR)==initial_reanalysis_date)
pos_end <- length(OutputsModel$DatesR)
discharge_data <- OutputsModel$Qsim[pos_ini:pos_end]*area/(1000*86400)
attr(discharge_data, "dimensions") <- c("time")

#set discharge results in netcdf
discharge_reanalysis <- reanalysis_int$tp
discharge_reanalysis <- subsetGrid(discharge_reanalysis, years = year(initial_reanalysis_date):year(OutputsModel$DatesR[pos_end]))
discharge_reanalysis$Data <- discharge_data

save(discharge_reanalysis, file=paste(dir_out ,"discharge_reanalysis.RData", sep=""))
