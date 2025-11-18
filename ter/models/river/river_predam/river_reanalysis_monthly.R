library(airGR);
library(imputeTS)
library(lubridate); library(hydroGOF)
#para tener hydroGOF: 
#paso 1: install.packages("raster"):
#paso2:install.packages("hydroGOF")
#library(hydroGOF)

dir <- "~/Documents/NEREIDA/"
area <- 1802660000.000 #all basin: 3171660000.000 #m², alt: 1802660000.000

#load climate
meteo <- read.csv(paste0(dir, "ter/data/meteo/meteo-soil_dailydata_up.csv"))
meteo$date <- as.Date(meteo$date)
initial_reanalysis_date <- meteo$date[1] #initial date after calibration, 2 years for spin-up, since seasonal forecast stars in 1993, there is no problem

#load observed discharge (it must be changed every new forecast)
#Q <- read.csv(paste0(dir,"ter/data/water/predam/q_aca_predam_ter.csv"))
Q <- read.csv(paste0(dir,"ter/data/water/predam/q_atl_predam_ter.csv"))
colnames(Q) <- c("date", "Q")
Q$date <- as.Date(Q$date, format="%m/%d/%Y")
#Q$date <- as.Date(Q$date, format="%Y/%m/%d")
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

datairgr_obs <- data_airGR[complete.cases(data_airGR),]

## conversion of observed discharge to monthly time step
datairgr_obs$month <- format(datairgr_obs$DatesR, "%Y-%m")

Qobs <- datairgr_obs[c("month","Qmm")]

Qobs <- Qobs %>%
  group_by(month) %>%
  summarise(
    across(c(Qmm), ~sum(.x, na.rm = TRUE))
  )

#Qobs <- SeriesAggreg(data.frame(as.POSIXct(data_airGR$DatesR), data_airGR$Qmm),
#                     Format = "%Y-%m",
#                     ConvertFun = "sum")
#

#create daily inputs
data_airGR_month <- datairgr_obs[c("month","P", "E")]
data_airGR_month <- data_airGR_month %>%
  group_by(month) %>%
  summarise(
    across(c(P,E), ~sum(.x, na.rm = TRUE))
  )

InputsModel <- list(DatesR = as.POSIXct(paste0(data_airGR_month$month, "-01")),
     Precip = data_airGR_month$P, PotEvap=data_airGR_month$E)

class(InputsModel) <- c("InputsModel", "monthly", "GR")

#InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M(), 
#                                 DatesR = as.POSIXct(paste0(data_airGR_month$month, "-01")),
#                                 Precip = data_airGR$P, PotEvap = data_airGR$E)
                                 
#InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = as.POSIXct(data_airGR$DatesR),
#                                 Precip = data_airGR$P, PotEvap = data_airGR$E)

## conversion of InputsModel to monthly time step
#InputsModel <- SeriesAggreg(InputsModel, Format = "%Y-%m")

#Select time period to run and warm-up 
WarmUp_ini <- as.Date(InputsModel$DatesR[1])+1
WarmUp_end <- WarmUp_ini+500
Run_ini <- as.Date(format(WarmUp_end + 31, "%Y-%m-01")) #WarmUp_end+1
Run_end <- data_airGR$DatesR[length(data_airGR$DatesR)]

## run period selection
Ind_Run <- seq(which(format(InputsModel$DatesR, format = "%Y-%m")==format(Run_ini, format="%Y-%m")),
               which(format(InputsModel$DatesR, format = "%Y-%m")==format(Run_end, format="%Y-%m")))
IndPeriod_WarmUp <- seq(which(format(InputsModel$DatesR, format = "%Y-%m") == format(WarmUp_ini, format="%Y-%m")), 
                        which(format(InputsModel$DatesR, format = "%Y-%m") == format(WarmUp_end, format="%Y-%m")))

## preparation of the RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IndPeriod_WarmUp = IndPeriod_WarmUp)

## efficiency criterion: Nash-Sutcliffe Efficiency
Qobs <- Qobs[Ind_Run, 2]
InputsCrit  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE, InputsModel = InputsModel,
                                RunOptions = RunOptions, Obs = Qobs)

cal_ranges <- as.matrix(data.frame(X1=c(50,2000), X2=c(-7.5,6)))
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR2M, 
                                   FUN_CALIB = Calibration_Michel,
                                   SearchRanges = cal_ranges)

OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, 
                                   RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, 
                                   CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR2M,
                                   verbose = TRUE)

Param <- OutputsCalib$ParamFinalR

#Simulation run
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)

pdf(paste0(dir, "ter/models/river/river_predam/Ter_GR2M.pdf"))
plot(OutputsModel, Qobs = Qobs$Qmm)
dev.off()

pdf(paste0(dir, "ter/models/river/river_predam/KGE_cal.pdf"))
ggof(OutputsModel$Qsim, data_airGR$Qmm[Ind_Run], dates=data_airGR$DatesR[Ind_Run], 
     ftype = "dm", FUN = mean)
dev.off()

#save parmeters calibration for future simualtion (forecast)
save(Param, file=paste0(dir, "ter/river_up/parameter_calibration.RData"))

#basura?

#IndPeriod_WarmUp <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_ini), 
#                        which(format(InputsModel$DatesR, format = "%Y-%m-%d") == WarmUp_end))
#Ind_Run <- seq(which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_ini), 
#               which(format(InputsModel$DatesR, format = "%Y-%m-%d") == Run_end))

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

pdf(paste0(dir, "ter/models/river/river_predam/Ter_GR4J_monthly.pdf"))
plot(OutputsModel, Qobs = data_airGR$Qmm[Ind_Run])
dev.off()

#write.csv(data.frame(obs=data_airGR$Qmm[Ind_Run], sim=OutputsModel$Qsim), 
#          paste0("/home/ry4902/Documents/Workflow/Complete/Ter_GR4J_dat.txt"),
#          quote = F, row.names = F)

#pdf(paste0(dir_out, "GR4J/KGE_cal.pdf"))
pdf(paste0(dir, "ter/models/river/river_predam/KGE_cal_monthly.pdf"))
ggof(OutputsModel$Qsim, data_airGR$Qmm[Ind_Run], dates=data_airGR$DatesR[Ind_Run], 
     ftype = "dm", FUN = mean)
dev.off()

#save parmeters calibration for future simualtion (forecast)
save(Param, file=paste0(dir, "ter/models/river/river_predam/parameter_calibration_monthly.RData"))
