library(airGR);
library(imputeTS)
library(lubridate); library(hydroGOF)
#para tener hydroGOF: 
#paso 1: install.packages("raster"):
#paso2:install.packages("hydroGOF")
#library(hydroGOF)

dir <- "~/Documents/NEREIDA/main/"
basin <- "ter" # tordera
area <- 1661.609*1e6 #all basin: 3171660000.000 #m², alt: 1802660000.000

bal_meteo_Q <- read.csv(paste0(dir, basin,"/data/water/iwrm/data_iwrm.csv"))
bal_meteo_Q$date <- as.Date(bal_meteo_Q$date)

#data ready to input in hydrologic model
data_airGR <- data.frame(DatesR=bal_meteo_Q$date, 
                         P=bal_meteo_Q$P_Sau, 
                         T=bal_meteo_Q$T_Sau, 
                         E=bal_meteo_Q$E_Sau, 
                         Qls=bal_meteo_Q$Qin_calc_sau*1000 , 
                         Qmm=bal_meteo_Q$Qin_calc_sau*(1000*86400)/(area) ) #m3/s   m/aream2    mm   86400s/d
#l 1m3        =  m  86400s 1000mm = mm      
#s 1000l Am2  =  s    d      1m   = d
data_airGR$Qls[data_airGR$Qls<0] <- NA
data_airGR$Qmm[data_airGR$Qmm<0] <- NA
#data_airGR$Qls <- na_ma(data_airGR$Qls)
#data_airGR$Qmm <- na_ma(data_airGR$Qmm)

summary(data_airGR)

#Creating inputs for model
InputsModel <- airGR::CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, DatesR = as.POSIXct(data_airGR$DatesR),
                                 Precip = data_airGR$P, PotEvap = data_airGR$E)

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
RunOptions <- airGR::CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IndPeriod_WarmUp = IndPeriod_WarmUp)
#str(RunOptions)

#Other statistical criteria for simulation (NSE)
#InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
#                               RunOptions = RunOptions, VarObs = "Q", Obs = data_airGR$Qmm[Ind_Run])
#str(InputsCrit)

#selected statistical criteria for simulation (KGE)
InputsCrit <- airGR::CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, VarObs = "Q", 
                               Obs = data_airGR$Qmm[Ind_Run])
#str(InputsCrit)

#Calibration options (the ranges were taken according to literature)
cal_ranges <- as.matrix(data.frame(X1=c(50,600), X2=c(-2,2), X3=c(10,450), X4=c(0.25,4.35)))
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

pdf(paste0(dir, "ter/river_up/Ter_GR4J.pdf"))
plot(OutputsModel, Qobs = data_airGR$Qmm[Ind_Run])
dev.off()

#write.csv(data.frame(obs=data_airGR$Qmm[Ind_Run], sim=OutputsModel$Qsim), 
#          paste0("/home/ry4902/Documents/Workflow/Complete/Ter_GR4J_dat.txt"),
#          quote = F, row.names = F)

#pdf(paste0(dir_out, "GR4J/KGE_cal.pdf"))
pdf(paste0(dir, "ter/river_up/KGE_cal.pdf"))
ggof(OutputsModel$Qsim, data_airGR$Qmm[Ind_Run], dates=data_airGR$DatesR[Ind_Run], 
     ftype = "dm", FUN = mean)
dev.off()

#save parmeters calibration for future simualtion (forecast)
save(Param, file=paste0(dir, "ter/river_up/parameter_calibration.RData"))

#RUn the CemaNeigeGR4J (with snow module)
## preparation of the InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, DatesR = as.POSIXct(data_airGR$DatesR),
                                 Precip = data_airGR$P, PotEvap = data_airGR$E, TempMean = data_airGR$T,
                                 ZInputs = median(BasinInfo$HypsoData),
                                 HypsoData = BasinInfo$HypsoData, NLayers = 5)

## run period selection
Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d")=="1990-01-01"),
               which(format(BasinObs$DatesR, format = "%Y-%m-%d")=="1999-12-31"))


## --- original version of CemaNeige

## preparation of the RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Run)

