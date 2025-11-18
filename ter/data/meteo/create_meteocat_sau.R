tp_meteocat <- read.table("~/Documents/NEREIDA/main/ter/data/meteo/meteocat_canbranques.txt", header=T)
tp_meteocat <- data.frame(date=as.Date(paste(tp_meteocat$year, tp_meteocat$month, tp_meteocat$day, sep="-")),
           tp=tp_meteocat$tp)
tp1 <- tp_meteocat$tp
write.csv(tp_meteocat, "~/Documents/NEREIDA/main/ter/data/meteo/meteocat_canbranques.csv",
          quote = F, row.names = F)

tp_meteocat <- read.table("~/Documents/NEREIDA/main/ter/data/meteo/meteocat_campdevanol.txt", header=T)
tp_meteocat <- data.frame(date=as.Date(paste(tp_meteocat$year, tp_meteocat$month, tp_meteocat$day, sep="-")),
                          tp=tp_meteocat$tp)
tp2 <- tp_meteocat$tp
write.csv(tp_meteocat, "~/Documents/NEREIDA/main/ter/data/meteo/meteocat_campdevanol.csv",
          quote = F, row.names = F)

tp_meteocat <- read.table("~/Documents/NEREIDA/main/ter/data/meteo/meteocat_valldenbas.txt", header=T)
tp_meteocat <- data.frame(date=as.Date(paste(tp_meteocat$year, tp_meteocat$month, tp_meteocat$day, sep="-")),
                          tp=tp_meteocat$tp)
tp3 <- tp_meteocat$tp
write.csv(tp_meteocat, "~/Documents/NEREIDA/main/ter/data/meteo/meteocat_valldenbas.csv",
          quote = F, row.names = F)

tp_meteocat <- read.table("~/Documents/NEREIDA/main/ter/data/meteo/meteocat_vic.txt", header=T)
tp_meteocat <- data.frame(date=as.Date(paste(tp_meteocat$year, tp_meteocat$month, tp_meteocat$day, sep="-")),
                          tp=tp_meteocat$tp)
tp4 <- tp_meteocat$tp
write.csv(tp_meteocat, "~/Documents/NEREIDA/main/ter/data/meteo/meteocat_vic.csv",
          quote = F, row.names = F)

tp_mean <- data.frame(date=tp_meteocat$date, tp=rowMeans(data.frame(tp1, tp2, tp3, tp4)))
write.csv(tp_mean, "~/Documents/NEREIDA/main/ter/data/meteo/meteocat_mean.csv",
          quote = F, row.names = F)
