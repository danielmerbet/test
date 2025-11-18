library(zoo)
sqd_balance <- read.csv("~/Documents/NEREIDA/main/ter/data/water/sqd/balance_original.csv")
sqd_balance$date <- as.Date(sqd_balance$date, format="%m/%d/%Y")
past_q <- read.csv("~/Documents/NEREIDA/main/ter/data/water/pasteral/q_aca_postdam_ter.csv")
past_q$date <- as.Date(past_q$date, format="%Y/%m/%d")

sqd_past <- merge(sqd_balance, past_q, by="date")
colnames(sqd_past) <- c("date", "V", "Qi", "Qo", "Q_past")
sqd_past$Qi[sqd_past$Qi<0] <- 0
sqd_past$Qi <- na.approx(sqd_past$Qi)
sqd_past$Qo <- na.approx(sqd_past$Qo)

#send to WTP Cardedeu
past_balance <- data.frame(date=sqd_past$date,
                          Qi=sqd_past$Qo,
                          Qdiv=(sqd_past$Qo-sqd_past$Q_past))
past_balance$Qi[past_balance$Qi<0] <- 0
#past_balance$Qi <- abs(past_balance$Qi)
#past_balance$Qdiv[past_balance$Qdiv<0] <- 0
write.csv(past_balance, "~/Documents/NEREIDA/main/ter/data/water/pasteral/balance.csv", 
          quote = F, row.names = F)
