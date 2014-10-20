# -----------------------------------------------------------
# SSN
library(SSN)
ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/lsn.ssn", o.write=FALSE)
obs<- getSSNdata.frame(ssn1, Name = "Obs")
obs.complete <- read.csv("ssn_RF_data.csv")
obs.fss2 <- read.csv('fss2_s2_data.csv')
obs.fss2 <- within(obs.fss2, rm(X))

vars <- c("STATION_KEY", "SITE_NAME", names(obs.fss2))

obs.vars <- obs.complete[,vars]

obs.vars <- merge(obs.vars, 
                  obs[,c("STATION_KE","rid", "ratio", "locID", "netID", "pid", "upDist",  "afvArea",
                                "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", "HU_12_NAME", "HU_08", "LONG_RAW", "LAT_RAW", "NHDHigh",
                                "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", "HU_12")],
                  by.x = 'STATION_KEY',
                  by.y = 'STATION_KE',
                  all.x = TRUE)

putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')
