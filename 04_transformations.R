# -----------------------------------------------------------
# SSN
library(SSN)
library(stringr)
library(MASS)
library(plyr)
library(reshape2)
obs.complete <- read.csv("ssn_RF_data.csv")
obs.complete$SVN <- str_trim(obs.complete$SVN)
#obs.bsti <- read.csv('fss2_s2_data_testing.csv')
#obs.bsti <- within(obs.bsti, rm(X))
#load('fss2_s2_scaled.Rdata')
obs.bsti <- bsti.s2
ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN06/lsn.ssn", 
                  predpts = "preds", o.write = TRUE)
#ssn1 <- importSSN('C:/users/pbryant/desktop/midcoasttmdl-gis/revisedssn/lsn05/lsn.ssn', o.write = TRUE)
obs<- getSSNdata.frame(ssn1, Name = "Obs")
vars <- c("STATION_KEY", "SVN", "DATE","YEAR",'BSTI')

obs.vars <- obs.complete[,vars]
obs.vars <- cbind(obs.vars, obs.bsti)

obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",
                         "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", 
                         "HU_12_NAME", "HU_08", "LONG_RAW", "NHDHigh",
                         "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", 
                         "HU_12")],
                  obs.vars, 
                  by = "SVN",
                  all.x = TRUE)
#obs.vars <- obs.vars[!is.na(obs.vars$STRMPWR),]

#RID crosswalk by SVN
#rid.cross <- merge(obs[,c('STATION_KE','SVN','rid','pid')],
#       obs.complete[,c('SVN','rid','pid')],by='SVN',suffixes=c('.NEW','.OLD'))

#Output for GIS display
#write.csv(obs.vars[,c("SVN","rid","sum_1095_days","FSS_26Aug14", 
#                      "PDISRSA_1YR", "POWNRCA_PRI", "PALITHERODRCA", 
#                      "PASILTRCA", "DAPOPRCA2010" )], 'final_model_vars.csv')

#don't run when going through. only for data exploration
#obs.vars <- arrange(obs.vars, STATION_KEY, desc(YEAR))
#obs.vars.sub <- obs.vars[!duplicated(obs.vars$STATION_KEY),]

#### Response TRANSFORMATION####
names(obs.vars)
####FSS_26Aug14####
# hist(obs.vars.sub$FSS_26Aug14)
# plot(density(obs.vars.sub$FSS_26Aug14))
# hist(log10(obs.vars.sub$FSS_26Aug14))
# plot(density(log10(obs.vars.sub$FSS_26Aug14)))
# shapiro.test(log10(obs.vars.sub$FSS_26Aug14))
# ks.test(log10(obs.vars.sub$FSS_26Aug14), "pnorm")
#log should be just fine
obs.vars$log10_BSTI <- log10(obs.vars$BSTI)

#### Variable Scaling #### 
# melted <- melt(obs.vars[,c(names(obs.bsti),'log10_FSS_26Aug14')])
# min.max <- ddply(melted, .(variable), 
#                  summarize, 
#                  min_val = min(value), 
#                  max_val = max(value))
# obs.vars[,c(names(obs.bsti),'log10_FSS_26Aug14')] <- as.data.frame(lapply(
#   obs.vars[,c(names(obs.bsti),'log10_FSS_26Aug14')], function(x) {(x-min(x))/
#       (max(x)-min(x))}))
# save(min.max, file = 'minmax.Rdata')
max_log10_bsti <- max(obs.vars$log10_BSTI, na.rm = TRUE)

obs.vars$log10_BSTI <- obs.vars$log10_BSTI / max_log10_bsti * 100

#### Put the data back into the ssn ####
#Now that we have the transformed variables we put them back in the 
#SSN object
obs.vars$HDWTR <- as.factor(obs.vars$HDWTR)
obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
row.names(obs.vars) <- obs.vars$pid
ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

obs <- getSSNdata.frame(ssn1, Name = 'Obs')
preds <- getSSNdata.frame(ssn1, Name = "preds")

pid.order <- preds$pid
preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
obs_sub <- obs[,c('STATION_KEY',pkeep,'log10_BSTI','HDWTR')]
obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_BSTI, decreasing = TRUE),]
obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]
preds <- merge(preds, obs_sub, by = 'STATION_KEY', all.x = TRUE)
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
ssn1 <- putSSNdata.frame(preds, ssn1, Name = "preds")

#save the ssn object to the github folder
#writeSSN(ssn1, filename = 'bugs.ssn')
