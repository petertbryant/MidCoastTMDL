library(reshape2)
library(plyr)

options(stringsAsFactors = FALSE)

#RUN1
vi_median_name <- "bsti_vi_median_20160701_1539.Rdata"
bsti_name <- "bsti_20160701_1539.Rdata"

#RUN2
# vi_median_name <- "bsti_vi_median_20160205_1147.Rdata"
# bsti_name <- "bsti_20160205_1147.Rdata"

load(vi_median_name)
load(bsti_name)

pcor <- cor(bsti[,setdiff(bsti.vi.median$var_name,c("DATE", "BSTI"))])


pcor[grep('FED',rownames(pcor)),grep('FED',rownames(pcor))]
#               OWN_FED_PRCA OWN_FED_PRSA OWN_FED_PARCA OWN_FED_PARSA
# OWN_FED_PRCA     1.0000000    0.9120892     0.8459522     0.6731482
# OWN_FED_PRSA     0.9120892    1.0000000     0.7368272     0.5924237
# OWN_FED_PARCA    0.8459522    0.7368272     1.0000000     0.8012213
# OWN_FED_PARSA    0.6731482    0.5924237     0.8012213     1.0000000
#All federal spatial scales are highly correlated let's try combining the
#RSA, RCA, ARSA and ARCA into a single variable and also try combining
#RSA and ARSA into a single variable. We will fit competing models and select the 
#one with the lowest RMSE and if the RMSE is the same we will use the complete ARCA

#I went back and added the SQM for each scale of influence so we can 
#recalculate different scales
obs.complete <- read.csv("ssn_RF_data.csv")

obs.complete$OWN_FED_RCA <- obs.complete$OWN_FED_PRCA * obs.complete$SQM_RCA
obs.complete$OWN_FED_RSA <- obs.complete$OWN_FED_PRSA * obs.complete$SQM_RSA
obs.complete$OWN_FED_ARCA <- obs.complete$OWN_FED_PARCA * obs.complete$SQM_ARCA
obs.complete$OWN_FED_ARSA <- obs.complete$OWN_FED_PARSA * obs.complete$SQM_ARSA

obs.complete$OWN_FED_PARCA_ALL <- rowSums(obs.complete[,
                                                      grep('OWN_FED_[^P]', 
                                                           names(obs.complete)
                                                           )]) / 
  rowSums(obs.complete[, grep('SQM_', names(obs.complete))])

obs.complete$OWN_FED_PARSA_ALL <- rowSums(obs.complete[,
                                                      grep('OWN_FED_[^P][^C][^C]', 
                                                           names(obs.complete)
                                                           )]) / 
  rowSums(obs.complete[, grep('SQM_.[^C][^C]', names(obs.complete))])

ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN06/lsn.ssn",
                  predpts = "preds", o.write = TRUE)
# ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN07/lsn.ssn",
#                   predpts = "preds", o.write = TRUE)
#only needs to be run once - RAN 07-01-2016
#createDistMat(ssn1, o.write = TRUE)
# ssn2 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN06/lsn.ssn",
#                   predpts = "preds", o.write = TRUE)
#ssn1 <- importSSN('C:/users/pbryant/desktop/midcoasttmdl-gis/revisedssn/lsn05/lsn.ssn', o.write = TRUE)
obs<- getSSNdata.frame(ssn1, Name = "Obs")

#Put together a data frame of the selected predictors and identifying information
vars <- c("STATION_KEY", "SVN", "DATE","YEAR",'BSTI','OWN_FED_PARCA_ALL',
          'OWN_FED_PARSA_ALL', 'sum_1095_days', 'POP_DARCA', 'XSLOPE_MAP',
          'MIN_Z', 'OWN_FED_PRCA', 'HDWTR')
obs.complete.vars <- obs.complete[,vars]

obs.complete.vars$log10_BSTI <- log10(obs.complete.vars$BSTI)

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
max_log10_bsti <- max(obs.complete.vars$log10_BSTI, na.rm = TRUE)

obs.complete.vars$log10_BSTI <- obs.complete.vars$log10_BSTI / max_log10_bsti * 100

#Normalize by maximum range
melted <- melt(obs.complete.vars[,names(obs.complete.vars[,-c(grep('_P',names(obs.complete.vars)),
                                                    which(names(obs.complete.vars) %in% 
                                                            c('DATE', 'STATION_KEY',
                                                              'SVN', 'YEAR', 'BSTI', 'HDWTR')))])])
min.max <- ddply(melted, .(variable), 
                 summarize, 
                 min_val = min(value), 
                 max_val = max(value))
obs.complete.vars[,-c(grep('_P',names(obs.complete.vars)),
                 which(names(obs.complete.vars) %in% 
                         c('DATE', 'STATION_KEY',
                           'SVN', 'YEAR', 'BSTI', 'HDWTR')))] <- as.data.frame(
                           lapply(obs.complete.vars[,-c(grep('_P', names(obs.complete.vars)), 
                                                   which(names(obs.complete.vars) %in% c('DATE', 'STATION_KEY',
                                                                                    'SVN', 'YEAR', 'BSTI', 'HDWTR')))], 
                                  function(x) {((x) / (max(x)))*100}))
#Merge the selected predictors with the critical columns by SVN so they match
#to the correct sample
obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",
                         "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", 
                         "HU_12_NAME", "HU_08", "NHDHigh", "LONG_RAW", 
                         "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", 
                         "HU_12")],
                  obs.complete.vars, 
                  by = "SVN",
                  all.x = TRUE)

#### Put the data back into the ssn ####
#Now that we have the transformed variables we put them back in the 
#SSN object
obs.vars$HDWTR <- as.factor(obs.vars$HDWTR)
obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
row.names(obs.vars) <- obs.vars$pid
ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

obs <- getSSNdata.frame(ssn1, Name = 'Obs')
preds <- getSSNdata.frame(ssn1, Name = "preds")
# preds <- getSSNdata.frame(ssn1, Name = "obs_prd")

pid.order <- preds$pid
preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
obs_sub <- obs[,c('STATION_KEY','OWN_FED_PARCA_ALL',
                  'OWN_FED_PARSA_ALL', 'sum_1095_days', 
                  'POP_DARCA', 'XSLOPE_MAP',
                  'MIN_Z', 'OWN_FED_PRCA', 'log10_BSTI','HDWTR')]
obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_BSTI, decreasing = TRUE),]
obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]
preds.vars <- merge(preds, obs_sub, by = 'STATION_KEY', all.x = TRUE)
# preds.vars <- merge(preds[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",
#                          "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", 
#                          "HU_12_NAME", "HU_08", "LONG_RAW", "NHDHigh",
#                          "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", 
#                          "HU_12")],
#                   obs.complete.vars, 
#                   by = "SVN",
#                   all.x = TRUE)
preds.vars <- preds.vars[match(pid.order,preds.vars$pid),]
row.names(preds.vars) <- preds.vars$pid
# ssn1 <- putSSNdata.frame(preds.vars, ssn1, Name = "obs_prd")
ssn1 <- putSSNdata.frame(preds.vars, ssn1, Name = "preds")

#Fit the models for comparison
ssn1_glmssn10 <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + 
                OWN_FED_PARCA_ALL + POP_DARCA + HDWTR,
              EstMeth = "ML",
              ssn1,
              CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
              addfunccol = "afvArea",
              family = "Gaussian")
save_name <- "ssn1_glmssn10_ML.Rdata"
save(ssn1_glmssn10, file = save_name)

ssn1_glmssn11 <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + 
                          OWN_FED_PARSA_ALL + POP_DARCA + HDWTR,
                        EstMeth = "ML",
                        ssn1,
                        CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
                        addfunccol = "afvArea",
                        family = "Gaussian")
save_name <- "ssn1_glmssn11_ML.Rdata"
save(ssn1_glmssn11, file = save_name)

FED_results <- InfoCritCompare2(model.list = list(ssn1_glmssn9, ssn1_glmssn10, ssn1_glmssn11))
#All RMSPE evaluate to 2
#Go with OWN_FED_ARCA_ALL
fit2 <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + 
                 OWN_FED_PARCA_ALL + POP_DARCA + HDWTR,
               EstMeth = "REML",
               ssn1,
               CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
               addfunccol = "afvArea",
               family = "Gaussian")