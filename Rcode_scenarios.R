library(SSN)
library(plyr)

options(scipen = 100)

#list of impaired stations
impaired <- data.frame(STATION_KEY = c(21842,34660,21792,33361,26818,33418,33417,34695,26822,33320,33333,30403,34665,26816,25297,26964,29906,33327),
                       TMDL_Target = c(14,14,14,3,7,rep(14,6),8,14,8,14,8,14,14))

#Get the model object
load('ssn1_glmssn_SSE.Rdata')

#Fill in the model variables into the prediction data frame
obs <- getSSNdata.frame(ssn1.glmssn.SSE, Name = 'Obs')
preds <- getSSNdata.frame(ssn1.glmssn.SSE, Name = "preds")
pid.order <- preds$pid
preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
preds <- merge(preds, ddply(obs, .(STATION_KEY), summarise, sum_1095_days = mean(sum_1095_days),
                                PALITHERODRCA = mean(PALITHERODRCA),
                                PADISRSA_1YR = mean(PADISRSA_1YR),
                                PASILTRCA = mean(PASILTRCA),
                                APOPRCA2010 = mean(APOPRCA2010),
                                POWNRCA_FED = mean(POWNRCA_FED),
                                FSS_26Aug14 = mean(FSS_26Aug14)), by = 'STATION_KEY')
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
ssn1.glmssn.SSE <- putSSNdata.frame(preds, ssn1.glmssn.SSE, Name = "preds")

#### Run a scenario with zero disturbance ####
ssn1.glmssn.SSE.0 <- ssn1.glmssn.SSE

#get the data frame and build the scenario
preds.0 <- getSSNdata.frame(ssn1.glmssn.SSE.0, Name = "preds")
preds.0$PADISRSA_1YR <- 0
preds.0 <- preds.0[match(pid.order,preds.0$pid),]
row.names(preds.0) <- preds.0$pid
ssn1.glmssn.SSE.0 <- putSSNdata.frame(preds.0, ssn1.glmssn.SSE.0, Name = "preds")

#Run the prediction
ssn1.glmssn.SSE.0.preds <- predict(ssn1.glmssn.SSE.0,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.0.dis <- getSSNdata.frame(ssn1.glmssn.SSE.0.preds, Name = 'preds')
preds.0.dis$FSS_26Aug14_untran <- 10^(preds.0.dis$log10_FSS_26Aug14)

impaired.0 <- merge(preds.0.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.0$target.met <- ifelse(impaired.0$FSS_26Aug14_untran < impaired.0$TMDL_Target,1,0)
impaired.0$pr <- (1-impaired.0$FSS_26Aug14_untran/impaired.0$FSS_26Aug14)*100
View(arrange(impaired.0,HU_8_NAME))

#### Run a scenario with zeros for dis, pop and 100 for fed_own ####
ssn1.glmssn.SSE.1 <- ssn1.glmssn.SSE

#get the data frame and build the scenario
preds.1 <- getSSNdata.frame(ssn1.glmssn.SSE.1, Name = "preds")
preds.1$PADISRSA_1YR <- 0
preds.1$APOPRCA2010 <- 0
preds.1$POWNRCA_FED <- 100
preds.1 <- preds.1[match(pid.order,preds.1$pid),]
row.names(preds.1) <- preds.1$pid
ssn1.glmssn.SSE.1 <- putSSNdata.frame(preds.1, ssn1.glmssn.SSE.1, Name = "preds")

#Run the prediction
ssn1.glmssn.SSE.1.preds <- predict(ssn1.glmssn.SSE.1,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.1.dis <- getSSNdata.frame(ssn1.glmssn.SSE.1.preds, Name = 'preds')
preds.1.dis$FSS_26Aug14_untran <- 10^(preds.1.dis$log10_FSS_26Aug14)

impaired.1 <- merge(preds.1.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.1$target.met <- ifelse(impaired.1$FSS_26Aug14_untran < impaired.1$TMDL_Target,1,0)
impaired.1$pr <- (1-impaired.1$FSS_26Aug14_untran/impaired.1$FSS_26Aug14)*100
View(arrange(impaired.1,HU_8_NAME))

#### Run a scenario with zeros for all human influence ####
ssn1.glmssn.SSE.2 <- ssn1.glmssn.SSE

#get the data frame and build the scenario
preds.2 <- getSSNdata.frame(ssn1.glmssn.SSE.2, Name = "preds")
preds.2$PADISRSA_1YR <- 0
preds.2$APOPRCA2010 <- 0
preds.2$POWNRCA_FED <- 0
preds.2 <- preds.2[match(pid.order,preds.2$pid),]
row.names(preds.2) <- preds.2$pid
ssn1.glmssn.SSE.2 <- putSSNdata.frame(preds.2, ssn1.glmssn.SSE.2, Name = "preds")

#Run the prediction
ssn1.glmssn.SSE.2.preds <- predict(ssn1.glmssn.SSE.2,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.2.dis <- getSSNdata.frame(ssn1.glmssn.SSE.2.preds, Name = 'preds')
preds.2.dis$FSS_26Aug14_untran <- 10^(preds.2.dis$log10_FSS_26Aug14)

impaired.2 <- merge(preds.2.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.2$target.met <- ifelse(impaired.2$FSS_26Aug14_untran < impaired.2$TMDL_Target,1,0)
impaired.2$pr <- (1-impaired.2$FSS_26Aug14_untran/impaired.2$FSS_26Aug14)*100
View(arrange(impaired.2,HU_8_NAME))

#### Run a scenario with completely unrealistic best case scenario ####
ssn1.glmssn.SSE.3 <- ssn1.glmssn.SSE

#get the data frame and build the scenario
preds.3 <- getSSNdata.frame(ssn1.glmssn.SSE.3, Name = "preds")
preds.3$sum_1095_days <- max(preds.3$sum_1095_days)
preds.3$PALITHERODRCA <- min(preds.3$PALITHERODRCA)
preds.3$PASILTRCA <- min(preds.3$PASILTRCA)
preds.3$PADISRSA_1YR <- 0
preds.3$APOPRCA2010 <- 0
preds.3$POWNRCA_FED <- 100
preds.3 <- preds.3[match(pid.order,preds.3$pid),]
row.names(preds.3) <- preds.3$pid
ssn1.glmssn.SSE.3 <- putSSNdata.frame(preds.3, ssn1.glmssn.SSE.3, Name = "preds")

#Run the prediction
ssn1.glmssn.SSE.3.preds <- predict(ssn1.glmssn.SSE.3,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.3.dis <- getSSNdata.frame(ssn1.glmssn.SSE.3.preds, Name = 'preds')
preds.3.dis$FSS_26Aug14_untran <- 10^(preds.3.dis$log10_FSS_26Aug14)

impaired.3 <- merge(preds.3.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.3$target.met <- ifelse(impaired.3$FSS_26Aug14_untran < impaired.3$TMDL_Target,1,0)
impaired.3$pr <- (1-impaired.3$FSS_26Aug14_untran/impaired.3$FSS_26Aug14)*100
View(arrange(impaired.3,HU_8_NAME))

#### Run a scenario with completely unrealistic best case scenario ####
ssn1.glmssn.SSE.4 <- ssn1.glmssn.SSE

#get the data frame and build the scenario
preds.4 <- getSSNdata.frame(ssn1.glmssn.SSE.4, Name = "preds")
preds.4$sum_1095_days <- 0
preds.4$PALITHERODRCA <- 0
preds.4$PASILTRCA <- 0
preds.4$PADISRSA_1YR <- 0
preds.4$APOPRCA2010 <- 0
preds.4$POWNRCA_FED <- 0
preds.4 <- preds.4[match(pid.order,preds.4$pid),]
row.names(preds.4) <- preds.4$pid
ssn1.glmssn.SSE.4 <- putSSNdata.frame(preds.4, ssn1.glmssn.SSE.4, Name = "preds")

#Run the prediction
ssn1.glmssn.SSE.4.preds <- predict(ssn1.glmssn.SSE.4,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.4.dis <- getSSNdata.frame(ssn1.glmssn.SSE.4.preds, Name = 'preds')
preds.4.dis$FSS_26Aug14_untran <- 10^(preds.4.dis$log10_FSS_26Aug14)

impaired.4 <- merge(preds.4.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.4$target.met <- ifelse(impaired.4$FSS_26Aug14_untran < impaired.4$TMDL_Target,1,0)
impaired.4$pr <- (1-impaired.4$FSS_26Aug14_untran/impaired.4$FSS_26Aug14)*100
View(arrange(impaired.4,HU_8_NAME))