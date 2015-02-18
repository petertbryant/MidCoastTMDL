library(SSN)
library(plyr)
library(hydroGOF)

options(scipen = 100)

#list of impaired stations
impaired <- data.frame(STATION_KEY = c(21842,34660,21792,33361,26818,33418,33417,34695,26822,33320,33333,30403,34665,26816,25297,26964,29906,33327),
                       TMDL_Target = c(14,14,14,3,7,rep(14,6),8,14,8,14,8,14,14))
impaired$TMDL_Target_Scaled_log <- (log10(impaired$TMDL_Target)-min.max[min.max$variable == 'log10_FSS_26Aug14','min_val'])/(min.max[min.max$variable == 'log10_FSS_26Aug14','max_val']-min.max[min.max$variable == 'log10_FSS_26Aug14','min_val'])

#Get the model object
load('ssn1_glmssn_EEE.Rdata')

#Fill in the model variables into the prediction data frame
obs <- getSSNdata.frame(ssn1.glmssn.EEE, Name = 'Obs')
preds <- getSSNdata.frame(ssn1.glmssn.EEE, Name = "preds")
pid.order <- preds$pid
preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
preds <- merge(preds, ddply(obs, .(STATION_KEY), summarise, sum_1095_days = mean(sum_1095_days),
                                PALITHERODRCA = mean(PALITHERODRCA),
                                PDISRSA_1YR = mean(PDISRSA_1YR),
                                PASILTRCA = mean(PASILTRCA),
                                DAPOPRCA2010 = mean(DAPOPRCA2010),
                                POWNRCA_PRI = mean(POWNRCA_PRI),
                                log10_FSS_26Aug14 = mean(log10_FSS_26Aug14)), by = 'STATION_KEY')
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
ssn1.glmssn.EEE <- putSSNdata.frame(preds, ssn1.glmssn.EEE, Name = "preds")

#Generate prediction intervals
pred.EEE <- predict(ssn1.glmssn.EEE, predpointsID = "preds", newdata = 'preds')
pred.int <- getSSNdata.frame(pred.EEE, Name = 'preds')

pred.int$untran_FSS_predict <- 10^(pred.int$log10_FSS_26Aug14*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val'])
pred.int$untran_FSS <- 10^(preds$log10_FSS_26Aug14*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val'])
rmse(pred.int$untran_FSS_predict,pred.int$untran_FSS)
#3.42724

critval <- qnorm(0.975)
# pred.int$uci <- pred.int$log10_FSS_26Aug14 + (critval * pred.int$log10_FSS_26Aug14.predSE)
# pred.int$lci <- pred.int$log10_FSS_26Aug14 - (critval * pred.int$log10_FSS_26Aug14.predSE)
# pred.int$fit <- pred.int$log10_FSS_26Aug14
# pred.int$log10_FSS_26Aug14 <- preds$log10_FSS_26Aug14
preds$unscaled_log10_FSS <- preds$log10_FSS_26Aug14*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']
preds$unscaled_fit <- pred.int$log10_FSS_26Aug14*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']
preds$unscaled_predSE <- pred.int$log10_FSS_26Aug14.predSE*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']
preds$uci <- preds$unscaled_log10_FSS + (critval * preds$unscaled_predSE)
preds$lci <- preds$unscaled_log10_FSS - (critval * preds$unscaled_predSE)
preds$untran_FSS <- 10^preds$unscaled_log10_FSS
preds$untran_fit <- 10^preds$unscaled_fit
preds$untran_uci <- 10^preds$uci
preds$untran_lci <- 10^preds$lci
# link <- gaussian(link = 'identity')
# fit2 <- link$linkinv(fit)
# upr2 <- link$linkinv(upr)
# lwr2 <- link$linkinv(lwr)

#Generate smooth confidence interval lines
mean.df <- ddply(preds, .(log10_FSS_26Aug14), summarize, mean_fit = mean(fit), mean_uci = mean(uci), mean_lci = mean(lci))
ul <- loess.smooth(mean.df$log10_FSS_26Aug14, mean.df$mean_uci)
ll <- loess.smooth(mean.df$log10_FSS_26Aug14, mean.df$mean_lci)
i.for <- order(ul$y)
i.back <- order(ul$y, decreasing = TRUE)
x.p <- c(ul$x[i.for],ll$x[i.back])
y.p <- c(ul$y[i.for],ll$y[i.back])

#Make the plot
png('predInterval.png',width = 6, height = 6, res = 100, units = "in")
plot.new()
polygon(x.p,y.p,col = "light grey", border = FALSE)
par(new = TRUE)
plot(preds$log10_FSS_26Aug14, preds$fit, pch = 19, ylim = c(0,1), xlab = 'Scaled log10 FSS', ylab = 'Scaled log10 FSS Predicted')
lines(loess.smooth(mean.df$log10_FSS_26Aug14, mean.df$mean_fit))
lines(loess.smooth(mean.df$log10_FSS_26Aug14, mean.df$mean_uci))
lines(loess.smooth(mean.df$log10_FSS_26Aug14, mean.df$mean_lci))
dev.off()

#with untran
mean.df <- ddply(preds, .(untran_FSS), summarize, mean_fit = mean(untran_fit), mean_uci = mean(untran_uci), mean_lci = mean(untran_lci))
ul <- loess.smooth(mean.df$untran_FSS, mean.df$mean_uci)
ll <- loess.smooth(mean.df$untran_FSS, mean.df$mean_lci)
i.for <- order(ul$y)
i.back <- order(ul$y, decreasing = TRUE)
x.p <- c(ul$x[i.for],ll$x[i.back])
y.p <- c(ul$y[i.for],ll$y[i.back])

png('predInterval_untran.png',width = 6, height = 6, res = 100, units = "in")
plot.new()
par(new = TRUE)
plot(preds$untran_FSS, preds$untran_fit, pch = 19, ylim = c(0,70), xlab = 'FSS', ylab = 'FSS Predicted')
polygon(x.p,y.p,col="grey", border = FALSE)
par(new = TRUE)
plot(preds$untran_FSS, preds$untran_fit, pch = 19, ylim = c(0,70), xlab = 'FSS', ylab = 'FSS Predicted')
lines(loess.smooth(mean.df$untran_FSS, mean.df$mean_fit))
lines(loess.smooth(mean.df$untran_FSS, mean.df$mean_uci))
lines(loess.smooth(mean.df$untran_FSS, mean.df$mean_lci))
dev.off()


#### Look at lci in relation to TMDL_Target ####
tmp <- merge(preds, impaired, all.y = T, by ="STATION_KEY")
tmp$imp <- ifelse(tmp$lci < tmp$TMDL_Target_Scaled_log,0,1)
tmp


#### Solve for human variables ####
A <-preds$fit[415] - ssn1.glmssn.EEE$estimates$betahat[1] + ssn1.glmssn.EEE$estimates$betahat[2]*preds$sum_1095_days[415] - 
      ssn1.glmssn.EEE$estimates$betahat[5]*preds$PALITHERODRCA[415] - ssn1.glmssn.EEE$estimates$betahat[6]*preds$PASILTRCA[415]
b <- array(c(ssn1.glmssn.EEE$estimates$betahat[3],ssn1.glmssn.EEE$estimates$betahat[4],ssn1.glmssn.EEE$estimates$betahat[7]), dim=c(1,3))

#### Run a scenario with zero disturbance ####
ssn1.glmssn.EEE.0 <- ssn1.glmssn.EEE

#get the data frame and build the scenario
preds.0 <- getSSNdata.frame(ssn1.glmssn.EEE.0, Name = "preds")
preds.0$PADISRSA_1YR <- 0
preds.0 <- preds.0[match(pid.order,preds.0$pid),]
row.names(preds.0) <- preds.0$pid
ssn1.glmssn.EEE.0 <- putSSNdata.frame(preds.0, ssn1.glmssn.EEE.0, Name = "preds")

#Run the prediction
ssn1.glmssn.EEE.0.preds <- predict(ssn1.glmssn.EEE.0,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.0.dis <- getSSNdata.frame(ssn1.glmssn.EEE.0.preds, Name = 'preds')
preds.0.dis$FSS_26Aug14_untran <- 10^(preds.0.dis$log10_FSS_26Aug14)

impaired.0 <- merge(preds.0.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.0$target.met <- ifelse(impaired.0$FSS_26Aug14_untran < impaired.0$TMDL_Target,1,0)
impaired.0$pr <- (1-impaired.0$FSS_26Aug14_untran/impaired.0$FSS_26Aug14)*100
View(arrange(impaired.0,HU_8_NAME))

#### Run a scenario with zeros for human influence ####
ssn1.glmssn.EEE.1 <- ssn1.glmssn.EEE

#get the data frame and build the scenario
preds.1 <- getSSNdata.frame(ssn1.glmssn.EEE.1, Name = "preds")
preds.1$PDISRSA_1YR <- 0
preds.1$DAPOPRCA2010 <- 0
preds.1$POWNRCA_PRI <- 0
preds.1 <- preds.1[match(pid.order,preds.1$pid),]
row.names(preds.1) <- preds.1$pid
ssn1.glmssn.EEE.1 <- putSSNdata.frame(preds.1, ssn1.glmssn.EEE.1, Name = "preds")

#Run the prediction
ssn1.glmssn.EEE.1.preds <- predict(ssn1.glmssn.EEE.1,predpointsID = "preds", newdata = 'preds')

#Check the results - NEED TO RESCALE!
preds.1.dis <- getSSNdata.frame(ssn1.glmssn.EEE.1.preds, Name = 'preds')
preds.1.dis$FSS_26Aug14_untran <- 10^(preds.1.dis $log10_FSS_26Aug14*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val'])


impaired.1 <- merge(preds.1.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.1$target.met <- ifelse(impaired.1$FSS_26Aug14_untran < impaired.1$TMDL_Target,1,0)
impaired.1$pr <- (1-impaired.1$TMDL_Target/impaired.1$FSS_26Aug14_untran)*100
View(arrange(impaired.1,HU_8_NAME))

#### Run a scenario with zeros for all human influence ####
ssn1.glmssn.EEE.2 <- ssn1.glmssn.EEE

#get the data frame and build the scenario
preds.2 <- getSSNdata.frame(ssn1.glmssn.EEE.2, Name = "preds")
preds.2$PADISRSA_1YR <- 0
preds.2$APOPRCA2010 <- 0
preds.2$POWNRCA_FED <- 0
preds.2$sum_1095_days <- 0.9*preds.2$sum_1095_days
preds.2 <- preds.2[match(pid.order,preds.2$pid),]
row.names(preds.2) <- preds.2$pid
ssn1.glmssn.EEE.2 <- putSSNdata.frame(preds.2, ssn1.glmssn.EEE.2, Name = "preds")

#Run the prediction
ssn1.glmssn.EEE.2.preds <- predict(ssn1.glmssn.EEE.2,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.2.dis <- getSSNdata.frame(ssn1.glmssn.EEE.2.preds, Name = 'preds')
preds.2.dis$FSS_26Aug14_untran <- 10^(preds.2.dis$log10_FSS_26Aug14)

impaired.2 <- merge(preds.2.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.2$target.met <- ifelse(impaired.2$FSS_26Aug14_untran < impaired.2$TMDL_Target,1,0)
impaired.2$pr <- (1-impaired.2$FSS_26Aug14_untran/impaired.2$FSS_26Aug14)*100
View(arrange(impaired.2,HU_8_NAME))

#### Run a scenario with completely unrealistic best case scenario ####
ssn1.glmssn.EEE.3 <- ssn1.glmssn.EEE

#get the data frame and build the scenario
preds.3 <- getSSNdata.frame(ssn1.glmssn.EEE.3, Name = "preds")
preds.3$sum_1095_days <- max(preds.3$sum_1095_days)
preds.3$PALITHERODRCA <- min(preds.3$PALITHERODRCA)
preds.3$PASILTRCA <- min(preds.3$PASILTRCA)
preds.3$PADISRSA_1YR <- 0
preds.3$APOPRCA2010 <- 0
preds.3$POWNRCA_FED <- 100
preds.3 <- preds.3[match(pid.order,preds.3$pid),]
row.names(preds.3) <- preds.3$pid
ssn1.glmssn.EEE.3 <- putSSNdata.frame(preds.3, ssn1.glmssn.EEE.3, Name = "preds")

#Run the prediction
ssn1.glmssn.EEE.3.preds <- predict(ssn1.glmssn.EEE.3,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.3.dis <- getSSNdata.frame(ssn1.glmssn.EEE.3.preds, Name = 'preds')
preds.3.dis$FSS_26Aug14_untran <- 10^(preds.3.dis$log10_FSS_26Aug14)

impaired.3 <- merge(preds.3.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.3$target.met <- ifelse(impaired.3$FSS_26Aug14_untran < impaired.3$TMDL_Target,1,0)
impaired.3$pr <- (1-impaired.3$FSS_26Aug14_untran/impaired.3$FSS_26Aug14)*100
View(arrange(impaired.3,HU_8_NAME))

#### Run a scenario with completely unrealistic best case scenario ####
ssn1.glmssn.EEE.4 <- ssn1.glmssn.EEE

#get the data frame and build the scenario
preds.4 <- getSSNdata.frame(ssn1.glmssn.EEE.4, Name = "preds")
preds.4$sum_1095_days <- 0
preds.4$PALITHERODRCA <- 0
preds.4$PASILTRCA <- 0
preds.4$PADISRSA_1YR <- 0
preds.4$APOPRCA2010 <- 0
preds.4$POWNRCA_FED <- 0
preds.4 <- preds.4[match(pid.order,preds.4$pid),]
row.names(preds.4) <- preds.4$pid
ssn1.glmssn.EEE.4 <- putSSNdata.frame(preds.4, ssn1.glmssn.EEE.4, Name = "preds")

#Run the prediction
ssn1.glmssn.EEE.4.preds <- predict(ssn1.glmssn.EEE.4,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.4.dis <- getSSNdata.frame(ssn1.glmssn.EEE.4.preds, Name = 'preds')
preds.4.dis$FSS_26Aug14_untran <- 10^(preds.4.dis$log10_FSS_26Aug14)

impaired.4 <- merge(preds.4.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.4$target.met <- ifelse(impaired.4$FSS_26Aug14_untran < impaired.4$TMDL_Target,1,0)
impaired.4$pr <- (1-impaired.4$FSS_26Aug14_untran/impaired.4$FSS_26Aug14)*100
View(arrange(impaired.4,HU_8_NAME))

#### Run a scenario for a specific site ####
ssn1.glmssn.EEE.5 <- ssn1.glmssn.EEE

#get the data frame and build the scenario
preds.5 <- getSSNdata.frame(ssn1.glmssn.EEE.5, Name = "preds")
# preds.5$sum_1095_days <- 0
# preds.5$PALITHERODRCA <- 0
# preds.5$PASILTRCA <- 0
preds.5$PDISRSA_1YR <- 0.1379148
preds.5$DAPOPRCA2010 <- 0
preds.5$POWNRCA_PRI <- 0.1725864
preds.5 <- preds.5[match(pid.order,preds.5$pid),]
row.names(preds.5) <- preds.5$pid
ssn1.glmssn.EEE.5 <- putSSNdata.frame(preds.5, ssn1.glmssn.EEE.5, Name = "preds")

#Run the prediction
ssn1.glmssn.EEE.5.preds <- predict(ssn1.glmssn.EEE.5,predpointsID = "preds", newdata = 'preds')

#Check the results
preds.5.dis <- getSSNdata.frame(ssn1.glmssn.EEE.5.preds, Name = 'preds')
preds.5.dis[preds.5.dis$pid == 1030,]
preds.5.dis$FSS_26Aug14_untran <- 10^(preds.5.dis$log10_FSS_26Aug14)

impaired.5 <- merge(preds.5.dis, impaired, by = 'STATION_KEY', all.y = TRUE)
impaired.5$target.met <- ifelse(impaired.5$FSS_26Aug14_untran < impaired.5$TMDL_Target,1,0)
impaired.5$pr <- (1-impaired.5$FSS_26Aug14_untran/impaired.5$FSS_26Aug14)*100
View(arrange(impaired.5,HU_8_NAME))