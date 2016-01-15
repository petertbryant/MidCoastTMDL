library(SSN)
library(plyr)
library(RODBC)
con <- odbcConnectAccess('//deqlab1/biomon/Databases/Biomon_Phoenix.mdb')
refOG <- sqlFetch(con, 'STATION 2015_calculated')
ref <- refOG[!is.na(refOG$F2014_REF),]
ref <- ref[ref$F2014_REF == 'Y',]
#load('ssn1_glmssn9_HWFAC_ML_20151216.Rdata')
#fit <- ssn1_glmssn_9_REML
fit <- glmssn(log10_FSS_26Aug14 ~ sum_1095_days + XSLOPE_MAP + MIN_Z + OWN_PRI_PRCA + 
                POP_DARCA + HDWTR, 
              EstMeth = "REML",
              ssn1,
              CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
              addfunccol = "afvArea",
              family = "Gaussian")
impaired <- read.csv('midcoast_new_status.csv')
impaired <- impaired[grep('Imp',impaired$biocriteria_status), ]
# impaired <- read.csv('midcoast_Updated_Status_Table.csv')
obs <- getSSNdata.frame(fit, Name = 'Obs')
preds <- getSSNdata.frame(fit, Name = "preds_up")

impaired$log10_Target <- log10(impaired$TMDL.FSS.Target)
impaired$pr_target <- abs(((impaired$TMDL.FSS.Target - impaired$FSS)/impaired$FSS) * 100)
impaired$pr_Q90_target <- abs(((impaired$Sediment.Stressor.Benchmark - impaired$FSS)/impaired$FSS) * 100)
#Fill in the model variables into the prediction data frame
pid.order <- preds$pid
preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
obs_sub <- obs[,c('STATION_KEY',all.vars(fit$args$formula))]
obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_FSS_26Aug14, decreasing = TRUE),]
obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]
preds <- merge(preds, obs_sub, by = 'STATION_KEY', all.x = TRUE)
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
fit <- putSSNdata.frame(preds, fit, Name = "preds_up")
preds_obs <- preds[,c('pid','log10_FSS_26Aug14')]
preds_obs <- rename(preds_obs, c('log10_FSS_26Aug14' = 'FSS'))

fit_preds <- predict.glmssn.na(fit, predpointsID = "preds_up", 
                              newdata = 'preds_up')
preds <- getSSNdata.frame(fit_preds, Name = 'preds_up')
preds <- merge(preds, preds_obs, by = 'pid')
critval <- qnorm(0.975)
preds$uci <- preds$log10_FSS_26Aug14 + (critval * preds$log10_FSS_26Aug14.predSE)
preds$lci <- preds$log10_FSS_26Aug14 - (critval * preds$log10_FSS_26Aug14.predSE)
preds$FSS_u <- 10^(preds$FSS/100 * max_log10_fss)
preds$fit_u <- 10^(preds$log10_FSS_26Aug14/100 * max_log10_fss)
preds$uci_u <- 10^(preds$uci/100 * max_log10_fss)
preds$lci_u <- 10^(preds$lci/100 * max_log10_fss)

ggplot(data = preds, aes(x = FSS_u, y = fit_u)) + 
  geom_point() + 
  xlim(0, 75) + 
  ylim(0, 75) +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(aes(x = FSS_u, y = uci_u), se = FALSE) +
  stat_smooth(aes(x = FSS_u, y = lci_u), se = FALSE) + 
  scale_y_continuous(limits = c(-10,100))

#get the data frame and build the scenario
preds.0 <- getSSNdata.frame(fit, Name = "preds_up")
preds.0$POP_DARCA <- 0
preds.0$OWN_PRI_PRCA <- 0
preds.0$sum_1095_days <- median(obs$sum_1095_days, na.rm = TRUE)
preds.0$XSLOPE_MAP <- preds.0$XSLOPE_MAP * 1.05
#preds.0$OWN_AGR_PARCA <- 0
preds.0 <- preds.0[match(pid.order,preds.0$pid),]
row.names(preds.0) <- preds.0$pid
fit_0 <- putSSNdata.frame(preds.0, fit, Name = "preds_up")

#Run the prediction
fit_0_preds <- predict.glmssn.na(fit_0, predpointsID = "preds_up", 
                                newdata = 'preds_up')

#Check the results
preds.0 <- getSSNdata.frame(fit_0_preds, Name = 'preds_up')
#preds.0$FSS_26Aug14_untran <- 10^(preds.0$log10_FSS_26Aug14)
preds.0$uci <- preds.0$log10_FSS_26Aug14 + (critval * preds.0$log10_FSS_26Aug14.predSE)
preds.0$lci <- preds.0$log10_FSS_26Aug14 - (critval * preds.0$log10_FSS_26Aug14.predSE)


obs_sub$FSS <- 10^(obs_sub$log10_FSS_26Aug14/100 * max_log10_fss)
preds.0 <- merge(preds.0, obs_sub[,c('STATION_KEY', 'FSS')], by = 'STATION_KEY', all.x = TRUE)
preds.0$FSS_target <- 10^(preds.0$log10_FSS_26Aug14/100 * max_log10_fss)
preds.0$Sed_Stressor <- ifelse(preds.0$FSS > preds.0$FSS_target,TRUE,FALSE)
preds.0$pr_target <- abs(((preds.0$FSS_target - preds.0$FSS)/preds.0$FSS) * 100)
preds.0$predSE_untran <- 10^(preds.0$log10_FSS_26Aug14.predSE/100 * max_log10_fss)
preds.0 <- preds.0[!is.na(preds.0$FSS),]
ss <- preds.0[preds.0$STATION_KEY %in% impaired$STATION_KEY & preds.0$Sed_Stressor,]
ss[order(ss$pr_target),]


impaired.0 <- merge(preds.0, impaired, by.x = 'STATION_KEY', by.y = 'Station.Key', all.y = TRUE)
#impaired.0 <- merge(impaired.0, preds_obs, by = 'pid')



impaired.0$log10_Target <- impaired.0$log10_Target / max_log10_fss * 100

#impaired.0$target_met <- ifelse(impaired.0$log10_FSS_26Aug14 < (log10(impaired.0$Sediment.Stressor.Benchmark)/ max_log10_fss * 100),1,0)
impaired.0$untran_FSS <- 10^(impaired.0$log10_FSS_26Aug14/100 * max_log10_fss)
impaired.0$pr_achieved <- abs(((impaired.0$untran_FSS - impaired.0$FSS)/impaired.0$FSS) * 100)
impaired.0$pr_target_met <- ifelse(impaired.0$pr_achieved >= impaired.0$pr_target,1,0)
impaired.0$pr_Q90_target_met <- ifelse(impaired.0$pr_achieved >+ impaired.0$pr_Q90_target,1,0)
impaired.0$untran_uci <- 10^(impaired.0$uci/100 * max_log10_fss)
impaired.0$untran_lci <- 10^(impaired.0$lci/100 * max_log10_fss)
impaired.0$lci_meets <- ifelse(impaired.0$lci < impaired.0$log10_Target, 1, 0)
#impaired.0[,c('STATION_KEY','log10_FSS_26Aug14','uci','lci','log10_Target','untran_FSS','TMDL.FSS.Target','target_met',all.vars(ssn1_glmssn5$args$formula))]
impaired.0 <- merge(impaired.0, preds_obs, by = 'pid', all.x = TRUE)

df <- melt(impaired.0[,], id.vars = "STATION_KEY",
     measure.vars = c('untran_FSS', 'untran_uci', 'untran_lci', 'FSS.y'))

hline.data <- data.frame(z = impaired.0[,"TMDL.FSS.Target"], 
                         STATION_KEY = impaired.0[,c('STATION_KEY')])
impaired.0$STATION_KEY %in% c(21792,37179,25297)
ggplot(df, 
       aes(x = 1, y = value, color = variable)) + 
  geom_point(size = c(3, 8, 8, 3), shape = c(19, 95, 95, 19)) + 
  scale_color_discrete(labels = c("Predicted BSTI at 0 Anthro", "UCI", "LCI", "Observed BSTI")) +
  facet_wrap(~ STATION_KEY) +
  geom_hline(aes(yintercept = z), hline.data) 

impaired.0[,c('STATION_KEY',grep('untran',names(impaired.0), value = TRUE),
              'TMDL.FSS.Target','Sediment.Stressor.Benchmark','target_met')]
