library(SSN)
library(plyr)
library(RODBC)
library(ggplot2)
library(reshape2)

#Gather reference site info for determining reference condition
con <- odbcConnectAccess('//deqlab1/biomon/Databases/Biomon_Phoenix.mdb')
refOG <- sqlFetch(con, 'STATION 2015_calculated')
odbcCloseAll()
ref <- refOG[!is.na(refOG$F2014_REF),]
ref <- ref[ref$F2014_REF == 'Y',]

#Identify biocriteria impaired station-samples
impaired <- read.csv('midcoast_new_status.csv')
impaired <- impaired[grep('Imp',impaired$biocriteria_status), ]

#Per Ver Hoef, re-fit using REML
fit <- glmssn(log10_FSS_26Aug14 ~ STRMPWR + PPT_1981_2010 + MIN_Z + 
                OWN_PRI_PRCA + SQM_ARCA + OWN_URB_PARCA + HDWTR, 
              EstMeth = "REML",
              ssn1,
              CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
              addfunccol = "afvArea",
              family = "Gaussian")

#Extract observation data frame to get variables to populate prediction data frame
obs <- getSSNdata.frame(fit, Name = 'Obs')
preds <- getSSNdata.frame(fit, Name = "preds")

#Fill in the model variables into the prediction data frame
# pid.order <- preds$pid
# preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
# obs_sub <- obs[,c('STATION_KEY',all.vars(fit$args$formula))]
# obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_FSS_26Aug14, decreasing = TRUE),]
# obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]
# preds <- merge(preds, obs_sub, by = 'STATION_KEY', all.x = TRUE)
# preds <- preds[match(pid.order,preds$pid),]
# row.names(preds) <- preds$pid
# fit <- putSSNdata.frame(preds, fit, Name = "preds")

#get the data frame and build the scenario
preds.0 <- getSSNdata.frame(fit, Name = "preds")
preds.0[preds.0$STATION_KEY %in% impaired$STATION_KEY, 
        'OWN_PRI_PRCA'] <- quantile(
          preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY, 'OWN_PRI_PRCA'], 
          seq(0, 1, 0.25))[4]

preds.0[preds.0$STATION_KEY %in% impaired$STATION_KEY, 
        'OWN_URB_PARCA'] <- quantile(
          preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY, 'OWN_URB_PARCA'], 
          seq(0, 1, 0.25))[4]

#Put the prediction data frame modified for the scenario into the SSN object
fit_0 <- putSSNdata.frame(preds.0, fit, Name = "preds")

#Run the prediction using the scenario
fit_0_preds <- predict.glmssn(fit_0, predpointsID = "preds", 
                                 newdata = 'preds')

#Get the prediction results
preds.0 <- getSSNdata.frame(fit_0_preds, Name = 'preds')

#Bring in the unsacaled-untransformed observed FSS for comparison
obs_sub$FSS <- 10^(obs_sub$log10_FSS_26Aug14/100 * max_log10_fss)
preds.0 <- merge(preds.0, obs_sub[,c('STATION_KEY', 'FSS')], 
                 by = 'STATION_KEY', all.x = TRUE)

#Generate prediction intervals
critval <- qnorm(0.975)
preds.0$uci <- preds.0$log10_FSS_26Aug14 + (critval * preds.0$log10_FSS_26Aug14.predSE)
preds.0$lci <- preds.0$log10_FSS_26Aug14 - (critval * preds.0$log10_FSS_26Aug14.predSE)

#Unscale and untransform the prediction for comparison to observed fss
preds.0$FSS_target <- 10^(preds.0$log10_FSS_26Aug14/100 * max_log10_fss)

#Using unscaled/untransformed FSS determine sediment stressor stations 
#based on scenario
preds.0$Sed_Stressor <- ifelse(preds.0$FSS > preds.0$FSS_target,TRUE,FALSE)

#Identify the percent reduction target
preds.0$pr_target <- abs(((preds.0$FSS_target - preds.0$FSS)/preds.0$FSS) * 100)

#Untransform the standard error of the prediction and prediction intervals
preds.0$predSE_untran <- 10^(preds.0$log10_FSS_26Aug14.predSE/100 * max_log10_fss)
preds.0$untran_uci <- 10^(preds.0$uci/100 * max_log10_fss)
preds.0$untran_lci <- 10^(preds.0$lci/100 * max_log10_fss)

#Isolate only those stations with sediment identified as stressor
preds.0 <- preds.0[!is.na(preds.0$FSS),]
ss <- preds.0[preds.0$STATION_KEY %in% impaired$STATION_KEY & preds.0$Sed_Stressor,]
ss[order(ss$pr_target),]

df <- melt(ss, id.vars = "STATION_KEY", 
           measure.vars = c('FSS', 'FSS_target', 'untran_uci', 'untran_lci'))
hline.data <- df[df$variable == 'FSS_target',]
df <- df[df$variable != 'FSS_target',]
ggplot(df, 
       aes(x = 1, y = value, color = variable)) + 
  geom_point(size = c(3, 8, 8), shape = c(19, 95, 95)) + 
  scale_color_discrete(labels = c("Observed BSTI",
                                  "UCI", "LCI")) +
  facet_wrap(~ STATION_KEY) +
  geom_hline(aes(yintercept = value), hline.data) 
