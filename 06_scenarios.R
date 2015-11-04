load('ssn1_glmssn5_20151019.Rdata')
impaired <- read.csv('midcoast_Updated_Status_Table.csv')
obs <- getSSNdata.frame(ssn1_glmssn5, Name = 'Obs')
preds <- getSSNdata.frame(ssn1_glmssn5, Name = "preds_up")

impaired$log10_Target <- log10(impaired$TMDL.FSS.Target)
impaired$pr_target <- abs(((impaired$TMDL.FSS.Target - impaired$FSS)/impaired$FSS) * 100)
#Fill in the model variables into the prediction data frame
pid.order <- preds$pid
preds <- rename(preds, c('STATION_KE' = "STATION_KEY"))
obs_sub <- obs[,c('STATION_KEY',all.vars(ssn1_glmssn5$args$formula))]
obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_FSS_26Aug14, decreasing = TRUE),]
obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]
preds <- merge(preds, obs_sub, by = 'STATION_KEY')
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
ssn1_glmssn5 <- putSSNdata.frame(preds, ssn1_glmssn5, Name = "preds_up")

ssn1_glmssn5_preds <- predict(ssn1_glmssn5, predpointsID = "preds_up", 
                              newdata = 'preds_up')
preds <- getSSNdata.frame(ssn1_glmssn5_preds, Name = 'preds_up')

critval <- qnorm(0.975)
preds$uci <- preds$log10_FSS_26Aug14 + (critval * preds$log10_FSS_26Aug14.predSE)
preds$lci <- preds$log10_FSS_26Aug14 - (critval * preds$log10_FSS_26Aug14.predSE)


#get the data frame and build the scenario
preds.0 <- getSSNdata.frame(ssn1_glmssn5, Name = "preds_up")
preds.0$DIS_1YR_PARSA <- 0
preds.0$OWN_PRI_PRCA <- 0
preds.0$OWN_AGR_PARCA <- 0
preds.0 <- preds.0[match(pid.order,preds.0$pid),]
row.names(preds.0) <- preds.0$pid
ssn1_glmssn5_0 <- putSSNdata.frame(preds.0, ssn1_glmssn5, Name = "preds_up")

#Run the prediction
ssn1_glmssn5_0_preds <- predict(ssn1_glmssn5_0, predpointsID = "preds_up", 
                                newdata = 'preds_up')

#Check the results
preds.0 <- getSSNdata.frame(ssn1_glmssn5_0_preds, Name = 'preds_up')
preds.0$FSS_26Aug14_untran <- 10^(preds.0.dis$log10_FSS_26Aug14)
preds.0$uci <- preds.0$log10_FSS_26Aug14 + (critval * preds.0$log10_FSS_26Aug14.predSE)
preds.0$lci <- preds.0$log10_FSS_26Aug14 - (critval * preds.0$log10_FSS_26Aug14.predSE)

impaired.0 <- merge(preds.0, impaired, by.x = 'STATION_KEY', by.y = 'Station.Key', all.y = TRUE)
impaired.0$target_met <- ifelse(impaired.0$log10_FSS_26Aug14 < impaired.0$log10_Target,1,0)
impaired.0$untran_FSS <- 10^impaired.0$log10_FSS_26Aug14
impaired.0$untran_uci <- 10^impaired.0$uci
impaired.0$untran_lci <- 10^impaired.0$lci
impaired.0[,c('STATION_KEY','log10_FSS_26Aug14','uci','lci','log10_Target','untran_FSS','TMDL.FSS.Target','target_met',all.vars(ssn1_glmssn5$args$formula))]

df <- melt(impaired.0[impaired.0$STATION_KEY %in% c(21792,37179,25297),], id.vars = "STATION_KEY",
     measure.vars = c('untran_FSS', 'untran_uci', 'untran_lci'))

hline.data <- data.frame(z = impaired.0[impaired.0$STATION_KEY %in% c(21792,37179,25297),c('TMDL.FSS.Target')], 
                         STATION_KEY = impaired.0[impaired.0$STATION_KEY %in% c(21792,37179,25297),c('STATION_KEY')])

ggplot(df, 
       aes(x = 1, y = value, color = variable)) + 
  geom_point(size = c(3, 8, 8), shape = c(19, 95, 95)) + 
  facet_wrap(~ STATION_KEY) +
  geom_hline(aes(yintercept = z), hline.data)

impaired.0[,c('STATION_KEY',grep('untran',names(impaired.0), value = TRUE),
              'TMDL.FSS.Target','Sediment.Stressor.Benchmark','target_met')]
