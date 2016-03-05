library(SSN)
library(plyr)
library(RODBC)
library(ggplot2)
library(reshape2)

# #This workspace makes it so you don't have to re-run steps 02-04.
# load('06_scenarios_inputs_02192016.Rdata')
# 
# rm(list = setdiff(ls(), c('ssn1', 'max_log10_bsti', 'min.max')))
# 
# #Gather reference site info for determining reference condition
# con <- odbcConnectAccess('//deqlab1/biomon/Databases/Biomon_Phoenix.mdb')
# refOG <- sqlFetch(con, 'STATION 2015_calculated')
# odbcCloseAll()
# ref <- refOG[!is.na(refOG$F2014_REF),]
# ref <- ref[ref$F2014_REF == 'Y',]
# 
# #Saving this here makes it so we don't have to run 32 bit windows but can still
# #use the ref information from Biomon_Phoenix
# save.image('06_scenarios_inputs_02192106_0848.Rdata')
# load('06_scenarios_inputs_02192106_0848.Rdata')
# 
# # RUN 2
# fit <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + KFACT_MARCA + 
#                 OWN_FED_PRCA + DIS_3YR_PRSA + ROADLEN_DRSA + OWN_URB_PARCA + HDWTR, 
#               EstMeth = "REML",
#               ssn1,
#               CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
#               addfunccol = "afvArea",
#               family = "Gaussian")
#  
# save.image('06_scenarios_post_fit_02192016_0848.Rdata')
load('06_scenarios_post_fit_02192016_0848.Rdata')

# #RUN 3
# fit <- glmssn(log10_BSTI ~ STRMPWR + EROD_PARCA + MIN_Z + OWN_FED_PRCA + 
#                 SQM_ARCA + OWN_URB_PARCA + ROADLEN_DARCA, 
#               EstMeth = "REML",
#               ssn1,
#               CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
#               addfunccol = "afvArea",
#               family = "Gaussian")
 
CART_imp <- read.csv('midcoast_Updated_Status_Table.csv')
impaired <- read.csv('midcoast_new_status.csv')
impaired <- impaired[grep('Imp',impaired$biocriteria_status), ]
impaired <- impaired[order(impaired$STATION_KEY, decreasing = TRUE),]
#write.csv(impaired, 'mc_biocrite_impaired.csv')
obs <- getSSNdata.frame(fit, Name = 'Obs')
preds <- getSSNdata.frame(fit, Name = "preds")

#Takes the highest BSTI score from each station
obs_sub <- obs[,c('STATION_KEY',all.vars(fit$args$formula))]
obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_BSTI, decreasing = TRUE),]
obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]

#Preserve observed BSTI values in the prediction data set
preds_obs <- preds[,c('pid','log10_BSTI')]
preds_obs <- rename(preds_obs, c('log10_BSTI' = 'log10_BSTI_obs'))

#Check model fit with unmodified prediction variables
fit_preds <- predict.glmssn(fit, predpointsID = "preds", 
                              newdata = 'preds')
preds <- getSSNdata.frame(fit_preds, Name = 'preds')
preds <- merge(preds, preds_obs, by = 'pid')
critval <- qnorm(0.975)
preds$uci <- preds$log10_BSTI + (critval * preds$log10_BSTI.predSE)
preds$lci <- preds$log10_BSTI - (critval * preds$log10_BSTI.predSE)
preds$BSTI_u <- 10^(preds$log10_BSTI_obs/100 * max_log10_bsti)
preds$fit_u <- 10^(preds$log10_BSTI/100 * max_log10_bsti)
preds$uci_u <- 10^(preds$uci/100 * max_log10_bsti)
preds$lci_u <- 10^(preds$lci/100 * max_log10_bsti)

ggplot(data = preds, aes(x = BSTI_u, y = fit_u)) + 
  geom_point() + 
  xlim(0, 75) + 
  ylim(0, 75) +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(aes(x = BSTI_u, y = uci_u), se = FALSE) +
  stat_smooth(aes(x = BSTI_u, y = lci_u), se = FALSE) + 
  scale_y_continuous(limits = c(-10,100))

#Generate predictions at critical conditions for human influence
#and at the observed rainfall amounts at each station individually
for (i in 1:length(unique(preds$STATION_KEY))) {
  preds.0 <- getSSNdata.frame(fit, Name = "preds")
  
  stn <- as.character(unique(preds$STATION_KEY)[i])
  
  #RUN1
  preds.0[preds.0$STATION_KEY == stn, 
          'ROADLEN_DRSA'] <- quantile(
            preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'ROADLEN_DRSA'], 
            seq(0,1,.25))[4]
  preds.0[preds.0$STATION_KEY == stn, 
          'OWN_URB_PARCA'] <- quantile(
            preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'OWN_URB_PARCA'], 
            seq(0,1,.25))[2]
  
  ref_DIS_3YR_PRSA <- quantile(
    preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'DIS_3YR_PRSA'], 
    seq(0,1,.25))[4]
  preds.0[preds.0$STATION_KEY == stn & preds.0$DIS_3YR_PRSA > ref_DIS_3YR_PRSA, 
          'DIS_3YR_PRSA'] <- ref_DIS_3YR_PRSA
    
  ref_OWN_FED_PRCA <- quantile(
    preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'OWN_FED_PRCA'], 
    seq(0,1,.25))[2]
  preds.0[preds.0$STATION_KEY == stn & preds.0$OWN_FED_PRCA < ref_OWN_FED_PRCA, 
          'OWN_FED_PRCA'] <- ref_OWN_FED_PRCA
  
  #Run the fit
  fit_0 <- putSSNdata.frame(preds.0, fit, Name = "preds")
  
  #Run the prediction
  fit_0_preds <- predict.glmssn(fit_0, predpointsID = "preds", 
                                newdata = 'preds')
  
  #Check the results
  preds.0 <- getSSNdata.frame(fit_0_preds, Name = 'preds')
  
  #Extract prediction at selected station
  tmp <- preds.0[preds.0$STATION_KEY == stn,]
  
  #Build data frame to use for sediment stressor identification
  if (i == 1) {
    ssid <- tmp
  } else {
    ssid <- rbind(ssid, tmp)
  }
}

#SSN for target Identification
ssid <- merge(ssid, preds_obs, by = 'pid')
ssid$BSTI <- as.integer(10^(ssid$log10_BSTI_obs/100 * max_log10_bsti))
ssid$BSTI_target <- as.integer(round(10^(ssid$log10_BSTI/100 * max_log10_bsti)))
ssid$Sed_Stressor <- ifelse(ssid$BSTI > ssid$BSTI_target,TRUE,FALSE)
ssid$pr_target <- round(abs(((ssid$BSTI_target - ssid$BSTI)/ssid$BSTI) * 100),1)
ssid$predSE_untran <- 10^(ssid$log10_BSTI.predSE/100 * max_log10_bsti)
ssid$STATION_KEY <- as.character(ssid$STATION_KEY)
ssid <- ssid[order(ssid$STATION_KEY, decreasing = TRUE),]
write.csv(ssid, 'sedstressor_ssn_identified.csv', row.names = FALSE)
# ssid$untran_uci <- 10^(ssid$uci/100 * max_log10_bsti)
# ssid$untran_lci <- 10^(ssid$lci/100 * max_log10_bsti)
#ssid <- ssid[!is.na(ssid$BSTI),]
ss <- ssid[ssid$STATION_KEY %in% impaired$STATION_KEY & ssid$Sed_Stressor,]
ss[order(ss$pr_target),]

lapply(list('sum_1095_days', 'XSLOPE_MAP', 'MIN_Z', 'KFACT_MARCA',
             'ROADLEN_DRSA', 'HDWTR'), function(x) {preds[,x] <- (preds[,x]/100) * 
               min.max[min.max$variable == x, 'max_val']})

unique(ss$HU_10_NAME)

CART_HUCs <- preds[preds$STATION_KEY %in% CART_imp$Station.Key, 
                   c('STATION_KEY',grep("HU",names(preds),value = TRUE))]
unique(ss$HU_12_NAME)[!unique(ss$HU_12_NAME) %in% unique(CART_HUCs$HU_12_NAME)]

#### Generate target range ####
load('C:/users/pbryant/desktop/midcoasttmdl-gis/precip_daily_sum_1095_days.Rdata')
dfdall$sum_1095_days <- dfdall$sum_1095_days / 
  (min.max[min.max$variable == 'sum_1095_days', 'max_val'])*100

#Calculate rainfall quantiles at 10th percentiles
#rfq <- quantile(obs_sub$sum_1095_days, seq(0,1,.1))

#Run predictions using critical condition values for human variables and at 
#each 10th percentile of rainfall
for (j in 1:length(unique(ss$STATION_KEY))) {
  stn <- unique(ss$STATION_KEY)[j]

  rfq <- quantile(dfdall[dfdall$STATION_KEY == stn, 'sum_1095_days'], 
                  seq(0, 1, 0.01))

  #Set up list of dataframes to score prediction outputs based on percentile distribution
  #df_list <- list('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')
  df_list <- as.list(names(rfq))
  names(df_list) <- names(rfq)
  df_list <- lapply(df_list, function(x) X <- NULL)
  
    for (i in 1:length(df_list)) {
    #get the data frame and build the scenario
    preds.0 <- getSSNdata.frame(fit, Name = "preds")
    
    #RUN1
    preds.0[preds.0$STATION_KEY == stn, 
            'ROADLEN_DRSA'] <- quantile(
              preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'ROADLEN_DRSA'], 
              seq(0,1,.25))[4]
    preds.0[preds.0$STATION_KEY == stn, 
            'OWN_URB_PARCA'] <- quantile(
              preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'OWN_URB_PARCA'], 
              seq(0,1,.25))[2]
    
    ref_DIS_3YR_PRSA <- quantile(
      preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'DIS_3YR_PRSA'], 
      seq(0,1,.25))[4]
    preds.0[preds.0$STATION_KEY == stn & preds.0$DIS_3YR_PRSA > ref_DIS_3YR_PRSA, 
            'DIS_3YR_PRSA'] <- ref_DIS_3YR_PRSA
    
    ref_OWN_FED_PRCA <- quantile(
      preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'OWN_FED_PRCA'], 
      seq(0,1,.25))[2]
    preds.0[preds.0$STATION_KEY == stn & 
              preds.0$OWN_FED_PRCA < ref_OWN_FED_PRCA, 
            'OWN_FED_PRCA'] <- ref_OWN_FED_PRCA
    
    preds.0[preds.0$STATION_KEY == stn, 'sum_1095_days'] <- rfq[i]
    
    #Run the fit
    fit_0 <- putSSNdata.frame(preds.0, fit, Name = "preds")
    
    #Run the prediction
    fit_0_preds <- predict.glmssn(fit_0, predpointsID = "preds", 
                                  newdata = 'preds')
    
    #Check the results
    preds.0 <- getSSNdata.frame(fit_0_preds, Name = 'preds')
    
    #Put it into dataframe list
    df_list[i][[1]] <- preds.0[preds.0$STATION_KEY == stn,]
  }
  #Untransform bsti and keep track of which percentile the prediction was made with
  for (i in 1:length(df_list)) {
    quant <- names(df_list[i])
    df_list[[i]]$quant <- quant 
    df_list[[i]]$bsti_untran <- 10^(df_list[[i]]$log10_BSTI/100 * max_log10_bsti)
  }
  
  #Convert list to a dataframe to prepare for plotting
  tmp_quant <- ldply(df_list)
  
  if (j == 1) {
    df_quant <- tmp_quant
  } else {
    df_quant <- rbind(df_quant, tmp_quant)
  }
}

#Percent reduction 
rf <- obs[obs$STATION_KEY == 21792, 'sum_1095_days']
bsti_target <- (10^(df_bm[31,3] + df_bm[31,4] * rf))
bsti_obs <- 10^(obs[obs$STATION_KEY == 21792, 'log10_BSTI']/100*max_log10_bsti)
(1 - (bsti_target / bsti_obs)) * 100

# #Untransform bsti and keep track of which percentile the prediction was made with
# for (i in 1:length(df_list)) {
#   quant <- names(df_list[i])
#   df_list[[i]]$quant <- quant 
#   df_list[[i]]$bsti_untran <- 10^(df_list[[i]]$log10_BSTI/100 * max_log10_bsti)
# }
# 
# #Convert list to a dataframe to prepare for plotting
# df_quant <- ldply(df_list)
#Set the percentiles as a factor for plotting
df_quant$quant <- factor(df_quant$quant, levels = c('0%', '10%', '20%', '30%', 
                                                    '40%', '50%', '60%', '70%', 
                                                    '80%', '90%', '100%'))
save(df_quant, file = "tmdl_target_curve_data.Rdata")
to_plot <- df_quant[df_quant$STATION_KEY %in% ss$STATION_KEY,]
plot(to_plot$quant, ss$bsti_untran)

betahat <-dcast(data.frame(variable = rownames(fit$estimates$betahat), 
                           betahat = fit$estimates$betahat), 
                . ~ variable, 
                value.var = "betahat")[,-1]
ss$SITE_NAME <- as.character(ss$SITE_NAME)
for (i in 1:nrow(ss)) {
  to_plot <- df_quant[df_quant$STATION_KEY == ss[i, 'STATION_KEY'],]
  # g = ggplot(data = to_plot, aes(x = sum_1095_days, y = bsti_untran)) + 
  #         geom_point() + ggtitle(ss$STATION_KEY[i]) + ylim(0, 50) 
  #g = g + geom_smooth(method = 'loess') 
  #g = g + geom_line()
  rf_range <- range(dfdall[dfdall$STATION_KEY == ss[i, 'STATION_KEY'], 'sum_1095_days'])
  bm_list <- simplify_target_equation(betahat, ss, ss[i, 'STATION_KEY'])
  g = ggplot() + stat_function(data = data.frame(x = rf_range), 
                        aes(x, color = 'red'), 
                        fun = function(x) 10^(bm_list$b_u + bm_list$m_u * x))
  g = g + geom_point(data = to_plot, aes(x = sum_1095_days, y = bsti_untran)) + 
    ylim(0, 50) + ggtitle(paste(ss[i, c('STATION_KEY','SITE_NAME')], collapse = " - "))
  print(g)
  
  tmp_df_bm <- as.data.frame(bm_list)
  tmp_df_bm <- cbind(ss[i ,c('STATION_KEY','SITE_NAME')], tmp_df_bm)
  
  if (i == 1) {
    df_bm <- tmp_df_bm
  } else {
    df_bm <- rbind(df_bm, tmp_df_bm)
  }
  # df_ob <- data.frame('s' = ss$STATION_KEY[i],
  #            'o' = ssid[ssid$STATION_KEY == ss$STATION_KEY[i],'BSTI'],
  #            'p' = preds[preds$STATION_KEY %in% ss$STATION_KEY[i],'sum_1095_days'])
  # print(g + geom_point(data = df_ob, aes(x = p, y = o), color = "orange"))
}

#preds.0$BSTI_untran <- 10^(preds.0$log10_BSTI)
preds.0$uci <- preds.0$log10_BSTI + (critval * preds.0$log10_BSTI.predSE)
preds.0$lci <- preds.0$log10_BSTI - (critval * preds.0$log10_BSTI.predSE)


obs_sub$BSTI <- 10^(obs_sub$log10_BSTI/100 * max_log10_bsti)
preds.0 <- merge(preds.0, obs_sub[,c('STATION_KEY', 'BSTI')], by = 'STATION_KEY', all.x = TRUE)



#CART for Target identification
impaired$log10_Target <- log10(impaired$TMDL.FSS.Target)
impaired$pr_target <- abs(((impaired$TMDL.FSS.Target - impaired$FSS)/impaired$FSS) * 100)
impaired$pr_Q90_target <- abs(((impaired$Sediment.Stressor.Benchmark - impaired$FSS)/impaired$FSS) * 100)

impaired.0 <- merge(preds.0, impaired[,names(impaired) != 'FSS'], by.x = 'STATION_KEY', by.y = 'Station.Key')
impaired.0$log10_Target <- impaired.0$log10_Target / max_log10_bsti * 100
#impaired.0$target_met <- ifelse(impaired.0$log10_BSTI < (log10(impaired.0$Sediment.Stressor.Benchmark)/ max_log10_bsti * 100),1,0)
impaired.0$untran_BSTI <- round(10^(impaired.0$log10_BSTI/100 * max_log10_bsti))
impaired.0$pr_achieved <- abs(((impaired.0$untran_BSTI - impaired.0$BSTI)/impaired.0$BSTI) * 100)
impaired.0$pr_target_met <- ifelse(impaired.0$pr_achieved >= impaired.0$pr_target,1,0)
impaired.0$pr_Q90_target_met <- ifelse(impaired.0$pr_achieved >= impaired.0$pr_Q90_target,1,0)
impaired.0$untran_uci <- 10^(impaired.0$uci/100 * max_log10_bsti)
impaired.0$untran_lci <- 10^(impaired.0$lci/100 * max_log10_bsti)
impaired.0$lci_meets <- ifelse(impaired.0$lci < impaired.0$log10_Target, 1, 0)
#impaired.0[,c('STATION_KEY','log10_BSTI','uci','lci','log10_Target','untran_BSTI','TMDL.BSTI.Target','target_met',all.vars(ssn1_glmssn5$args$formula))]
impaired.0 <- merge(impaired.0, preds_obs, by = 'pid', all.x = TRUE)
impaired.0$BSTI_obs <- 10^(obs_sub$log10_BSTI_obs/100 * max_log10_bsti)

df <- melt(impaired.0[,c('STATION_KEY','BSTI','untran_BSTI',
                         'untran_uci','untran_lci')], id.vars = "STATION_KEY",
     measure.vars = c('untran_BSTI', 'untran_uci', 'untran_lci', 'BSTI'))

# df2 <- melt(ss, id.vars = "STATION_KEY", measure.vars = c('BSTI', 'BSTI_target',
#                                                          'untran_uci', 'untran_lci'))

hline.data <- data.frame(z = impaired.0[,"TMDL.FSS.Target"], 
                         STATION_KEY = impaired.0[,c('STATION_KEY')])
ggplot(df, 
       aes(x = 1, y = value, color = variable)) + 
  geom_point(size = c(3, 8, 8, 3), shape = c(19, 95, 95, 19)) + 
  scale_color_discrete(labels = c("Predicted BSTI at 0 Anthro", "UCI", "LCI", "Observed BSTI")) +
  facet_wrap(~ STATION_KEY) +
  geom_hline(aes(yintercept = z), hline.data) 

# impaired.0[,c('STATION_KEY',grep('untran',names(impaired.0), value = TRUE),
#               'TMDL.BSTI.Target','Sediment.Stressor.Benchmark','target_met')]
