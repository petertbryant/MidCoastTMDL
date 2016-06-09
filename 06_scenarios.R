library(SSN)
library(plyr)
library(RODBC)
library(ggplot2)
library(reshape2)

options(stringsAsFactors = FALSE)

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
preds_obs <- getSSNdata.frame(fit, Name = "preds")

#Takes the highest BSTI score from each station
obs_sub <- obs[,c('STATION_KEY',all.vars(fit$args$formula))]
obs_sub <- obs_sub[order(obs_sub$STATION_KEY, obs_sub$log10_BSTI, decreasing = TRUE),]
obs_sub <- obs_sub[!duplicated(obs_sub$STATION_KEY),]

#Preserve observed BSTI values in the prediction data set
preds_obs <- preds_obs[,c('pid','log10_BSTI')]
preds_obs <- rename(preds_obs, c('log10_BSTI' = 'log10_BSTI_obs'))


#### ####
#### ####
#### ####
#### ####
#### ####

source('C:/Users/pbryant/Desktop/MidCoastTMDL/fun_simplify_target_equation.R')

betahat <-dcast(data.frame(variable = rownames(fit$estimates$betahat), 
                           betahat = fit$estimates$betahat), 
                . ~ variable, 
                value.var = "betahat")[,-1]

#Test to see if changing each one individually affects prediction estimation
#RESULT: No difference in predicitions when each all are modified at the same time
#Generate predictions at TMDL Target conditions and at aobserved rainfall amounts
preds.0 <- getSSNdata.frame(fit, Name = "preds")
preds.0[, 'ROADLEN_DRSA'] <- quantile(
          preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'ROADLEN_DRSA'], 
          seq(0,1,.25))[4]
preds.0[, 'OWN_URB_PARCA'] <- quantile(
          preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'OWN_URB_PARCA'], 
          seq(0,1,.25))[2]

ref_DIS_3YR_PRSA <- quantile(
  preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'DIS_3YR_PRSA'], 
  seq(0,1,.25))[4]
preds.0[preds.0$DIS_3YR_PRSA > ref_DIS_3YR_PRSA, 
        'DIS_3YR_PRSA'] <- ref_DIS_3YR_PRSA

ref_OWN_FED_PRCA <- quantile(
  preds.0[preds.0$STATION_KEY %in% ref$STATION_KEY,'OWN_FED_PRCA'], 
  seq(0,1,.25))[2]
preds.0[preds.0$OWN_FED_PRCA < ref_OWN_FED_PRCA, 
        'OWN_FED_PRCA'] <- ref_OWN_FED_PRCA

#Run the fit
fit_0 <- putSSNdata.frame(preds.0, fit, Name = "preds")

#Run the prediction
fit_0_preds <- predict.glmssn(fit_0, predpointsID = "preds", 
                              newdata = 'preds')

#Check the results
ssid_all <- getSSNdata.frame(fit_0_preds, Name = 'preds')

for (i in 1:nrow(ssid_all)) {
  bm_list <- simplify_target_equation(betahat, ssid_all, ssid_all[i, 'STATION_KEY'])
  tmp_df_bm <- as.data.frame(bm_list)
  tmp_df_bm <- cbind(ssid_all[i ,c('STATION_KEY','SITE_NAME')], tmp_df_bm)

  if (i == 1) {
    df_bm <- tmp_df_bm
  } else {
    df_bm <- rbind(df_bm, tmp_df_bm)
  }
}
ssid <- ssid_all
ssid <- merge(ssid, preds_obs, by = 'pid')
ssid$BSTI <- as.integer(10^(ssid$log10_BSTI_obs/100 * max_log10_bsti))
ssid$BSTI_target <- as.integer(round(10^(ssid$log10_BSTI/100 * max_log10_bsti)))
ssid$Sed_Stressor <- ifelse(ssid$BSTI > ssid$BSTI_target,TRUE,FALSE)
ssid$pr_target <- round(abs(((ssid$BSTI - ssid$BSTI_target)/ssid$BSTI) * 100),1)
#ssid$predSE_untran <- 10^(ssid$log10_BSTI.predSE/100 * max_log10_bsti)
ssid$STATION_KEY <- as.character(ssid$STATION_KEY)
ssid <- ssid[order(ssid$STATION_KEY, decreasing = TRUE),]
ss <- ssid[ssid$STATION_KEY %in% impaired$STATION_KEY & ssid$Sed_Stressor,]

#CART vs SSN condition compare
ss2 <- ssid[ssid$STATION_KEY %in% impaired$STATION_KEY,]
cc <- merge(impaired, ss2, by = 'STATION_KEY')
cc_sub <- cc[,c('STATION_KEY','SITE_NAME.x','FSS','Q90TH','sediment_resid_status','BSTI','BSTI_target','Sed_Stressor')]
cc_sub$CART_Sed_Stressor <- ifelse(cc_sub$BSTI > cc_sub$Q90TH,TRUE,FALSE)
cc_sub$agree <- ifelse(cc_sub$Sed_Stressor & cc_sub$CART_Sed_Stressor, TRUE, 
                       ifelse(!cc_sub$Sed_Stressor & !cc_sub$CART_Sed_Stressor, TRUE, FALSE))
cc_sub <- within(cc_sub, rm(FSS, sediment_resid_status))
cc_sub <- cc_sub[!duplicated(cc_sub$STATION_KEY),]
cc_sub <- plyr::rename(cc_sub, c('SITE_NAME.x' = 'SITE_NAME', 'Q90TH'= 'CART_Target', 'BSTI_target' = 'SSN_Target', 'Sed_Stressor' = 'SSN_Sed_Stressor'))
cc_sub <- cc_sub[,c('STATION_KEY','SITE_NAME','BSTI','CART_Target','CART_Sed_Stressor','SSN_Target','SSN_Sed_Stressor','agree')]
write.csv(cc_sub, file = 'ssn_cart_compare.csv', row.names = FALSE)

df_bm_ss <- df_bm[df_bm$STATION_KEY %in% ss$STATION_KEY,]
write.csv(df_bm_ss, file = 'b_values_sediment_stressor_sites.csv', row.names = FALSE)

# lapply(list('sum_1095_days', 'XSLOPE_MAP', 'MIN_Z', 'KFACT_MARCA',
#              'ROADLEN_DRSA', 'HDWTR'), function(x) {preds[,x] <- (preds[,x]/100) * 
#                min.max[min.max$variable == x, 'max_val']})
# 
# unique(ss$HU_10_NAME)
# 
# CART_HUCs <- preds[preds$STATION_KEY %in% CART_imp$Station.Key, 
#                    c('STATION_KEY',grep("HU",names(preds),value = TRUE))]
# unique(ss$HU_12_NAME)[!unique(ss$HU_12_NAME) %in% unique(CART_HUCs$HU_12_NAME)]

#### Generate target range ####
load('C:/users/pbryant/desktop/midcoasttmdl-gis/precip_daily_sum_1095_days.Rdata')
dfdall$sum_1095_days <- dfdall$sum_1095_days / 
  (min.max[min.max$variable == 'sum_1095_days', 'max_val'])*100


#Percent reduction 
rf <- obs[obs$STATION_KEY == 21792, 'sum_1095_days']
bsti_target <- (10^(df_bm[31,3] + df_bm[31,4] * rf))
bsti_obs <- 10^(obs[obs$STATION_KEY == 21792, 'log10_BSTI']/100*max_log10_bsti)
(1 - (bsti_target / bsti_obs)) * 100

ss$SITE_NAME <- as.character(ss$SITE_NAME)
for (i in 1:nrow(ss)) {
  #Pull out the rainfall range for the ith site
  rf_range <- range(dfdall[dfdall$STATION_KEY == ss[i, 'STATION_KEY'], 'sum_1095_days'])
  rf_range <- min.max[min.max$variable == 'sum_1095_days', "max_val"]*(rf_range/100)
  
  #Pull out the simplified equation values to use for plotting
  bm_list <- simplify_target_equation(betahat, ss, ss[i, 'STATION_KEY'])
  
  #Build plot using the simplified equation
  g = ggplot() + stat_function(data = data.frame(x = rf_range), size = 1, 
                        aes(x), 
                        fun = function(x) 10^(bm_list$b_u + (-5.228741e-05 * x)))
  
  #Set formatting on plot
  g = g + ylim(0, 50) + 
    ggtitle(paste(ss[i, c('STATION_KEY','SITE_NAME')], collapse = " - ")) +
    xlab("3 year sum of rainfall (mm)") + ylab("BSTI")
  
  #Extract the observed BSTI and the rainfall at which it occurs
  df_ob <- data.frame('s' = ss$STATION_KEY[i],
             'o' = ssid[ssid$STATION_KEY == ss$STATION_KEY[i],'BSTI'],
             'p' = min.max[min.max$variable == 'sum_1095_days', "max_val"]*(preds[preds$STATION_KEY %in% ss$STATION_KEY[i],'sum_1095_days'])/100)
  
  #Add the observed BSTI to the plot
  print(g + geom_point(data = df_ob, aes(x = p, y = o), color = "orange") + theme(legend.position = "none"))
}

#Compare SSN sed stressor ID to CART Stressor ID 
#CART sites not in SSN
CART_imp[!CART_imp$Station.Key %in% ss$STATION_KEY, c('Station.Key','Site.Name')]
#    Station.Key                   Site.Name
# 10       34695        Fivemile Cr at Mouth
# 16    dfw_2464                Pittenger Cr
# 24   dfw_36277                   Sweet Cr.
# 26    dfw_2492 Wolf Cr at RM 3.12 (Umpqua)

#SSN sites not in CART
ss[!ss$STATION_KEY %in% CART_imp$Station.Key, c('STATION_KEY', 'SITE_NAME')]
#     STATION_KEY                 SITE_NAME
# 171   dfw_39988                    Haight
# 57    dfw_39983                 Drift Cr.
# 352       37188 Big Elk Creek at RM 12.77
# 157       37185   Morris Creek at RM 3.82
# 141       37165    McLeod Creek at RM 2.1
# 365       34659   Beaver Cr NF at RM 3.06
# 370       34521             Needle Branch
# 165       33331                  Jeans Cr
# 410       33325                   Wolf Cr

#Number of stations in agreement
nrow(ss[ss$STATION_KEY %in% CART_imp$Station.Key,])
# 22

#preds.0$BSTI_untran <- 10^(preds.0$log10_BSTI)
# preds.0$uci <- preds.0$log10_BSTI + (critval * preds.0$log10_BSTI.predSE)
# preds.0$lci <- preds.0$log10_BSTI - (critval * preds.0$log10_BSTI.predSE)
# 
# 
# obs_sub$BSTI <- 10^(obs_sub$log10_BSTI/100 * max_log10_bsti)
# preds.0 <- merge(preds.0, obs_sub[,c('STATION_KEY', 'BSTI')], by = 'STATION_KEY', all.x = TRUE)
# 
# 
# 
# #CART for Target identification
# impaired$log10_Target <- log10(impaired$TMDL.FSS.Target)
# impaired$pr_target <- abs(((impaired$TMDL.FSS.Target - impaired$FSS)/impaired$FSS) * 100)
# impaired$pr_Q90_target <- abs(((impaired$Sediment.Stressor.Benchmark - impaired$FSS)/impaired$FSS) * 100)
# 
# impaired.0 <- merge(preds.0, impaired[,names(impaired) != 'FSS'], by.x = 'STATION_KEY', by.y = 'Station.Key')
# impaired.0$log10_Target <- impaired.0$log10_Target / max_log10_bsti * 100
# #impaired.0$target_met <- ifelse(impaired.0$log10_BSTI < (log10(impaired.0$Sediment.Stressor.Benchmark)/ max_log10_bsti * 100),1,0)
# impaired.0$untran_BSTI <- round(10^(impaired.0$log10_BSTI/100 * max_log10_bsti))
# impaired.0$pr_achieved <- abs(((impaired.0$untran_BSTI - impaired.0$BSTI)/impaired.0$BSTI) * 100)
# impaired.0$pr_target_met <- ifelse(impaired.0$pr_achieved >= impaired.0$pr_target,1,0)
# impaired.0$pr_Q90_target_met <- ifelse(impaired.0$pr_achieved >= impaired.0$pr_Q90_target,1,0)
# impaired.0$untran_uci <- 10^(impaired.0$uci/100 * max_log10_bsti)
# impaired.0$untran_lci <- 10^(impaired.0$lci/100 * max_log10_bsti)
# impaired.0$lci_meets <- ifelse(impaired.0$lci < impaired.0$log10_Target, 1, 0)
# #impaired.0[,c('STATION_KEY','log10_BSTI','uci','lci','log10_Target','untran_BSTI','TMDL.BSTI.Target','target_met',all.vars(ssn1_glmssn5$args$formula))]
# impaired.0 <- merge(impaired.0, preds_obs, by = 'pid', all.x = TRUE)
# impaired.0$BSTI_obs <- 10^(obs_sub$log10_BSTI_obs/100 * max_log10_bsti)
# 
# df <- melt(impaired.0[,c('STATION_KEY','BSTI','untran_BSTI',
#                          'untran_uci','untran_lci')], id.vars = "STATION_KEY",
#      measure.vars = c('untran_BSTI', 'untran_uci', 'untran_lci', 'BSTI'))
# 
# # df2 <- melt(ss, id.vars = "STATION_KEY", measure.vars = c('BSTI', 'BSTI_target',
# #                                                          'untran_uci', 'untran_lci'))
# 
# hline.data <- data.frame(z = impaired.0[,"TMDL.FSS.Target"], 
#                          STATION_KEY = impaired.0[,c('STATION_KEY')])
# ggplot(df, 
#        aes(x = 1, y = value, color = variable)) + 
#   geom_point(size = c(3, 8, 8, 3), shape = c(19, 95, 95, 19)) + 
#   scale_color_discrete(labels = c("Predicted BSTI at 0 Anthro", "UCI", "LCI", "Observed BSTI")) +
#   facet_wrap(~ STATION_KEY) +
#   geom_hline(aes(yintercept = z), hline.data) 
# 
# impaired.0[,c('STATION_KEY',grep('untran',names(impaired.0), value = TRUE),
#               'TMDL.BSTI.Target','Sediment.Stressor.Benchmark','target_met')]
