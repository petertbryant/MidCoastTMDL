library(SSN)
library(dplyr)

options(stringsAsFactors = FALSE)

ref <- read.csv('reference_sites_03102017.csv')

#To adequately use and interpet the model there are some precursors that are needed
#so to that end we need to re-run some of the code with the dataset used to fit the
#selected model
ssn_RF_data <- read.csv("ssn_RF_data.csv")

bsti <- ssn_RF_data[, !colnames(ssn_RF_data) %in% c('SVN','STATION_KEY',
                                                    'rid','rid_LSN04')]

bsti$DATE <- as.POSIXct(bsti$DATE)

#This removes those variables where all the values are 0
bsti <- bsti[, setdiff(names(bsti), c("X2year_count_60_days",
                                      "X10year_count_60_days", 
                                      "X25year_count_60_days", 
                                      "X50year_count_60_days",
                                      "X100year_count_60_days", 
                                      "X10year_count_180_days", 
                                      "X25year_count_180_days",
                                      "X50year_count_180_days", 
                                      "X100year_count_180_days"))]
#Need to convert Inf values to 0
bsti[, grep("X", names(bsti))] <- as.data.frame(sapply(
  bsti[, grep("X",names(bsti))], function(x) {
    replace(x, is.infinite(x),0)
  }
))

#Soil characteristics are known complements of each other. We are going
#to select the size class representative of fine sediment <0.05mm
bsti <- bsti[, -grep('^SAND|^CLAY|^SILT_P', names(bsti))]

#COMP is a complement of EROD and has high correlation. Although EROD is 
#a component of the derivation of SUSCEP it maintains low correlation < 0.4
bsti <- bsti[, -grep('COMP', names(bsti))]

bsti.run <- within(bsti, rm(STRMPWR))

stdpreds <- function(newset,originalset) {
  xnames <- colnames(newset)
  sx <- matrix(rep(NA,ncol(newset)*nrow(newset)),nrow=nrow(newset))
  mx <- c(rep(NA,ncol(newset)))
  sdx <- c(rep(NA,ncol(newset)))
  for(i in 1:ncol(newset)) {
    var <- with(originalset,get(xnames[i]))
    sx[,i] <- (newset[,i]-mean(var))/(2*sd(var))
    mx[i] <- mean(var)
    sdx[i] <- sd(var)
  }
  colnames(sx) <- colnames(newset)
  names(mx) <- colnames(newset)
  names(sdx) <- colnames(newset)
  attr(sx, "mean") <- mx
  attr(sx, "sd") <- sdx
  return(sx)
}

bsti_matrix <- stdpreds(bsti.run, bsti.run)
var_means <- attr(bsti_matrix, "mean")
var_sd <- attr(bsti_matrix, "sd")
bsti.run <- (as.data.frame(bsti_matrix))


ssn1 <- importSSN("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/LSN05_Watersheds/LSN07_IndianCreek/lsn.ssn",
                  predpts = "preds", o.write = TRUE)
obs<- getSSNdata.frame(ssn1, Name = "Obs")
ssn_RF_data <- ssn_RF_data[,c("STATION_KEY", "SVN", "DATE","YEAR",'BSTI')]
ssn_RF_data <- cbind(ssn_RF_data,
                     bsti.run[,c('sum_1095_days','XSLOPE_MAP',
                                 "MIN_Z",'OWN_FED_PRCA','DIS_1YR_PARSA','HDWTR')])
ssn_RF_data$log10_BSTI <- log10(ssn_RF_data$BSTI)
obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",
                         "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME",
                         "HU_12_NAME", "HU_08", "NHDHigh", "LONG_RAW",
                         "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10",
                         "HU_12")],
                  ssn_RF_data,
                  by = "SVN",
                  all.x = TRUE)
obs.vars$HDWTR <- as.factor(obs.vars$HDWTR)
obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
row.names(obs.vars) <- obs.vars$pid
levels(obs.vars$HDWTR) <- c(0,1)
ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

preds <- getSSNdata.frame(ssn1, Name = "preds")
pid.order <- preds$pid
preds <- rename(preds, sum_1095_days = sum1095, OWN_FED_PRCA = own_fed_rc, DIS_1YR_PARSA = dis_1yr_20)
preds[,c('MIN_Z','XSLOPE_MAP', 'DIS_1YR_PARSA','OWN_FED_PRCA','sum_1095_days')] <- 
  stdpreds(newset = preds[,c('MIN_Z','XSLOPE_MAP', 'DIS_1YR_PARSA','OWN_FED_PRCA','sum_1095_days')],
         originalset = ssn_RF_data[,c('MIN_Z','XSLOPE_MAP', 'DIS_1YR_PARSA','OWN_FED_PRCA','sum_1095_days')])
preds$HDWTR <- as.factor(preds$HDWTR)
levels(preds$HDWTR) <- c(0,1)
preds <- preds[match(pid.order,preds$pid),]
row.names(preds) <- preds$pid
# ssn1 <- putSSNdata.frame(preds.vars, ssn1, Name = "obs_prd")
ssn1 <- putSSNdata.frame(preds, ssn1, Name = "preds")
createDistMat(ssn1, o.write = TRUE, predpts = "preds", amongpreds = TRUE)

rm(bsti, bsti_matrix, ssn_RF_data, bsti.run, obs, obs.vars, preds)


fit <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + OWN_FED_PRCA +
                DIS_1YR_PARSA + HDWTR,
              EstMeth = "REML",
              ssn1,
              CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
              addfunccol = "afvArea",
              family = "Gaussian")

predictions <- predict(fit_ic,predpointsID = "preds")

pred_df <- getSSNdata.frame(predictions, "preds")
as.integer(round(10^pred_df$log10_BSTI))
as.integer(round(10^pred_df$log10_BSTI.predSE))

ref_DIS_1YR_PARSA <- quantile(
  obs[obs$STATION_KEY %in% ref$STATION_KEY,'DIS_1YR_PARSA'], 
  seq(0,1,.25))[4]
ref_OWN_FED_PRCA <- quantile(
  obs[obs$STATION_KEY %in% ref$STATION_KEY,'OWN_FED_PRCA'], 
  seq(0,1,.25))[2]

target_df <- preds

target_df[,'DIS_1YR_PARSA'] <- ref_DIS_1YR_PARSA
target_df[target_df$OWN_FED_PRCA < ref_OWN_FED_PRCA,
          'OWN_FED_PRCA'] <- ref_OWN_FED_PRCA

fit_ic_0 <- putSSNdata.frame(target_df, fit_ic, Name = "preds")
targets <- predict(fit_ic_0, predpointsID = "preds")
targets_df <- getSSNdata.frame(targets, "preds")
as.integer(round(10^targets_df$log10_BSTI))
as.integer(round(10^targets_df$log10_BSTI.predSE))
targets_df <- rename(targets_df, target = log10_BSTI, target.SE = log10_BSTI.predSE)

pred_df <- merge(pred_df, targets_df[,c('pid','target','target.SE')], by = 'pid')
pred_df$BSTI <- as.integer(round(10^pred_df$log10_BSTI))
pred_df$BSTI.SE <- as.integer(round(10^pred_df$log10_BSTI.predSE))
pred_df$target <- as.integer(round(10^pred_df$target))
pred_df$target.SE <- as.integer(round(10^pred_df$target.SE))
pred_df$diff <- pred_df$BSTI - pred_df$target
View(pred_df[,c('BSTI','BSTI.SE','target','target.SE','diff')])
write.csv(pred_df[,c('pid','BSTI','BSTI.SE','target','target.SE','diff')], file = "IC_predictions_targets.csv")
