library(SSN)

options(stringsAsFactors = FALSE)

obs_df <- getSSNdata.frame(ssn1)

locs <- levels(obs_df$locID)

#Number of sites with number of samples
table(table(obs_df$locID))
#   1   2   3   4   5   6   7   8  10  11  13 
# 448  53  27  11   4   2   1   1   1   1   1 

nsites <- 5

for (i in 1:(length(locs))) {
  #Grab the individual station to predict on
  stn <- locs[i]
  
  #Grab the fresh data frame
  obs_NA <- obs_df
  
  #Set the observed values at this station to NA
  obs_NA[obs_NA$locID == stn, 'log10_BSTI'] <- NA
  
  #Put it in the object
  ssn_tmp <- putSSNdata.frame(obs_NA, ssn1)
  
  #Fit with the missing values
  tmp_fit <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + KFACT_MARCA +
                      OWN_FED_PRCA + DIS_3YR_PRSA + ROADLEN_DRSA + OWN_URB_PARCA + HDWTR,
                    EstMeth = "REML",
                    ssn_tmp,
                    CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
                    addfunccol = "afvArea",
                    family = "Gaussian")
  
  #Predict on the missing values TempFitMissingObs = tfmo
  tfmo <- predict(tmp_fit, "_MissingObs_")
  
  getPreds(tfmo, pred.type = "pred")
  
  with(getSSNdata.frame(ssn1), log10_BSTI[locID == stn])
}


