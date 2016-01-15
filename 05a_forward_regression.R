#### SubSSN Test ####
library(SSN)
library(stringr)

options('scipen' = 100)

##################################################
### Create the Distance Matrix
###################################################
#createDistMat(ssn1, predpts = 'preds_up', o.write = TRUE, amongpreds = TRUE)

###################################################
### Variable Selection Via Backward Deletion   ###
###################################################
vars <- pkeep

for (i in 1:length(vars)) {
  start.time <- Sys.time()
  print(paste("Starting model fit", i, "at", start.time))
  
  tmp <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14',vars[1:i])]),
                EstMeth = "ML",
                ssn1,
                CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
                addfunccol = "afvArea",
                family = "Gaussian")
  save_name <- paste0("ssn1_glmssn",i,'_forward_HWFAC_ML_20151230.Rdata')
  save(tmp, file = save_name)
  
  end.time <- Sys.time()
  print(paste("Completed model fit", i, "and summary at", end.time))
  print(paste("elapsed time was", end.time - start.time))
  cat("\n\n")
  
  rm(tmp)
}
