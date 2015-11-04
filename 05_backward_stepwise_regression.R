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
halt = FALSE
cntr = 1
vars <- names(obs.fss2)

while (halt == FALSE) {
  start.time <- Sys.time()
  print(paste("Starting model fit", cntr, "at", start.time))
  
  tmp <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14',vars)]),
                EstMeth = "REML",
                ssn1,
                CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
                addfunccol = "afvArea",
                family = "Gaussian")
  save_name <- paste0("ssn1_glmssn",cntr,'_20151019.Rdata')
  save(tmp, file = save_name)
  
  sum_tab <- summary(tmp)$fixed.effects.estimates
  var_to_remove <- sum_tab[which.max(sum_tab$prob.t),'FactorLevel']
  
  end.time <- Sys.time()
  print(paste("Completed model fit", cntr, "and summary at", end.time))
  print(paste("elapsed time was", end.time - start.time))
  cat("\n\n")
  
  if (max(sum_tab$prob.t) < 0.05) {halt = TRUE}
  
  rm(tmp)
  vars <- vars[!vars %in% var_to_remove]
  cntr = cntr + 1
}

ssn1_glmssn5 <- glmssn(log10_FSS_26Aug14 ~ MIN_Z + STRMPWR + DIS_1YR_PARSA +
                         EROD_PARCA + OWN_AGR_PARCA + OWN_PRI_PRCA,
                       EstMeth = "REML",
                       ssn1,
                       CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'), 
                       addfunccol = "afvArea")
save(ssn1_glmssn5, file = 'ssn1_glmssn5_20151019.Rdata')
