library(SSN)

options(stringsAsFactors = FALSE)

load('back_for_ward_results_20151019.Rdata')

#Extract the model variables we want to use in re-fitting the model
vars_to_use <- all.vars(as.formula(as.character(results$formula[15])))

obs <- getSSNdata.frame(ssn1, Name = "Obs")

#obs_sub <- obs[,c('STATION_KEY',vars_to_use)]
obs[obs$STRMPWR == 0, 'STRMPWR'] <- 0.00001
obs <- obs[obs$STRMPWR != 0.00001,]
obs$genius <- 1/(obs$MIN_Z * obs$STRMPWR)
obs$EROD_genius <- obs$EROD_PARCA * obs$genius
obs$DIS_genius <- obs$DIS_1YR_PARSA * obs$genius
obs$AGR_genius <- obs$OWN_AGR_PARCA * obs$genius
obs$PRI_genius <- obs$OWN_PRI_PRCA * obs$genius

cov(obs[,all.vars(as.formula(as.character(results$formula[15])))[2:7]])
cov(obs[,c('MIN_Z','STRMPWR','DIS_genius','EROD_genius','AGR_genius','PRI_genius')])

ssn1 <- putSSNdata.frame(DataFrame = obs, x = ssn1, Name = "Obs")

ssn1_glmssn5a <- glmssn(formula = log10_FSS_26Aug14 ~ MIN_Z + STRMPWR + 
                          DIS_genius + EROD_genius + AGR_genius + PRI_genius,
                        EstMeth = "REML",
                        ssn1,
                        CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'), 
                        addfunccol = "afvArea")
save(ssn1_glmssn5, file = 'ssn1_glmssn5a_11042015.Rdata')