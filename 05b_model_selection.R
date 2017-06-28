library(SSN)

source('functions_custom.R')

name_vec <- paste0("ssn1_glmssn_std_", 1:7, "_ML.Rdata")
#name_vec <- c(name_vec,paste0("ssn1_glmssn", 1:8, "_RUN2_ML.Rdata"))
#name_vec <- c(name_vec,paste0("ssn1_glmssn",1:8,"_ML_RUN3.Rdata"))
#name_vec <- c(name_vec, paste0("ssn1_glmssn", 1:9, "_HWFAC_ML_20151216.Rdata"))
for (i in 1:length(name_vec)) {
  load(name_vec[i])
  # assign(paste0(gsub("[0-9]+_ML.Rdata|{0-9}._RUN2_ML.Rdata",
  #                    "",name_vec[i]),i), tmp)
}
c(paste("ssn1_glmssn_std_RUN1_", 1:7, sep = ""), 
  paste("ssn1_glmssn_std_RUN2_", 1:7, sep = ""))
models <- list(ssn1_glmssn_std_RUN1_1, ssn1_glmssn_std_RUN1_2, ssn1_glmssn_std_RUN1_3, ssn1_glmssn_std_RUN1_4,
               ssn1_glmssn_std_RUN1_5, ssn1_glmssn_std_RUN1_6, ssn1_glmssn_std_RUN1_7, ssn1_glmssn_std_RUN2_1,
               ssn1_glmssn_std_RUN2_2, ssn1_glmssn_std_RUN2_3, ssn1_glmssn_std_RUN2_4, ssn1_glmssn_std_RUN2_5,
               ssn1_glmssn_std_RUN2_6, ssn1_glmssn_std_RUN2_7)
bhats <- lapply(models, function(x) {x$estimates$betahat})
results <- InfoCritCompare2(models)
results$dAIC <- min(results$AIC) - results$AIC
results$RMSPE_untran <- 10^results$RMSPE
results[order(results$RMSPE, decreasing = TRUE),]

#The last model from both runs have similar support 

log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + OWN_FED_PRCA + DIS_1YR_PARSA + HDWTR

log10_BSTI ~ STRMPWR + PPT_1981_2010 + MIN_Z + KFACT_MARCA + DIS_1YR_PARSA + OWN_FED_PRCA + TYPEF_PARCA + SQM_ARSA