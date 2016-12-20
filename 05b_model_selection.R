library(SSN)

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
results[order(results$RMSPE, decreasing = TRUE),]

save(results, file = "back_results_20160715.Rdata")
load('back_results_20160715.Rdata')

# results[1:13, 'method'] <- 'forward'
# results[14:22, 'method'] <- 'backward'

# load('back_for_ward_results_SLOPEQ_20151106.Rdata')
# results$track <- paste("SLOPEQ",c(paste0("f",c(1:5,7:13)), paste0("b",1:8)))
# results.slopeq <- results
# load('back_for_ward_results_STRMPWR_20151106.Rdata')
# results$track <- paste("STRMPWR",c(paste0("f",c(1:13)), paste0("b",1:6)))
# results_all <- rbind(results, results.slopeq)
# results_all <- results_all[!duplicated(results_all$AIC),]
# results_all$dAIC <- min(results_all$AIC) - results_all$AIC
# results_all[order(results_all$dAIC, decreasing = TRUE),]
# 
# 
# fit1 <- glmssn(log10_FSS_26Aug14 ~ sum_1095_days + XSLOPE_MAP + MIN_Z + KFACT_MARCA + OWN_PRI_PRCA + SUSCEP4_PRCA + POP_DARCA,
#                EstMeth = "ML",
#                ssn1,
#                CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
#                addfunccol = "afvArea",
#                family = "Gaussian")
# results_fit1 <- InfoCritCompare(list(fit1))
# results_fit1$method <- "extra"
# results_fit1$dAIC <- NA
# results_plus <- rbind(results, results_fit1)
# results_plus$dAIC <- min(results_plus$AIC) - results_plus$AIC
# results_plus[order(results_plus$dAIC, decreasing = TRUE),]
# 
# fit2 <- glmssn(log10_FSS_26Aug14 ~ sum_1095_days + XSLOPE_MAP + MIN_Z + OWN_PRI_PRCA + POP_DARCA,
#                EstMeth = "ML",
#                ssn1,
#                CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
#                addfunccol = "afvArea",
#                family = "Gaussian")
# results_fit2 <- InfoCritCompare(list(fit2))
# results_fit2$method <- "extra"
# results_fit2$dAIC <- NA
# results_plus <- rbind(results_plus, results_fit2)
# results_plus$dAIC <- min(results_plus$AIC) - results_plus$AIC
# results_plus[order(results_plus$dAIC, decreasing = TRUE),]
# 
# fit2 <-  glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14','sum_1095_days',
#                                        'MIN_Z','SUSCEP5_PARCA','KFACT_MARCA',
#                                        'OWN_PRI_PRCA','ROADLEN_DARSA',
#                                        'OWN_URB_PARCA', 'DIS_1YR_PARSA')]),
#                 EstMeth = "REML",
#                 ssn1,
#                 CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
#                 addfunccol = "afvArea",
#                 family = "Gaussian")
# 
# results_sub <- InfoCritCompare(list(ssn1_glmssn5, ssn1_glmssn6, fit1, fit2))
# results_sub$dAIC <- min(results_sub$AIC) - results_sub$AIC
# results_sub[order(results_sub$dAIC, decreasing = TRUE),]
# save(results_sub, file = 'manual_stepwise_results_201511106.Rdata')
# 
# load('manual_stepwise_results_201511106.Rdata')
# results_sub$track <- paste("manual", paste0("m",1:4))
# results_all <- rbind(results_all, results_sub)
# results_all <- results_all[!duplicated(results_all$formula),]
# results_all$dAIC <- min(results_all$AIC) - results_all$AIC
# results_all[order(results_all$dAIC, decreasing = TRUE),]
# 
# #4 overlapping AIC, all SLOPEQ based. Lowest RMSPE is manual m4
# #however, SLOPEQ b5 seems to be similar in RMSPE. Let's do scenarios with this one
# #first and then look at the others.
# load('ssn1_glmssn5_SLOPEQ_20151106.Rdata')
# assign(x = "ssn1_glmssn5", value = tmp)
# rm(tmp)
