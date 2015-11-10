library(SSN)

name_vec <- paste0("ssn1_glmssn_forward", 1:13, "_SLOPEQ_20151106.Rdata")
name_vec <- c(name_vec, paste0("ssn1_glmssn", 1:8, "_SLOPEQ_20151106.Rdata"))
for (i in 1:length(name_vec)) {
  load(name_vec[i])
  assign(gsub("_SLOPEQ_20151106.Rdata","",name_vec[i]), tmp)
}
#load("ssn1_glmssn5_20151019.Rdata")
results <- InfoCritCompare(list(ssn1_glmssn_forward1, ssn1_glmssn_forward2, 
                                ssn1_glmssn_forward3, ssn1_glmssn_forward4, 
                                ssn1_glmssn_forward5, #ssn1_glmssn_forward6,
                                ssn1_glmssn_forward7, ssn1_glmssn_forward8,
                                ssn1_glmssn_forward9, ssn1_glmssn_forward10,
                                ssn1_glmssn_forward11, ssn1_glmssn_forward12,
                                ssn1_glmssn_forward13,
                                ssn1_glmssn1, ssn1_glmssn2, 
                                ssn1_glmssn3, ssn1_glmssn4,
                                ssn1_glmssn5, ssn1_glmssn6,
                                ssn1_glmssn7, ssn1_glmssn8))

results$dAIC <- min(results$AIC) - results$AIC
results[order(results$dAIC, decreasing = TRUE),]
save(results,file = "back_for_ward_results_SLOPEQ_20151106.Rdata")

load('back_for_ward_results_SLOPEQ_20151106.Rdata')
results$track <- paste("SLOPEQ",c(paste0("f",c(1:5,7:13)), paste0("b",1:8)))
results.slopeq <- results
load('back_for_ward_results_STRMPWR_20151106.Rdata')
results$track <- paste("STRMPWR",c(paste0("f",c(1:13)), paste0("b",1:6)))
results_all <- rbind(results, results.slopeq)
results_all <- results_all[!duplicated(results_all$AIC),]
results_all$dAIC <- min(results_all$AIC) - results_all$AIC
results_all[order(results_all$dAIC, decreasing = TRUE),]


fit1 <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14','sum_1095_days',
                                      'MIN_Z','SUSCEP5_PARCA','KFACT_MARCA',
                                      'OWN_PRI_PRCA','ROADLEN_DARSA',
                                      'OWN_URB_PARCA')]),
               EstMeth = "REML",
               ssn1,
               CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
               addfunccol = "afvArea",
               family = "Gaussian")

fit2 <-  glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14','sum_1095_days',
                                       'MIN_Z','SUSCEP5_PARCA','KFACT_MARCA',
                                       'OWN_PRI_PRCA','ROADLEN_DARSA',
                                       'OWN_URB_PARCA', 'DIS_1YR_PARSA')]),
                EstMeth = "REML",
                ssn1,
                CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
                addfunccol = "afvArea",
                family = "Gaussian")

results_sub <- InfoCritCompare(list(ssn1_glmssn5, ssn1_glmssn6, fit1, fit2))
results_sub$dAIC <- min(results_sub$AIC) - results_sub$AIC
results_sub[order(results_sub$dAIC, decreasing = TRUE),]
save(results_sub, file = 'manual_stepwise_results_201511106.Rdata')

load('manual_stepwise_results_201511106.Rdata')
results_sub$track <- paste("manual", paste0("m",1:4))
results_all <- rbind(results_all, results_sub)
results_all <- results_all[!duplicated(results_all$formula),]
results_all$dAIC <- min(results_all$AIC) - results_all$AIC
results_all[order(results_all$dAIC, decreasing = TRUE),]

#4 overlapping AIC, all SLOPEQ based. Lowest RMSPE is manual m4
#however, SLOPEQ b5 seems to be similar in RMSPE. Let's do scenarios with this one
#first and then look at the others.
load('ssn1_glmssn5_SLOPEQ_20151106.Rdata')
assign(x = "ssn1_glmssn5", value = tmp)
rm(tmp)
