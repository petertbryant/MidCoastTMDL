library(SSN)

name_vec <- paste0("ssn1_glmssn_forward", 1:10, "_20151019.Rdata")
name_vec <- c(name_vec, paste0("ssn1_glmssn", 1:4, "_20151019.Rdata"))
for (i in 1:length(name_vec)) {
  load(name_vec[i])
  assign(gsub("_20151019.Rdata","",name_vec[i]), tmp)
}
load("ssn1_glmssn5_20151019.Rdata")
results <- InfoCritCompare(list(ssn1_glmssn_forward1, ssn1_glmssn_forward2, 
                                ssn1_glmssn_forward3, ssn1_glmssn_forward4, 
                                ssn1_glmssn_forward5, ssn1_glmssn_forward6,
                                ssn1_glmssn_forward7, ssn1_glmssn_forward8,
                                ssn1_glmssn_forward9, ssn1_glmssn_forward10,
                                ssn1_glmssn1, ssn1_glmssn2, 
                                ssn1_glmssn3, ssn1_glmssn4,
                                ssn1_glmssn5))

results$dAIC <- min(results$AIC) - results$AIC
results[order(results$dAIC, decreasing = TRUE),]
save(results,file = "back_for_ward_results_20151019.Rdata")
