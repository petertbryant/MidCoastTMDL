library(SSN)

options(scipen = 100, stringsAsFactors = FALSE)

acs <- expand.grid(c("Exponential.Euclid", "Spherical.Euclid", 
                     "Gaussian.Euclid", "Cauchy.Euclid", NA), 
                   c("Exponential.tailup", "LinearSill.tailup", 
                     "Spherical.tailup", "Mariah.tailup", "Epanech.tailup", NA), 
                   c("Exponential.taildown", "LinearSill.taildown", 
                     "Spherical.taildown", "Mariah.taildown", "Epanech.taildown", NA),
                   stringsAsFactors = FALSE)

for (i in 1:(nrow(acs)/10)) {
  start.time <- Sys.time()
  print(paste("Starting block", i, "at", start.time))
  
  ind <- (10*i - 9):(10*i)
  
  object_names <- paste0("fit", ind)
  ssn_list <- as.list(object_names)
  names(ssn_list) <- object_names
  
  for (j in 1:10) {
    if (any(is.na(acs[ind[j],]))) {
      cm <- c("locID", unlist(acs[ind[j],!is.na(acs[ind[j],])], use.names = FALSE))
    } else {
      cm <- c("locID", acs[ind[j],'Var1'], acs[ind[j],"Var2"], acs[ind[j],"Var3"])
    }
    
    model.start.time <- Sys.time()
    print(paste("Starting model fit", ind[j], "at", model.start.time))
    ssn_list[[j]] <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + KFACT_MARCA + 
                              OWN_FED_PRCA + DIS_3YR_PRSA + ROADLEN_DRSA + OWN_URB_PARCA + HDWTR,
                            EstMeth = "ML",
                            ssn1,
                            CorModels = cm,
                            addfunccol = "afvArea",
                            family = "Gaussian")
    met <- Sys.time()
    print(paste("Completed model fit", ind[j], "at", met))
    print(paste("Elapsed time was", met - model.start.time))
    cat("\n")
  }
  
  sbs <- Sys.time()
  print(paste("Starting summary of block", i, "at", sbs))
  cat("\n")
  
  tmp_results <- InfoCritCompare(ssn_list)
  
  if (i == 1) {
    results <- tmp_results
  } else {
    results <- rbind(results, tmp_results)
  }
  
  rm(tmp_results, ssn_list)
  end.time <- Sys.time()
  print(paste("Completed block", i, "at", end.time))
  print(paste("elapsed time was", end.time - start.time))
  cat("\n\n")
}

results$dAIC <- min(results$AIC) - results$AIC
results <- results[order(results$dAIC, decreasing = TRUE),]
save(results, file = "autocorrelation_function_compare.Rdata")
write.csv(results, "autocorrelation_function_compare.Rdata")
