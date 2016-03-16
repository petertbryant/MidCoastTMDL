###################################################
### check the residuals   
###################################################
fit_resid <- residuals(fit)
names( getSSNdata.frame(fit_resid) )
plot(fit_resid)

resids <- getSSNdata.frame(fit_resid, Name = "Obs")


###################################################
### plot the residuals
###################################################
png('residuals.png', width = 6, height = 6, units = 'in', res = 200)
par(mfrow = c(2, 2))
hist(fit_resid, xlab = "Residuals")
#hist(ssn1, "log10_FSS_26Aug14", xlab = 'Observed log10 FSS')
plot(resids$"_fit_",resids$"_resid_", xlab = 'Predicted log10 BSTI', ylab = 'Raw residuals')
plot(resids$"_fit_",resids$"_resid.stand_", xlab = 'Predicted log10 BSTI', ylab = 'Standardized residuals')
qqnorm(resids$"_resid.stand_", ylab = 'Standardized residuals')
abline(0,1)
dev.off()

###################################################
### cross validation
###################################################
cv.out <- CrossValidationSSN(fit)
png('LOOCV.png', width = 6, height = 4, units = 'in', res = 100)
par(mfrow = c(1, 2))
plot(fit$sampinfo$z,
     cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction", ylim = c(0,105))
abline(0, 1)
plot( na.omit( getSSNdata.frame(ssn1)[, "BSTI"]),
      cv.out[, "cv.se"], pch = 19,
      xlab = "Observed Data", ylab = "LOOCV Prediction SE")
dev.off()

###################################################
### likelihood ratio test
###################################################
source('lrtSSN.R')
fit_null <- glmssn(log10_BSTI ~ 1,
                   EstMeth = "REML",
                   ssn1,
                   CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
                   addfunccol = "afvArea",
                   family = "Gaussian")
fit_lrt <- lrtSSN(fit, fit_null)

###################################################
### non-spatial model
###################################################
fit_nonspatial <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + KFACT_MARCA + 
                         OWN_FED_PRCA + DIS_3YR_PRSA + ROADLEN_DRSA + OWN_URB_PARCA + HDWTR,
                       EstMeth = "REML",
                       ssn1,
                       CorModels = c("locID"),
                       addfunccol = "afvArea",
                       family = "Gaussian")
results_nsp <- InfoCritCompare(list(fit, fit_nonspatial))
results_nsp$dAIC <- min(results_nsp$AIC) - results_nsp$AIC
