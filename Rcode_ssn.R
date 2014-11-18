library(SSN)
library(plyr)

ssn1 <- importSSN('bugs.ssn')
obs <- getSSNdata.frame(ssn1, Name = "Obs")
obs <- rename(obs, c(
              "STATION_KE" = "STATION_KEY",
              "APOPRCA201" = "APOPRCA2010",
              "sum_1095_d" = "sum_1095_days",
              "FSS_26Aug1" = "FSS_26Aug14",
              "PALITHEROD" = "PALITHERODRCA",
              "PADISRSA_1" = "PADISRSA_1YR",
              "PASUSCEP5_" = "PASUSCEP5_DE", 
              "POWNRCA_PR" = "POWNRCA_PRI",
              "POWNRCA_FE" = "POWNRCA_FED",
              "PAOWNRCA_A" = "PAOWNRCA_AGR",
              "log10_FSS_" = "log10_FSS_26Aug14",
              "log10_sum_" = "log10_sum_1095_days",
              "sqrt_PADIS" = "sqrt_PADISRSA_1YR",
              "bin_PALITH" = "bin_PALITHERODRCA",
              "log10_XSLO" = "log10_XSLOPE_MAP",
              "log10_PASI" = "log10_PASILTRCA",
              "log10_MIN_" = "log10_MIN_Z",
              "log10_STRM" = "log10_STRMPWR",
              "sqrt_upDis" = "sqrt_upDist",
              "log10_APOP" = "log10_APOPRCA2010",
              "log10_PASU" = "log10_PASUSCEP5_DE",
              "log10_POWN" = "log10_POWNRCA_PRI",
            "log10_POWN.1" = "log10_POWNRCA_FED",
              "bin_PAOWNR" = "bin_PAOWNRCA_AGR"))
ssn1 <- putSSNdata.frame(obs, ssn1, Name = 'Obs')

###################################################
### Create the Distance Matrix
###################################################
createDistMat(ssn1, o.write = TRUE)

###################################################
# map the stream network with observed sites
###################################################
# plot(ssn1, lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "blue",
#      pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)",
#      asp = 1)

###################################################
### plot the Torgegram
###################################################
# ssn1.Torg <- Torgegram(ssn1, "log10_FSS_26Aug14", nlag = 15, maxlag = 50000, nlagcutoff=5)
# plot(ssn1.Torg)

###################################################
### run the model 
###################################################
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1 <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + bin_PALITHERODRCA +  
                        sqrt_PADISRSA_1YR, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Exponential.tailup", "Exponential.taildown"),
                       addfunccol = "afvArea",
                       family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1)

ssn1.glmssn2 <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + 
                         sqrt_PADISRSA_1YR + bin_PALITHERODRCA + log10_XSLOPE_MAP + 
                         log10_PASILTRCA + log10_MIN_Z, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Exponential.tailup", "Exponential.taildown",
                                     "Exponential.Euclid"), addfunccol = "afvArea",family = "Gaussian")

summary(ssn1.glmssn2)

###################################################
### check the residuals
###################################################
ssn1.resid1 <- residuals(ssn1.glmssn1)
names( getSSNdata.frame(ssn1.resid1) )
plot(ssn1.resid1)

###################################################
### plot the residuals
###################################################
par(mfrow = c(1, 2))
hist(ssn1.resid1)
hist(ssn1, "log10_FSS_26Aug14")

###################################################
### cross validation
###################################################
cv.out <- CrossValidationSSN(ssn1.glmssn1)
par(mfrow = c(1, 2))
plot(ssn1.glmssn1$sampinfo$z,
     cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(ssn1)[, "FSS_26Aug14"]),
      cv.out[, "cv.se"], pch = 19,
      xlab = "Observed Data", ylab = "LOOCV Prediction SE")


cv.out <- CrossValidationSSN(ssn1.glmssn2)
par(mfrow = c(1, 2))
plot(ssn1.glmssn2$sampinfo$z,
     cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(ssn1)[, "log10_FSS_26Aug14"]),
      cv.out[, "cv.se"], pch = 19,
      xlab = "Observed Data", ylab = "LOOCV Prediction SE")
###################################################
### cross validation stats
###################################################
CrossValidationStatsSSN(ssn1.glmssn1)

###################################################
### R2
###################################################
GR2(ssn1.glmssn1)
varcomp(ssn1.glmssn1)

###################################################
### run the model
###########
########### THIS IS WHERE IT BREAKS ###############
###########
###################################################
ssn1.glmssn2 <- glmssn(FSS_26Aug14  ~ sum_1095_days + 
                         PADISRSA_1YR + PALITHERODRCA + XSLOPE_MAP + 
                         PASILTRCA + MIN_Z, 
                       ssn1,
                       CorModels = c("Exponential.tailup", "Exponential.taildown"), 
                       addfunccol = "afvArea",
                       family = "Poisson")

summary(ssn1.glmssn2)

###################################################
### check the residuals
###################################################
ssn1.resid2 <- residuals(ssn1.glmssn2)
names( getSSNdata.frame(ssn1.resid2) )
#plot(ssn1.resid1)

###################################################
### plot the residuals
###################################################
par(mfrow = c(1, 2))
hist(ssn1.resid2)
hist(ssn1, "log10_FSS_26Aug14")

###################################################
### cross validation
###################################################
cv.out2 <- CrossValidationSSN(ssn1.glmssn2)
par(mfrow = c(1, 2))
plot(ssn1.glmssn2$sampinfo$z,
     cv.out2[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(ssn1)[, "log10_FSS_26Aug14"]),
      cv.out2[, "cv.se"], pch = 19,
      xlab = "Observed Data", ylab = "LOOCV Prediction SE")

###################################################
### cross validation stats
###################################################
CrossValidationStatsSSN(ssn1.glmssn2)

###################################################
### R2
###################################################
GR2(ssn1.glmssn2)
varcomp(ssn1.glmssn2)
