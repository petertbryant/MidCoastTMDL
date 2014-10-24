library(SSN)
library(plyr)

ssn1 <- importSSN('bugs.ssn')
obs <- getSSNdata.frame(ssn1, Name = "Obs")
obs <- rename(obs, c(
              "STATION_KE" = "STATION_KEY",
              "APOPRCA201" = "APOPRCA2010",
              "sum_60_day" = "sum_60_days",
              "sum_180_da" = "sum_180_days",
              "sum_365_da" = "sum_365_days",
              "sum_1095_d" = "sum_1095_days",
              "PPT_1981_2" = "PPT_1981_2010",
              "FSS_26Aug1" = "FSS_26Aug14",
              "PLITHERODR" = "PLITHERODRCA",
              "PLITHERODR.1" = "PLITHERODRSA",
              "PSILTRCA" = "PSILTRCA",
              "PDISRCA_1Y" = "PDISRCA_1YR",
              "PDISRCA_3Y" = "PDISRCA_3YR",
              "PDISRSA_1Y" = "PDISRSA_1YR",
              "POWNRCA_PR" = "POWNRCA_PRI",
              "POWNRCA_FE" = "POWNRCA_FED",
              "POWNRSA_FE" = "POWNRSA_FED",
              "PALITHEROD" = "PALITHERODRCA",
              "PALITHEROD.1" = "PALITHERODRSA",
              "PASILT_CLA" = "PASILT_CLAYRCA",
              "PADISRSA_1" = "PADISRSA_1YR",
              "PAOWNRCA_F" = "PAOWNRCA_FED",
              "PAOWNRCA_A" = "PAOWNRCA_AGR",
              "PAOWNRSA_F" = "PAOWNRSA_FED",
              "PAOWNRSA_A" = "PAOWNRSA_AGR",
              "log10_FSS_" = "log10_FSS_26Aug14",
              "log10_sum_" = "log10_sum_1095_days",
              "sqrt_PADIS" = "sqrt_PADISRSA_1YR",
              "bin_PALITH" = "bin_PALITHERODRCA",
              "log10_XSLO" = "log10_XSLOPE_MAP",
              "log10_PASI" = "log10_PASILTRCA",
              "log10_MIN_" = "log10_MIN_Z"))
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
ssn1.glmssn1 <- glmssn(FSS_26Aug14  ~ sum_1095_days + 
                         PADISRSA_1YR + PALITHERODRCA + XSLOPE_MAP + 
                         PASILTRCA + MIN_Z, 
                       ssn1,
                       EstMeth = "REML",
                        CorModels = c("Exponential.tailup", "Exponential.taildown",
                                      "Exponential.Euclid"), addfunccol = "afvArea",family = "Gaussian")

summary(ssn1.glmssn1)

###################################################
### check the residuals
###################################################
ssn1.resid1 <- residuals(ssn1.glmssn1)
names( getSSNdata.frame(ssn1.resid1) )
#plot(ssn1.resid1)

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
