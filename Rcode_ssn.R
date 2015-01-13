#### SubSSN Test ####
library(SSN)
library(stringr)

options('scipen' = 100)

##################################################
### Create the Distance Matrix
###################################################
createDistMat(ssn1, predpts = 'preds', o.write = TRUE, amongpreds = TRUE)

###################################################
# map the stream network with observed sites
###################################################
# plot(ssn1, lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "blue",
#      pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)",
#      asp = 1)

###################################################
### plot the Torgegram
###################################################
ssn1.Torg <- Torgegram(ssn1, "log10_FSS_26Aug14", nlag = 15, maxlag = 50000, nlagcutoff=5)
plot(ssn1.Torg)

###################################################
### Variable Selection Via Backward Deletion   ###
###################################################

#### ssn1.glmssn1.G ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1.G <- glmssn(FSS_26Aug14  ~ LAT_RAW + sum_1095_days + XSLOPE_MAP + MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA +
                         PACLAYRCA + PASILTRCA + PASUSCEP5_DE + DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "REML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1.G)
varcomp(ssn1.glmssn1.G)
#Plot the residuals
ssn1.resid1.G <- residuals(ssn1.glmssn1.G)
par(mfrow = c(1,2))
hist(ssn1.resid1.G)
hist(ssn1, "FSS_26Aug14")

qqnorm(ssn1.resid1.G)

resids.df <- getSSNdata.frame(ssn1.resid1.G)
plot(resids.df[,"_fit_"],resids.df[,"_resid_"])
#wedge shape means not constant variance. per p.560 in R Book second edition we will log transform the response.
mean(obs.vars$FSS_26Aug14)
var(obs.vars$FSS_26Aug14)
#variance is greater than mean. Our data is similar to count data, however the variance is greater than the mean NOT equal
#to it as assumed by the poisson distribution (p. 569). so we will use the gaussian distribution.

save(ssn1.glmssn1.G, file = 'ssn1_glmssn1_G.Rdata.Rdata')

#### ssn1.glmssn1.G Summary ####
# Call:
#   glmssn(formula = FSS_26Aug14 ~ LAT_RAW + sum_1095_days + XSLOPE_MAP + 
#            MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
#            PACLAYRCA + PASILTRCA + PASUSCEP5_DE + DAROADX + DAPOPRCA2010 + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -16.97  -4.58  -0.91   2.99  54.10 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)   44.427402  57.688254    0.77               0.4415    
# LAT_RAW       -1.006621   1.309691   -0.77               0.4424    
# sum_1095_days -0.000977   0.000221   -4.43              0.00001 ***
#   XSLOPE_MAP    -0.122992   0.090818   -1.35               0.1761    
# MIN_Z          0.000269   0.001371    0.20               0.8446    
# STRMPWR       -0.001419   0.001359   -1.04               0.2970    
# PDISRSA_1YR    0.056909   0.020917    2.72               0.0067 ** 
#   POWNRCA_PRI    0.022052   0.010038    2.20               0.0283 *  
#   PALITHERODRCA  0.078008   0.014355    5.43 < 0.0000000000000002 ***
#   PACLAYRCA      0.007501   0.098835    0.08               0.9395    
# PASILTRCA      0.244496   0.081294    3.01               0.0027 ** 
#   PASUSCEP5_DE   0.046812   0.104323    0.45               0.6538    
# DAROADX       -0.150326   0.981782   -0.15               0.8783    
# DAPOPRCA2010  72.502280  12.970459    5.59 < 0.0000000000000002 ***
#   PAOWNRCA_AGR   0.391062   0.255309    1.53               0.1260    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter      Estimate
# Exponential.tailup   parsill      0.006647
# Exponential.tailup     range    785.518222
# Exponential.taildown   parsill     22.369338
# Exponential.taildown     range   7999.619997
# Exponential.Euclid   parsill     21.818488
# Exponential.Euclid     range 106065.304782
# locID   parsill      0.000985
# Nugget   parsill     25.932001
# 
# Residual standard error: 8.374
# Generalized R-squared: 0.1923
# 
# varcomp(ssn1.glmssn1.G)
# VarComp Proportion
# 1    Covariates (R-sq) 0.19233360
# 2   Exponential.tailup 0.00007655
# 3 Exponential.taildown 0.25763036
# 4   Exponential.Euclid 0.25128616
# 5                locID 0.00001134
# 6               Nugget 0.29866198

#### ssn1.glmssn3.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now we run it with ML so we can use AIC for comparison.
start.time <- Sys.time()
print(start.time)
ssn1.glmssn3.G <- glmssn(log10_FSS_26Aug14  ~ LAT_RAW + sum_1095_days + XSLOPE_MAP + MIN_Z + STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + PASUSCEP5_DE + 
                           DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn3.G)
varcomp(ssn1.glmssn3.G)
#Plot the residuals
ssn1.resid2.G <- residuals(ssn1.glmssn3.G)
par(mfrow = c(1,2))
hist(ssn1.resid2.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid2.G)

resids2.df <- getSSNdata.frame(ssn1.resid2.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn3.G, file = 'ssn1_glmssn3_G.Rdata')

AIC(ssn1.glmssn3.G)
# -1216.322

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ LAT_RAW + sum_1095_days + 
#            XSLOPE_MAP + MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + 
#            PALITHERODRCA + PACLAYRCA + PASILTRCA + PASUSCEP5_DE + DAROADX + 
#            DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41270 -0.06425  0.01227  0.09220  0.47498 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.39063    0.06499   6.011 < 0.0000000000000002 ***
#   LAT_RAW        0.07282    0.09004   0.809              0.41891    
# sum_1095_days -0.30144    0.03674  -8.205 < 0.0000000000000002 ***
#   XSLOPE_MAP    -0.01542    0.04372  -0.353              0.72441    
# MIN_Z         -0.03382    0.05617  -0.602              0.54734    
# STRMPWR       -0.10247    0.10249  -1.000              0.31772    
# PDISRSA_1YR    0.08261    0.02903   2.846              0.00455 ** 
#   POWNRCA_PRI    0.05642    0.01554   3.632              0.00030 ***
#   PALITHERODRCA  0.13426    0.02147   6.254 < 0.0000000000000002 ***
#   PACLAYRCA      0.08025    0.04825   1.663              0.09664 .  
# PASILTRCA      0.09275    0.04367   2.124              0.03401 *  
#   PASUSCEP5_DE   0.02404    0.04556   0.528              0.59794    
# DAROADX        0.10404    0.12990   0.801              0.42342    
# DAPOPRCA2010   0.18214    0.04928   3.696              0.00023 ***
#   PAOWNRCA_AGR   0.07511    0.07969   0.943              0.34619    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0014070
# Exponential.tailup     range    614.2147098
# Exponential.taildown   parsill      0.0048045
# Exponential.taildown     range  12009.8212620
# Exponential.Euclid   parsill      0.0042566
# Exponential.Euclid     range 172375.7611170
# locID   parsill      0.0000187
# Nugget   parsill      0.0061414
# 
# Residual standard error: 0.1289502
# Generalized R-squared: 0.2582293

#Scaling the data left p-values the same. Residual standard error is halved. GR2 is the same. AIC way lower.
#Below, see that untransformed RMSE is improved to 3.6
#Effect = No diff for variable selection.

# > varcomp(ssn1.glmssn3.G)
# VarComp   Proportion
# 1    Covariates (R-sq) 0.2582292526
# 2   Exponential.tailup 0.0627656026
# 3 Exponential.taildown 0.2143258485
# 4   Exponential.Euclid 0.1898836200
# 5                locID 0.0008336874
# 6               Nugget 0.2739619889

#### ssn1.glmssn4.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#LAT_RAW not significant. Removed.
start.time <- Sys.time()
print(start.time)
ssn1.glmssn4.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + XSLOPE_MAP + MIN_Z + STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + PASUSCEP5_DE + 
                           DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn4.G)
varcomp(ssn1.glmssn4.G)
#Plot the residuals
# ssn1.resid4.G <- residuals(ssn1.glmssn4.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid4.G)
# hist(ssn1, "log10_FSS_26Aug14")

# qqnorm(ssn1.resid4.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid4.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn4.G, file = 'ssn1_glmssn4_G.Rdata')

AIC(ssn1.glmssn4.G)
# -1217.678

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + XSLOPE_MAP + 
#            MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
#            PACLAYRCA + PASILTRCA + PASUSCEP5_DE + DAROADX + DAPOPRCA2010 + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.39605 -0.06372  0.01655  0.09277  0.45899 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.42941    0.04448   9.655 < 0.0000000000000002 ***
#   sum_1095_days -0.30238    0.03680  -8.216 < 0.0000000000000002 ***
#   XSLOPE_MAP    -0.01588    0.04373  -0.363              0.71657    
# MIN_Z         -0.03224    0.05624  -0.573              0.56657    
# STRMPWR       -0.10343    0.10247  -1.009              0.31315    
# PDISRSA_1YR    0.08028    0.02894   2.774              0.00568 ** 
#   POWNRCA_PRI    0.05737    0.01551   3.699              0.00023 ***
#   PALITHERODRCA  0.13325    0.02152   6.193 < 0.0000000000000002 ***
#   PACLAYRCA      0.07366    0.04763   1.546              0.12240    
# PASILTRCA      0.10118    0.04237   2.388              0.01719 *  
#   PASUSCEP5_DE   0.02174    0.04559   0.477              0.63364    
# DAROADX        0.08845    0.12863   0.688              0.49187    
# DAPOPRCA2010   0.18290    0.04927   3.712              0.00022 ***
#   PAOWNRCA_AGR   0.07555    0.07968   0.948              0.34336    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter         Estimate
# Exponential.tailup   parsill      0.001445902
# Exponential.tailup     range    609.274554634
# Exponential.taildown   parsill      0.004741695
# Exponential.taildown     range  11913.925055205
# Exponential.Euclid   parsill      0.004397020
# Exponential.Euclid     range 169286.468577436
# locID   parsill      0.000000837
# Nugget   parsill      0.006142826
# 
# Residual standard error: 0.1293379
# Generalized R-squared: 0.2559038
# > varcomp(ssn1.glmssn4.G)
# VarComp    Proportion
# 1    Covariates (R-sq) 0.25590380841
# 2   Exponential.tailup 0.06431563996
# 3 Exponential.taildown 0.21091690146
# 4   Exponential.Euclid 0.19558531596
# 5                locID 0.00003724791
# 6               Nugget 0.27324108630

#### ssn1.glmssn5.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without XSLOPE_MAP
start.time <- Sys.time()
print(start.time)
ssn1.glmssn5.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + MIN_Z + STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + PASUSCEP5_DE + 
                           DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn5.G)
varcomp(ssn1.glmssn5.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn5.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn5.G, file = 'ssn1_glmssn5_G.Rdata')

AIC(ssn1.glmssn5.G)
#-1219.545

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + MIN_Z + 
#            STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + 
#            PASILTRCA + PASUSCEP5_DE + DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.39532 -0.06381  0.01632  0.09253  0.45947 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.42707    0.04375   9.761 < 0.0000000000000002 ***
#   sum_1095_days -0.30212    0.03680  -8.210 < 0.0000000000000002 ***
#   MIN_Z         -0.03701    0.05463  -0.677              0.49831    
# STRMPWR       -0.11282    0.09912  -1.138              0.25538    
# PDISRSA_1YR    0.08163    0.02868   2.847              0.00454 ** 
#   POWNRCA_PRI    0.05798    0.01543   3.758              0.00018 ***
#   PALITHERODRCA  0.13376    0.02148   6.228 < 0.0000000000000002 ***
#   PACLAYRCA      0.07351    0.04763   1.543              0.12318    
# PASILTRCA      0.09953    0.04214   2.362              0.01843 *  
#   PASUSCEP5_DE   0.02120    0.04559   0.465              0.64201    
# DAROADX        0.07821    0.12576   0.622              0.53419    
# DAPOPRCA2010   0.18159    0.04917   3.693              0.00024 ***
#   PAOWNRCA_AGR   0.07797    0.07943   0.982              0.32660    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0014455
# Exponential.tailup     range    594.5771837
# Exponential.taildown   parsill      0.0047112
# Exponential.taildown     range  12035.0996124
# Exponential.Euclid   parsill      0.0043714
# Exponential.Euclid     range 166537.0691033
# locID   parsill      0.0000271
# Nugget   parsill      0.0061437
# 
# Residual standard error: 0.1292245
# Generalized R-squared: 0.255661
# > varcomp(ssn1.glmssn5.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.255660974
# 2   Exponential.tailup 0.064433324
# 3 Exponential.taildown 0.209997888
# 4   Exponential.Euclid 0.194852381
# 5                locID 0.001207459
# 6               Nugget 0.273847974

#### ssn1.glmssn6.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without PASUSCEP5_DE
start.time <- Sys.time()
print(start.time)
ssn1.glmssn6.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + MIN_Z + STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn6.G)
varcomp(ssn1.glmssn6.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn6.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn6.G, file = 'ssn1_glmssn6_G.Rdata')

AIC(ssn1.glmssn6.G)
#-1221.331

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + MIN_Z + 
#            STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + 
#            PASILTRCA + DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.39687 -0.06433  0.01446  0.09419  0.45974 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.43466    0.04057  10.714 < 0.0000000000000002 ***
#   sum_1095_days -0.30137    0.03678  -8.193 < 0.0000000000000002 ***
#   MIN_Z         -0.04131    0.05376  -0.768              0.44251    
# STRMPWR       -0.10880    0.09876  -1.102              0.27097    
# PDISRSA_1YR    0.08133    0.02868   2.836              0.00469 ** 
#   POWNRCA_PRI    0.05721    0.01534   3.729              0.00021 ***
#   PALITHERODRCA  0.13194    0.02117   6.232 < 0.0000000000000002 ***
#   PACLAYRCA      0.06571    0.04471   1.470              0.14207    
# PASILTRCA      0.09665    0.04181   2.311              0.02107 *  
#   DAROADX        0.07653    0.12582   0.608              0.54321    
# DAPOPRCA2010   0.17948    0.04900   3.663              0.00027 ***
#   PAOWNRCA_AGR   0.07977    0.07938   1.005              0.31522    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0014291
# Exponential.tailup     range    593.3115380
# Exponential.taildown   parsill      0.0046824
# Exponential.taildown     range  12002.8962859
# Exponential.Euclid   parsill      0.0043734
# Exponential.Euclid     range 160238.4666463
# locID   parsill      0.0000417
# Nugget   parsill      0.0061461
# 
# Residual standard error: 0.129123
# Generalized R-squared: 0.2546098
# > varcomp(ssn1.glmssn6.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.254609791
# 2   Exponential.tailup 0.063892169
# 3 Exponential.taildown 0.209334689
# 4   Exponential.Euclid 0.195522394
# 5                locID 0.001864812
# 6               Nugget 0.274776144

#### ssn1.glmssn7.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without DAROADX
start.time <- Sys.time()
print(start.time)
ssn1.glmssn7.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + MIN_Z + STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn7.G)
varcomp(ssn1.glmssn7.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn7.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn7.G, file = 'ssn1_glmssn7_G.Rdata')

AIC(ssn1.glmssn7.G)
#-1222.975

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + MIN_Z + 
#            STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + 
#            PASILTRCA + DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.39781 -0.06590  0.01410  0.09299  0.45842 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.43690    0.03955  11.047 < 0.0000000000000002 ***
#   sum_1095_days -0.30248    0.03675  -8.231 < 0.0000000000000002 ***
#   MIN_Z         -0.04112    0.05374  -0.765              0.44438    
# STRMPWR       -0.10992    0.09880  -1.113              0.26626    
# PDISRSA_1YR    0.08066    0.02868   2.813              0.00504 ** 
#   POWNRCA_PRI    0.05798    0.01530   3.790              0.00016 ***
#   PALITHERODRCA  0.13101    0.02109   6.212 < 0.0000000000000002 ***
#   PACLAYRCA      0.06653    0.04460   1.492              0.13616    
# PASILTRCA      0.09590    0.04167   2.301              0.02164 *  
#   DAPOPRCA2010   0.18260    0.04880   3.742              0.00020 ***
#   PAOWNRCA_AGR   0.07902    0.07940   0.995              0.31996    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00146915
# Exponential.tailup     range    589.47313947
# Exponential.taildown   parsill      0.00467840
# Exponential.taildown     range  12035.63406844
# Exponential.Euclid   parsill      0.00419659
# Exponential.Euclid     range 150454.28068083
# locID   parsill      0.00000344
# Nugget   parsill      0.00615042
# 
# Residual standard error: 0.1284445
# Generalized R-squared: 0.2547728
# > varcomp(ssn1.glmssn7.G)
# VarComp   Proportion
# 1    Covariates (R-sq) 0.2547727624
# 2   Exponential.tailup 0.0663624945
# 3 Exponential.taildown 0.2113270009
# 4   Exponential.Euclid 0.1895631615
# 5                locID 0.0001555511
# 6               Nugget 0.2778190295

#### ssn1.glmssn8.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without MIN_Z 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn8.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn8.G)
varcomp(ssn1.glmssn8.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn8.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn8.G, file = 'ssn1_glmssn8_G.Rdata')

AIC(ssn1.glmssn8.G)
#[1] -1222.786

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41541 -0.06222  0.01672  0.09593  0.46353 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.43034    0.04143  10.388 < 0.0000000000000002 ***
#   sum_1095_days -0.30272    0.03683  -8.219 < 0.0000000000000002 ***
#   STRMPWR       -0.11134    0.09859  -1.129              0.25911    
# PDISRSA_1YR    0.08292    0.02855   2.904              0.00378 ** 
#   POWNRCA_PRI    0.05810    0.01532   3.792              0.00016 ***
#   PALITHERODRCA  0.13165    0.02130   6.181 < 0.0000000000000002 ***
#   PACLAYRCA      0.06614    0.04492   1.472              0.14131    
# PASILTRCA      0.09656    0.04220   2.288              0.02240 *  
#   DAROADX        0.07628    0.12626   0.604              0.54594    
# DAPOPRCA2010   0.17958    0.04901   3.664              0.00026 ***
#   PAOWNRCA_AGR   0.07932    0.07933   1.000              0.31769    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00157901
# Exponential.tailup     range    582.71052394
# Exponential.taildown   parsill      0.00452800
# Exponential.taildown     range  12192.93060979
# Exponential.Euclid   parsill      0.00480547
# Exponential.Euclid     range 165412.95661634
# locID   parsill      0.00000328
# Nugget   parsill      0.00612703
# 
# Residual standard error: 0.130548
# Generalized R-squared: 0.2506058
# > varcomp(ssn1.glmssn8.G)
# VarComp   Proportion
# 1    Covariates (R-sq) 0.2506057658
# 2   Exponential.tailup 0.0694314005
# 3 Exponential.taildown 0.1991021287
# 4   Exponential.Euclid 0.2113029196
# 5                locID 0.0001440143
# 6               Nugget 0.2694137711

#### ssn1.glmssn9.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without DAROADX 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn9.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn9.G)
varcomp(ssn1.glmssn9.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn9.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn9.G, file = 'ssn1_glmssn9_G.Rdata')

AIC(ssn1.glmssn9.G)
#-1224.419

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41630 -0.06364  0.01547  0.09445  0.46216 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.43263    0.04034  10.724 < 0.0000000000000002 ***
#   sum_1095_days -0.30387    0.03680  -8.256 < 0.0000000000000002 ***
#   STRMPWR       -0.11238    0.09862  -1.140              0.25484    
# PDISRSA_1YR    0.08223    0.02854   2.882              0.00407 ** 
#   POWNRCA_PRI    0.05887    0.01528   3.853              0.00013 ***
#   PALITHERODRCA  0.13072    0.02122   6.161 < 0.0000000000000002 ***
#   PACLAYRCA      0.06694    0.04481   1.494              0.13556    
# PASILTRCA      0.09588    0.04207   2.279              0.02294 *  
#   DAPOPRCA2010   0.18266    0.04880   3.743              0.00020 ***
#   PAOWNRCA_AGR   0.07856    0.07935   0.990              0.32247    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0014862
# Exponential.tailup     range    590.4341375
# Exponential.taildown   parsill      0.0045251
# Exponential.taildown     range  12202.7781303
# Exponential.Euclid   parsill      0.0046173
# Exponential.Euclid     range 155813.6299332
# locID   parsill      0.0000913
# Nugget   parsill      0.0061345
# 
# Residual standard error: 0.1298247
# Generalized R-squared: 0.250746
# > varcomp(ssn1.glmssn9.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.250746008
# 2   Exponential.tailup 0.066066458
# 3 Exponential.taildown 0.201161058
# 4   Exponential.Euclid 0.205259237
# 5                locID 0.004060033
# 6               Nugget 0.272707206

#### ssn1.glmssn10.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without PAOWNRCA_AGR 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn10.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAROADX + DAPOPRCA2010, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn10.G)
varcomp(ssn1.glmssn10.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn10.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn10.G, file = 'ssn1_glmssn10_G.Rdata')

AIC(ssn1.glmssn10.G)
#-1223.791

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAROADX + DAPOPRCA2010, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41496 -0.06142  0.01644  0.09660  0.46569 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.42519    0.04121  10.317 < 0.0000000000000002 ***
#   sum_1095_days -0.30273    0.03679  -8.229 < 0.0000000000000002 ***
#   STRMPWR       -0.10876    0.09866  -1.102              0.27061    
# PDISRSA_1YR    0.08878    0.02804   3.166              0.00160 ** 
#   POWNRCA_PRI    0.05719    0.01530   3.739              0.00020 ***
#   PALITHERODRCA  0.13292    0.02122   6.264 < 0.0000000000000002 ***
#   PACLAYRCA      0.07341    0.04425   1.659              0.09753 .  
# PASILTRCA      0.10305    0.04167   2.473              0.01362 *  
#   DAROADX        0.07398    0.12617   0.586              0.55782    
# DAPOPRCA2010   0.18666    0.04858   3.842              0.00013 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0015404
# Exponential.tailup     range    588.6529185
# Exponential.taildown   parsill      0.0045922
# Exponential.taildown     range  12220.6824017
# Exponential.Euclid   parsill      0.0047031
# Exponential.Euclid     range 170100.7003100
# locID   parsill      0.0000412
# Nugget   parsill      0.0061317
# 
# Residual standard error: 0.1304167
# Generalized R-squared: 0.2512154
# > varcomp(ssn1.glmssn10.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.251215396
# 2   Exponential.tailup 0.067812695
# 3 Exponential.taildown 0.202167840
# 4   Exponential.Euclid 0.207048256
# 5                locID 0.001814067
# 6               Nugget 0.269941746

#### ssn1.glmssn11.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without DAROADX 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn11.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAPOPRCA2010 + PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn11.G)
varcomp(ssn1.glmssn11.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn11.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn11.G, file = 'ssn1_glmssn11_G.Rdata')

AIC(ssn1.glmssn11.G)
# -1224.419

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41630 -0.06364  0.01547  0.09445  0.46216 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.43263    0.04034  10.724 < 0.0000000000000002 ***
#   sum_1095_days -0.30387    0.03680  -8.256 < 0.0000000000000002 ***
#   STRMPWR       -0.11238    0.09862  -1.140              0.25484    
# PDISRSA_1YR    0.08223    0.02854   2.882              0.00407 ** 
#   POWNRCA_PRI    0.05887    0.01528   3.853              0.00013 ***
#   PALITHERODRCA  0.13072    0.02122   6.161 < 0.0000000000000002 ***
#   PACLAYRCA      0.06694    0.04481   1.494              0.13556    
# PASILTRCA      0.09588    0.04207   2.279              0.02294 *  
#   DAPOPRCA2010   0.18266    0.04880   3.743              0.00020 ***
#   PAOWNRCA_AGR   0.07856    0.07935   0.990              0.32247    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0014862
# Exponential.tailup     range    590.4341375
# Exponential.taildown   parsill      0.0045251
# Exponential.taildown     range  12202.7781303
# Exponential.Euclid   parsill      0.0046173
# Exponential.Euclid     range 155813.6299332
# locID   parsill      0.0000913
# Nugget   parsill      0.0061345
# 
# Residual standard error: 0.1298247
# Generalized R-squared: 0.250746
# > varcomp(ssn1.glmssn11.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.250746008
# 2   Exponential.tailup 0.066066458
# 3 Exponential.taildown 0.201161058
# 4   Exponential.Euclid 0.205259237
# 5                locID 0.004060033
# 6               Nugget 0.272707206

#### ssn1.glmssn12.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without PAOWNRCA_AGR 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn12.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  STRMPWR + 
                           PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                           DAPOPRCA2010, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn12.G)
varcomp(ssn1.glmssn12.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn12.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn12.G, file = 'ssn1_glmssn12_G.Rdata')

AIC(ssn1.glmssn12.G)
#-1225.454

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAPOPRCA2010, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41584 -0.06307  0.01445  0.09537  0.46433 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.42743    0.04009  10.661 < 0.0000000000000002 ***
#   sum_1095_days -0.30382    0.03675  -8.267 < 0.0000000000000002 ***
#   STRMPWR       -0.10982    0.09870  -1.113              0.26621    
# PDISRSA_1YR    0.08810    0.02803   3.144              0.00173 ** 
#   POWNRCA_PRI    0.05795    0.01525   3.799              0.00016 ***
#   PALITHERODRCA  0.13205    0.02113   6.250 < 0.0000000000000002 ***
#   PACLAYRCA      0.07405    0.04413   1.678              0.09374 .  
# PASILTRCA      0.10239    0.04152   2.466              0.01388 *  
#   DAPOPRCA2010   0.18965    0.04838   3.920              0.00010 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0015210
# Exponential.tailup     range    581.1283878
# Exponential.taildown   parsill      0.0045996
# Exponential.taildown     range  12188.3334097
# Exponential.Euclid   parsill      0.0045026
# Exponential.Euclid     range 160554.5083216
# locID   parsill      0.0000557
# Nugget   parsill      0.0061352
# 
# Residual standard error: 0.129669
# Generalized R-squared: 0.2515232
# > varcomp(ssn1.glmssn12.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.251523162
# 2   Exponential.tailup 0.067705195
# 3 Exponential.taildown 0.204749399
# 4   Exponential.Euclid 0.200434188
# 5                locID 0.002479688
# 6               Nugget 0.273108369

#### ssn1.glmssn13.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without STRMPWR 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn13.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "ML",
                          CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn13.G)
varcomp(ssn1.glmssn13.G)
#Plot the residuals
# ssn1.resid5.G <- residuals(ssn1.glmssn13.G)
# par(mfrow = c(1,2))
# hist(ssn1.resid5.G)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid5.G)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid5.G)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn13.G, file = 'ssn1_glmssn13_G.Rdata')

AIC(ssn1.glmssn13.G)
#-1226.24

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + 
#            POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + DAPOPRCA2010, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41905 -0.06406  0.01493  0.09443  0.46434 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.42313    0.04044  10.463 < 0.0000000000000002 ***
#   sum_1095_days -0.30621    0.03676  -8.331 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.08868    0.02804   3.163              0.00162 ** 
#   POWNRCA_PRI    0.05658    0.01520   3.722              0.00021 ***
#   PALITHERODRCA  0.13315    0.02113   6.301 < 0.0000000000000002 ***
#   PACLAYRCA      0.07617    0.04417   1.725              0.08500 .  
# PASILTRCA      0.10199    0.04166   2.448              0.01459 *  
#   DAPOPRCA2010   0.18762    0.04834   3.881              0.00011 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00166019
# Exponential.tailup     range    570.54496898
# Exponential.taildown   parsill      0.00447740
# Exponential.taildown     range  11950.96404640
# Exponential.Euclid   parsill      0.00469486
# Exponential.Euclid     range 162851.86377026
# locID   parsill      0.00000134
# Nugget   parsill      0.00613858
# 
# Residual standard error: 0.130278
# Generalized R-squared: 0.2494496
# > varcomp(ssn1.glmssn13.G)
# VarComp    Proportion
# 1    Covariates (R-sq) 0.24944955133
# 2   Exponential.tailup 0.07341681353
# 3 Exponential.taildown 0.19799898266
# 4   Exponential.Euclid 0.20761557899
# 5                locID 0.00005941125
# 6               Nugget 0.27145966225

#### ssn1.glmssn14.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without PACLAYRCA 
start.time <- Sys.time()
print(start.time)
ssn1.glmssn14.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "ML",
                          CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn14.G)
varcomp(ssn1.glmssn14.G)
#Plot the residuals
ssn1.resid5.G <- residuals(ssn1.glmssn14.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn14.G, file = 'ssn1_glmssn14_G.Rdata')

AIC(ssn1.glmssn14.G)
#-1225.341

#This one had the next to highest AIC but was not significantly different. In the interest of parsimony this model will be the
#final model.

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + 
#            POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.43014 -0.06432  0.01135  0.09206  0.44924 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.44297    0.03642  12.161 < 0.0000000000000002 ***
#   sum_1095_days -0.31639    0.03629  -8.719 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.08801    0.02813   3.129              0.00182 ** 
#   POWNRCA_PRI    0.05862    0.01520   3.857              0.00012 ***
#   PALITHERODRCA  0.13621    0.02107   6.466 < 0.0000000000000002 ***
#   PASILTRCA      0.11911    0.04023   2.961              0.00316 ** 
#   DAPOPRCA2010   0.19232    0.04841   3.973              0.00008 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0016560
# Exponential.tailup     range    599.7843407
# Exponential.taildown   parsill      0.0045089
# Exponential.taildown     range  11580.1162327
# Exponential.Euclid   parsill      0.0042054
# Exponential.Euclid     range 142016.4350760
# locID   parsill      0.0000354
# Nugget   parsill      0.0061395
# 
# Residual standard error: 0.1286279
# Generalized R-squared: 0.2490099
# > varcomp(ssn1.glmssn14.G)
# VarComp  Proportion
# 1    Covariates (R-sq) 0.249009895
# 2   Exponential.tailup 0.075166931
# 3 Exponential.taildown 0.204661601
# 4   Exponential.Euclid 0.190883757
# 5                locID 0.001605222
# 6               Nugget 0.278672595

#### ssn1.glmssn2.NULL ####
#Generate a NULL model in order to perform Likelihood Ratio Test
start.time <- Sys.time()
print(start.time)
ssn1.glmssn2.NULL <- glmssn(log10_FSS_26Aug14 ~ 1, 
                            ssn1,
                            EstMeth = "ML",
                            CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                            addfunccol = "afvArea",
                            family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn2.NULL)
varcomp(ssn1.glmssn2.NULL)
#Plot the residuals
# ssn1.resid2.NULL <- residuals(ssn1.glmssn2.NULL)
# par(mfrow = c(1,2))
# hist(ssn1.resid2.NULL)
# hist(ssn1, "log10_FSS_26Aug14")
# 
# qqnorm(ssn1.resid2.NULL)
# 
# resids2.df <- getSSNdata.frame(ssn1.resid2.NULL)
# plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn2.NULL, file = 'ssn1_glmssn2_NULL.Rdata')

AIC(ssn1.glmssn2.NULL)
# -119.0271

#When no explanatory variables included we see that Euclidean autocorrelation comprises most of the variance explanation
#while tail down and the nugget come in second and third

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ 1, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.99240 -0.21425 -0.03816  0.18369  0.84644 
# 
# Coefficients:
#   Estimate Std. Error t value            Pr(>|t|)    
# (Intercept)   0.9924     0.1183   8.389 <0.0000000000000002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00563612
# Exponential.tailup     range   1035.94038161
# Exponential.taildown   parsill      0.01659172
# Exponential.taildown     range  16233.75320956
# Exponential.Euclid   parsill      0.10033440
# Exponential.Euclid     range 141107.96923950
# locID   parsill      0.00000513
# Nugget   parsill      0.02057526
# 
# Residual standard error: 0.378342
# Generalized R-squared: -0.00000000000000288658
# > varcomp(ssn1.glmssn2.NULL)
# VarComp              Proportion
# 1    Covariates (R-sq) -0.00000000000000288658
# 2   Exponential.tailup  0.03937415495179414499
# 3 Exponential.taildown  0.11591039811080697253
# 4   Exponential.Euclid  0.70094005434283590539
# 5                locID  0.00003583304288513766
# 6               Nugget  0.14373955955168074561

###################################################
###     Autocovariance function Selection       ###
###################################################

#### ssn1.glmssn.ELG ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ELG <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "LinearSill.taildown", "Gaussian.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.ELG, file = 'ssn1_glmssn_ELG.Rdata')

#### ssn1.glmssn.ELE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ELE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "LinearSill.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.ELE, file = 'ssn1_glmssn_ELE.Rdata')

#### ssn1.glmssn.EME ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.EME <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Mariah.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.EME, file = 'ssn1_glmssn_EME.Rdata')

#### ssn1.glmssn.ESE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ESE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Spherical.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.ESE, file = 'ssn1_glmssn_ESE.Rdata')

#### ssn1.glmssn.SSE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.SSE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.SSE, file = 'ssn1_glmssn_SSE.Rdata')

# start.time <- Sys.time()
# print(start.time)
# ssn1.glmssn.SSE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
#                             PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
#                             DAPOPRCA2010, 
#                           ssn1,
#                           EstMeth = "ML",
#                           CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
#                           addfunccol = "afvArea",
#                           family = "Gaussian")
# end.time <- Sys.time()
# print(end.time)
# print(end.time - start.time)
# save(ssn1.glmssn.SSE, file = 'ssn1_glmssn_SSE_scaled_ML.Rdata')


#### ssn1.glmssn.LLE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.LLE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","LinearSill.tailup", "LinearSill.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.LLE, file = 'ssn1_glmssn_LLE.Rdata')

#### ssn1.glmssn.MME ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.MME <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Mariah.tailup", "Mariah.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.MME, file = 'ssn1_glmssn_MME.Rdata')

#### ssn1.glmssn.MEE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.MEE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Mariah.tailup", "Exponential.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.MEE, file = 'ssn1_glmssn_MEE.Rdata')

#### ssn1.glmssn.SSS ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.SSS <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Spherical.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.SSS, file = 'ssn1_glmssn_SSS.Rdata')

#### ssn1.glmssn.ESG ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ESG <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Spherical.taildown", "Gaussian.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.ESG, file = 'ssn1_glmssn_ESG.Rdata')

#### ssn1.glmssn.EEE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.EEE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)
save(ssn1.glmssn.EEE, file = 'ssn1_glmssn_EEE.Rdata')

#############################
#### Ecological Approach ####
#############################
#The above uses a mathematical approach to model selection. The below uses a set of 
#variables that include potential interactions based on our knowledge of how landscape
#affect sediment dynamics. 
#WARNING: Takes ~14 hours to run
em <- read.csv('ecological_models.csv')
formulas <- as.character(em$Formula)
formulas <- formulas[formulas != ""]
formulas <- paste('log10_FSS_26Aug14 ~',formulas)
ssn.ecomod.list <- list()
for(i in 1:length(formulas)) {
  print(paste("Started model", i, "at", Sys.time()))
  fit <- glmssn(as.formula(formulas[i]),
                            ssn1,
                            EstMeth = "ML",
                            CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                            addfunccol = "afvArea",
                            family = "Gaussian")
  save(fit, file = paste("model",i,".Rdata",sep=''))
  ssn.ecomod.list <- c(ssn.ecomod.list, list(fit))
  print(paste("Finished model", i, "at", Sys.time()))
}
ecomod.compare <- InfoCritCompare(ssn.ecomod.list)
write.csv(ecomod.compare, 'ecological_models_comparison_01072015.csv')
save(ssn.ecomod.list, file = 'allModelsList_01072015.Rdata')

#started at 12:58
print(Sys.time)

#NULL model
ssn.SSE.NULL <- glmssn(log10_FSS_26Aug14 ~ 1,
              ssn1,
              EstMeth = "ML",
              CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
              addfunccol = "afvArea",
              family = "Gaussian")
save(ssn.SSE.NULL, file = 'ssn_SSE_NULL.Rdata')

###################################################
### check the residuals   
###################################################
ssn1.resid1 <- residuals(ssn1.glmssn14.G)
names( getSSNdata.frame(ssn1.resid1) )
plot(ssn1.resid1)

resids <- getSSNdata.frame(ssn1.resid1, Name = "Obs")


###################################################
### plot the residuals
###################################################
png('residuals.png', width = 6, height = 6, units = 'in', res = 200)
par(mfrow = c(2, 2))
hist(ssn1.resid1, xlab = "Residuals")
#hist(ssn1, "log10_FSS_26Aug14", xlab = 'Observed log10 FSS')
plot(resids$"_fit_",resids$"_resid_", xlab = 'Predicted log10 FSS', ylab = 'Raw residuals')
plot(resids$"_fit_",resids$"_resid.stand_", xlab = 'Predicted log10 FSS', ylab = 'Standardized residuals')
qqnorm(resids$"_resid.stand_", ylab = 'Standardized residuals')
abline(0,1)
dev.off()
###################################################
### cross validation
###################################################
cv.out <- CrossValidationSSN(ssn1.glmssn.EEE)
png('LOOCV_AR.png', width = 4, height = 4, units = 'in', res = 100)
#par(mfrow = c(1, 2))
plot(ssn1.glmssn.EEE$sampinfo$z,
     cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction", ylim = c(0,1))
abline(0, 1)
# plot( na.omit( getSSNdata.frame(ssn1)[, "FSS_26Aug14"]),
#       cv.out[, "cv.se"], pch = 19,
#       xlab = "Observed Data", ylab = "LOOCV Prediction SE")
dev.off()

#Unscaled and Untransformed
df.un <- cbind(cv.out, getSSNdata.frame(ssn1.glmssn.EEE)[,c('log10_FSS_26Aug14','FSS_26Aug14')])
#df.un$FSS_26Aug14 <- df.un$FSS_26Aug14*(min.max[min.max$variable == "FSS_26Aug14",'max_val']-min.max[min.max$variable == "FSS_26Aug14",'min_val']) + min.max[min.max$variable == "FSS_26Aug14",'min_val']
df.un$log10_FSS_26Aug14 <- ssn1.glmssn.EEE$sampinfo$z*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']
df.un$cv.pred <- df.un$cv.pred*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']
df.un$cv.se <- df.un$cv.se*(min.max[min.max$variable == "log10_FSS_26Aug14",'max_val']-min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']) + min.max[min.max$variable == "log10_FSS_26Aug14",'min_val']
#df.un <- as.data.frame(lapply(df.un, function(x) {10^x}))
df.un$cv.pred <- 10^df.un$cv.pred

png('LOOCV_Un.png', width = 4, height = 4, units = 'in', res = 100)
plot(10^df.un$log10_FSS_26Aug14, 10^df.un$cv.pred, pch = 19, xlab = "Observed Data", ylab = "LOOCV Prediction", ylim = c(0,70)) #ylim = c(0,1.8), xlim = c(0,1.8)
abline(0,1)
dev.off()
###################################################
### cross validation stats
###################################################
#We can use InfoCritCompare to evaluate model performance instead of calling
#CrossValidationStatsSSN for each individual model 
##ML##
variable.compare <- InfoCritCompare(list(ssn1.glmssn3.G,ssn1.glmssn4.G,ssn1.glmssn5.G,
                                         ssn1.glmssn6.G,ssn1.glmssn7.G,ssn1.glmssn8.G,ssn1.glmssn9.G,
                                         ssn1.glmssn10.G,ssn1.glmssn11.G,ssn1.glmssn12.G,ssn1.glmssn13.G,
                                         ssn1.glmssn14.G))
write.csv(variable.compare, 'variablecompare.csv')

##REML##
auto.compare <- InfoCritCompare(list(ssn1.glmssn.ELG,ssn1.glmssn.ELE,ssn1.glmssn.EME,ssn1.glmssn.ESE,
                                       ssn1.glmssn.SSE,ssn1.glmssn.LLE,ssn1.glmssn.MME,ssn1.glmssn.MEE,
                                       ssn1.glmssn.SSS,ssn1.glmssn.ESG,ssn1.glmssn.EEE))
write.csv(auto.compare, 'autocompare.csv')

###################################################
#### Untransformed RMSE for the selected model ####
###################################################
ssn1.glmssn.SSE.preds <- as.data.frame(getPreds(ssn1.glmssn.SSE))
ssn1.glmssn.SSE.preds$cv.pred.untran <- 10^(ssn1.glmssn.SSE.preds$cv.pred)
ssn1.glmssn.SSE.preds <- merge(ssn1.glmssn.SSE.preds, obs.vars[,c('pid','STATION_KEY','FSS_26Aug14')], by = 'pid')

library(hydroGOF)
rmse(ssn1.glmssn.SSE.preds$cv.pred.untran, ssn1.glmssn.SSE.preds$FSS_26Aug14)
#[1] 3.609358

ssn.ecomod.top <- ssn.ecomod.list[40][[1]]
ssn.ecomod.top.preds <- as.data.frame(getPreds(ssn.ecomod.top))
ssn.ecomod.top.preds$cv.pred.untran <- 10^(ssn.ecomod.top.preds$cv.pred)
ssn.ecomod.top.preds <- merge(ssn.ecomod.top.preds, obs.vars[,c('pid','STATION_KEY','FSS_26Aug14')], by = 'pid')

library(hydroGOF)
ssn.ecomod.top.preds$cv.pred.untran.unscale <- ssn.ecomod.top.preds$cv.pred.untran*(max(ssn.ecomod.top.preds$FSS_26Aug14)-min(ssn.ecomod.top.preds$FSS_26Aug14)) + min(ssn.ecomod.top.preds$FSS_26Aug14)
ssn.ecomod.top.preds$FSS_26Aug14.unscale <- ssn.ecomod.top.preds$FSS_26Aug14*(max(ssn.ecomod.top.preds$FSS_26Aug14)-min(ssn.ecomod.top.preds$FSS_26Aug14)) + min(ssn.ecomod.top.preds$FSS_26Aug14)
rmse(ssn.ecomod.top.preds$cv.pred.untran.unscale, ssn.ecomod.top.preds$FSS_26Aug14.unscale)
#[1] 3.610497

##################################
#### Check model significance ####
##################################
#Loop to make data.frame with X2, df and p-values for each model. To merge with InfoCritCompare
load("C:/Users/pbryant/Desktop/MidCoastTMDL-ModelRuns/Eco_Scaled_Transformed/allModelsList.Rdata")
load("ssn_SSE_NULL.Rdata")
for(i in 1:length(ssn.ecomod.list)) {
  output <- lrtSSN(ssn.SSE.NULL,ssn.ecomod.list[i][[1]])
  df <- data.frame(source = Reduce(paste, deparse(ssn.ecomod.list[i][[1]]$args$formula)), 
                   t(sapply(output,c)),
                   Cor.Models = paste(models[i][[1]]$args$CorModels,collapse=", "))
  ifelse(i == 1, ecomod.lrt <- df, ecomod.lrt <- rbind(ecomod.lrt, df))
}

rm(ssn.ecomod.list)
load("ssn1_glmssn3_G.Rdata")
load("ssn1_glmssn4_G.Rdata")
load("ssn1_glmssn5_G.Rdata")
load("ssn1_glmssn6_G.Rdata")
load("ssn1_glmssn7_G.Rdata")
load("ssn1_glmssn8_G.Rdata")
load("ssn1_glmssn9_G.Rdata")
load("ssn1_glmssn10_G.Rdata")
load("ssn1_glmssn11_G.Rdata")
load("ssn1_glmssn12_G.Rdata")
load("ssn1_glmssn13_G.Rdata")
load("ssn1_glmssn14_G.Rdata")
load("ssn1_glmssn2_NULL.Rdata")
models <- list(ssn1.glmssn3.G,ssn1.glmssn4.G,ssn1.glmssn5.G,
                                  ssn1.glmssn6.G,ssn1.glmssn7.G,ssn1.glmssn8.G,ssn1.glmssn9.G,
                                  ssn1.glmssn10.G,ssn1.glmssn11.G,ssn1.glmssn12.G,ssn1.glmssn13.G,
                                  ssn1.glmssn14.G)
for(i in 1:length(models)) {
  output <- lrtSSN(ssn1.glmssn2.NULL,models[i][[1]])
  df <- data.frame(source = Reduce(paste, deparse(models[i][[1]]$args$formula)), 
                   t(sapply(output,c)), 
                   Cor.Models = paste(models[i][[1]]$args$CorModels,collapse=", "))
  ifelse(i == 1, varsel.lrt <- df, varsel.lrt <- rbind(varsel.lrt, df))
}

qchisq(0.95, 5) #11.0705
#Based on this the p-value indicates significance and the X2 test statistic is greater than the critical value so we reject the null hypothesis that there is 
#no difference between the models

rm(list = ls()[grep('ssn1\\.|SSE',ls())])

#### Append lrt info to InfoCritCompare dfs ####
ecomod.lrt$source <- gsub("      "," ",ecomod.lrt$source)
ecomods <- read.csv('ecological_models_comparison.csv')

ecomods <- merge(ecomods, ecomod.lrt, by.x = 'formula', by.y = 'source')

varsel <- read.csv('variablecompare.csv')
varsel.lrt$source <- gsub("      "," ",varsel.lrt$source)
varsel.lrt$X <- 1:12
varsel <- merge(varsel, varsel.lrt, by = "X")
varsel <- within(varsel, rm(source,Cor.Models))


model.all <- rbind(ecomods, varsel)

View(arrange(model.all,AIC))

#Check aspatial models
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.aspatial <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days +  
                            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                            DAPOPRCA2010, 
                          ssn1,
                          EstMeth = "ML",
                          CorModels = c("locID"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

aspatial <- InfoCritCompare(list(ssn1.glmssn.aspatial))
#AIC is way higher than when the spatial autocovariance is accounted for

#Check top ecomod with same autocovariance as top mathematical model
start.time <- Sys.time()
print(start.time)
ssn1.top.ecomod.EEE <- glmssn(log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + PASILTRCA + POWNRCA_PRI + PDISRSA_1YR, 
                               ssn1,
                               EstMeth = "ML",
                               CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                               addfunccol = "afvArea",
                               family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

AIC(ssn1.top.ecomod.EEE)