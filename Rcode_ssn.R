#### SubSSN Test ####
library(SSN)
library(stringr)

options('scipen' = 100)

# ssn1 <- importSSN('bugs.ssn')
# obs <- getSSNdata.frame(ssn1, Name = "Obs")
# obs <- rename(obs, c(
#   "STATION_KE" = "STATION_KEY",
#   "APOPRCA201" = "APOPRCA2010",
#   "sum_1095_d" = "sum_1095_days",
#   "FSS_26Aug1" = "FSS_26Aug14",
#   "PALITHEROD" = "PALITHERODRCA",
#   "PADISRSA_1" = "PADISRSA_1YR",
#   "PASUSCEP5_" = "PASUSCEP5_DE", 
#   "POWNRCA_PR" = "POWNRCA_PRI",
#   "POWNRCA_FE" = "POWNRCA_FED",
#   "PAOWNRCA_A" = "PAOWNRCA_AGR",
#   "log10_FSS_" = "log10_FSS_26Aug14",
#   "log10_sum_" = "log10_sum_1095_days",
#   "sqrt_PADIS" = "sqrt_PADISRSA_1YR",
#   "bin_PALITH" = "bin_PALITHERODRCA",
#   "log10_XSLO" = "log10_XSLOPE_MAP",
#   "log10_PASI" = "log10_PASILTRCA",
#   "log10_MIN_" = "log10_MIN_Z",
#   "log10_STRM" = "log10_STRMPWR",
#   "sqrt_upDis" = "sqrt_upDist",
#   "log10_APOP" = "log10_APOPRCA2010",
#   "log10_PASU" = "log10_PASUSCEP5_DE",
#   "log10_POWN" = "log10_POWNRCA_PRI",
#   "log10_POWN.1" = "log10_POWNRCA_FED",
#   "bin_PAOWNR" = "bin_PAOWNRCA_AGR"))
# ssn1 <- putSSNdata.frame(obs, ssn1, Name = 'Obs')

# #save the ssn object
# writeSSN(ssn1, filename = 'C:/users/pbryant/desktop/SubSSN/LSN/lsn.ssn')

##################################################
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

# #create data frame to store model scores
# comp <- data.frame(model = character(), family = character(),GR2 = numeric(), AIC = numeric())
# model.name <- paste('ssn1.glmssn',i,sep='')
# comp <- rbind(comp, data.frame(model = model.name, family = model$args$family, GR2 = GR2(model), AIC = AIC(model)))

#### ssn1.glmssn1.P ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1.P <- glmssn(FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                         XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                         APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                         PAOWNRCA_AGR, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                       addfunccol = "afvArea",
                       family = "Poisson")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1.P)

#### ssn1.glmssn1.P Summary ####
# Call:
#   glmssn(formula = FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Poisson", CorModels = c("locID", 
#                                                                               "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# 7667  23712  30235  41034 397286 
# 
# Coefficients:
#   Estimate       Std. Error  t value            Pr(>|t|)    
# (Intercept)   -32348.079323716      1.515386808 -21346.4 <0.0000000000000002 ***
#   sum_1095_days     -3.119775127      0.000140869 -22146.6 <0.0000000000000002 ***
#   PALITHERODRCA    206.273919514      0.008752608  23567.1 <0.0000000000000002 ***
#   PADISRSA_1YR     163.218507408      0.012698045  12853.8 <0.0000000000000002 ***
#   XSLOPE_MAP       -84.184512283      0.026917236  -3127.5 <0.0000000000000002 ***
#   PASILTRCA        191.992964442      0.017466533  10992.0 <0.0000000000000002 ***
#   MIN_Z             -4.945760512      0.000523828  -9441.6 <0.0000000000000002 ***
#   upDist            -0.012169281      0.000001762  -6905.4 <0.0000000000000002 ***
#   APOPRCA2010        9.823244657      0.000660344  14876.0 <0.0000000000000002 ***
#   PASUSCEP5_DE    -103.129967646      0.034064173  -3027.5 <0.0000000000000002 ***
#   POWNRCA_FED      -59.623742572      0.004959222 -12022.8 <0.0000000000000002 ***
#   POWNRCA_PRI        1.501873168      0.003594141    417.9 <0.0000000000000002 ***
#   PAOWNRCA_AGR    -185.596252579      0.046173136  -4019.6 <0.0000000000000002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter      Estimate
# Exponential.tailup   parsill       0.00791
# Exponential.tailup     range       0.94441
# Exponential.taildown   parsill       5.82149
# Exponential.taildown     range    9372.63034
# Exponential.Euclid   parsill      31.43260
# Exponential.Euclid     range 1855839.15414
# locID   parsill       0.14508
# Nugget   parsill       7.71346
# 
# Residual standard error: 6.717183
# Generalized R-squared: 0.2012738
# 
# varcomp(ssn1.glmssn1.P)
# VarComp   Proportion
# 1    Covariates (R-sq) 0.2012738265
# 2   Exponential.tailup 0.0001400769
# 3 Exponential.taildown 0.1030523249
# 4   Exponential.Euclid 0.5564215456
# 5                locID 0.0025682282
# 6               Nugget 0.1365439979

#### ssn1,glmssn2.P ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn2.P <- glmssn(FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "REML",
                         CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Spherical.Euclid"),
                         addfunccol = "afvArea",
                         family = "Poisson")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn2.P)
varcomp(ssn1.glmssn2.P)

sum(2*ssn1.glmssn2.P$sampinfo$rankX,
    2*length(covparms(ssn1.glmssn2.P)[,3]),
    ssn1.glmssn2.P$estimate$m2LL)


#### ssn1.glmssn2.P summary ####
# Call:
#   glmssn(formula = FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Poisson", CorModels = c("locID", 
#                                                                               "Spherical.tailup", "Spherical.taildown", "Spherical.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.8348 -0.3455 -0.1097  0.2799  5.2519 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    2.0937461137  0.3027896272   6.915 < 0.0000000000000002 ***
#   sum_1095_days -0.0001426263  0.0000205638  -6.936 < 0.0000000000000002 ***
#   PALITHERODRCA  0.0085079981  0.0012093449   7.035 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.0113573673  0.0035260561   3.221              0.00133 ** 
#   XSLOPE_MAP    -0.0076710420  0.0070249368  -1.092              0.27519    
# PASILTRCA      0.0155914644  0.0045029632   3.462              0.00056 ***
#   MIN_Z         -0.0001263228  0.0001146507  -1.102              0.27089    
# upDist        -0.0000005927  0.0000004324  -1.371              0.17087    
# APOPRCA2010    0.0005734273  0.0001790293   3.203              0.00142 ** 
#   PASUSCEP5_DE  -0.0018607305  0.0082866062  -0.225              0.82239    
# POWNRCA_FED   -0.0024792771  0.0011756569  -2.109              0.03528 *  
#   POWNRCA_PRI    0.0002202237  0.0010075880   0.219              0.82705    
# PAOWNRCA_AGR  -0.0074959239  0.0145567020  -0.515              0.60674    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter      Estimate
# Spherical.tailup   parsill     0.3531532
# Spherical.tailup     range   719.5558850
# Spherical.taildown   parsill     1.2234418
# Spherical.taildown     range  8475.0751669
# Spherical.Euclid   parsill     0.5859025
# Spherical.Euclid     range 89735.7125985
# locID   parsill     0.0000177
# Nugget   parsill     1.9439507
# 
# Residual standard error: 2.026442
# Generalized R-squared: 0.2031228
# > varcomp(ssn1.glmssn2.P)
# VarComp     Proportion
# 1  Covariates (R-sq) 0.203122753305
# 2   Spherical.tailup 0.068530887282
# 3 Spherical.taildown 0.237414105161
# 4   Spherical.Euclid 0.113696880150
# 5              locID 0.000003429584
# 6             Nugget 0.377231944519
# 
AIC(ssn1.glmssn2.P)

#### ssn1.glmssn1.G ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1.G <- glmssn(FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + STRMPWR + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR, 
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

save(ssn1.glmssn1.G, file = 'ssn1_glmssn1_G.Rdata')



#### ssn1.glmssn1.G Summary ####
# Call:
#   glmssn(formula = FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + STRMPWR + 
#            upDist + APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# NA -4.733 -1.883  1.804     NA 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    4.866772280  4.757936414   1.023              0.30669    
# sum_1095_days -0.001015342  0.000223751  -4.538              0.00001 ***
#   PALITHERODRCA  0.080914768  0.014885451   5.436 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.179376638  0.058489303   3.067              0.00224 ** 
#   XSLOPE_MAP    -0.071579399  0.096926333  -0.738              0.46044    
# PASILTRCA      0.224054968  0.080109363   2.797              0.00529 ** 
#   MIN_Z          0.000216243  0.001416546   0.153              0.87871    
# STRMPWR       -0.002311689  0.001402346  -1.648              0.09967 .  
# upDist        -0.000013174  0.000007087  -1.859              0.06343 .  
# APOPRCA2010    0.010376813  0.003320259   3.125              0.00184 ** 
#   PASUSCEP5_DE   0.032336816  0.104520168   0.309              0.75711    
# POWNRCA_FED   -0.029057149  0.016818600  -1.728              0.08445 .  
# POWNRCA_PRI    0.006150014  0.015159138   0.406              0.68508    
# PAOWNRCA_AGR   0.231576609  0.273762558   0.846              0.39787    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter      Estimate
# Exponential.tailup   parsill       0.00780
# Exponential.tailup     range 1207587.53828
# Exponential.taildown   parsill      22.35247
# Exponential.taildown     range    7470.36495
# Exponential.Euclid   parsill      23.72528
# Exponential.Euclid     range   84142.08191
# locID   parsill       0.00132
# Nugget   parsill      25.83852
# 
# Residual standard error: 8.480884
# Generalized R-squared: 0.1739565
# > varcomp(ssn1.glmssn1.G)
# VarComp    Proportion
# 1    Covariates (R-sq) 0.17395649382
# 2   Exponential.tailup 0.00008957718
# 3 Exponential.taildown 0.25671204844
# 4   Exponential.Euclid 0.27247831654
# 5                locID 0.00001519753
# 6               Nugget 0.29674836650

#### ssn1.glmssn2.G ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn2.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "REML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn2.G)
varcomp(ssn1.glmssn2.G)
#Plot the residuals
ssn1.resid2.G <- residuals(ssn1.glmssn2.G)
par(mfrow = c(1,2))
hist(ssn1.resid2.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid2.G)

resids2.df <- getSSNdata.frame(ssn1.resid2.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn2.G, file = 'ssn1_glmssn2_G.Rdata')

AIC(ssn1.glmssn2.G)
#-105.5518

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.73635 -0.12893  0.01022  0.15388  0.79842 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.9213661205  0.1607734132   5.731 < 0.0000000000000002 ***
#   sum_1095_days -0.0000526756  0.0000062505  -8.427 < 0.0000000000000002 ***
#   PALITHERODRCA  0.0024591410  0.0004020573   6.116 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.0040278403  0.0016621312   2.423              0.01561 *  
#   XSLOPE_MAP    -0.0000178327  0.0025361037  -0.007              0.99439    
# PASILTRCA      0.0059055871  0.0021891557   2.698              0.00714 ** 
#   MIN_Z         -0.0000268923  0.0000399062  -0.674              0.50059    
# upDist        -0.0000002423  0.0000001889  -1.283              0.19990    
# APOPRCA2010    0.0002085302  0.0000944620   2.208              0.02757 *  
#   PASUSCEP5_DE  -0.0015281608  0.0027803510  -0.550              0.58273    
# POWNRCA_FED   -0.0014553292  0.0004525459  -3.216              0.00135 ** 
#   POWNRCA_PRI    0.0001767840  0.0004128116   0.428              0.66859    
# PAOWNRCA_AGR   0.0039318593  0.0078508560   0.501              0.61664    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter      Estimate
# Exponential.tailup   parsill      0.000767
# Exponential.tailup     range  47485.747209
# Exponential.taildown   parsill      0.016826
# Exponential.taildown     range  11182.487976
# Exponential.Euclid   parsill      0.028016
# Exponential.Euclid     range 410423.667504
# locID   parsill      0.005254
# Nugget   parsill      0.020986
# 
# Residual standard error: 0.2680474
# Generalized R-squared: 0.254094
# 
# VarComp  Proportion
# 1    Covariates (R-sq) 0.254093965
# 2   Exponential.tailup 0.007966887
# 3 Exponential.taildown 0.174682084
# 4   Exponential.Euclid 0.290844619
# 5                locID 0.054547087
# 6               Nugget 0.217865359

#### ssn1.glmssn3.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#This is the first model that incorporates all the variables from the randomForest selection. The others have just been for testing up to this point.

start.time <- Sys.time()
print(start.time)
ssn1.glmssn3.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR + STRMPWR + LAT_RAW, 
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
# -261.2114

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR + STRMPWR + LAT_RAW, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.731794 -0.129756  0.009524  0.156384  0.800008 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.8414788007  1.7733305450   0.475              0.63526    
# sum_1095_days -0.0000522965  0.0000061853  -8.455 < 0.0000000000000002 ***
#   PALITHERODRCA  0.0024762179  0.0004005154   6.183 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.0041019505  0.0016548215   2.479              0.01340 *  
#   XSLOPE_MAP     0.0011391689  0.0026112807   0.436              0.66278    
# PASILTRCA      0.0056640424  0.0022155014   2.557              0.01076 *  
#   MIN_Z         -0.0000261882  0.0000396463  -0.661              0.50910    
# upDist        -0.0000002584  0.0000001849  -1.398              0.16253    
# APOPRCA2010    0.0002264328  0.0000943252   2.401              0.01661 *  
#   PASUSCEP5_DE  -0.0011060781  0.0027839406  -0.397              0.69125    
# POWNRCA_FED   -0.0014691880  0.0004539825  -3.236              0.00126 ** 
#   POWNRCA_PRI    0.0002151096  0.0004131586   0.521              0.60276    
# PAOWNRCA_AGR   0.0042068033  0.0077824816   0.541              0.58898    
# STRMPWR       -0.0000623496  0.0000399368  -1.561              0.11889    
# LAT_RAW        0.0019786470  0.0403950199   0.049              0.96095    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter     Estimate
# Exponential.tailup   parsill      0.00690
# Exponential.tailup     range    702.33004
# Exponential.taildown   parsill      0.01487
# Exponential.taildown     range  14305.21401
# Exponential.Euclid   parsill      0.01473
# Exponential.Euclid     range 188442.85614
# locID   parsill      0.00117
# Nugget   parsill      0.02013
# 
# Residual standard error: 0.2404168
# Generalized R-squared: 0.2566976
# > varcomp(ssn1.glmssn3.G)
# VarComp Proportion
# 1    Covariates (R-sq) 0.25669763
# 2   Exponential.tailup 0.08871207
# 3 Exponential.taildown 0.19121298
# 4   Exponential.Euclid 0.18941116
# 5                locID 0.01504784
# 6               Nugget 0.25891831

#### ssn1.glmssn4.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
start.time <- Sys.time()
print(start.time)
ssn1.glmssn4.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
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
ssn1.resid4.G <- residuals(ssn1.glmssn4.G)
par(mfrow = c(1,2))
hist(ssn1.resid4.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid4.G)

resids2.df <- getSSNdata.frame(ssn1.resid4.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn4.G, file = 'ssn1_glmssn4_G.Rdata')

AIC(ssn1.glmssn4.G)
# -271.8111
#This one has the best AIC. 

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + PASILTRCA + APOPRCA2010 + POWNRCA_FED, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.76655 -0.13678  0.01146  0.15282  0.76962 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    0.870171365  0.117844516   7.384 < 0.0000000000000002 ***
#   sum_1095_days -0.000053320  0.000006051  -8.811 < 0.0000000000000002 ***
#   PALITHERODRCA  0.002498675  0.000387002   6.456 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.004526428  0.001524776   2.969              0.00308 ** 
#   PASILTRCA      0.006511561  0.002039054   3.193              0.00146 ** 
#   APOPRCA2010    0.000226396  0.000089475   2.530              0.01159 *  
#   POWNRCA_FED   -0.001632472  0.000306690  -5.323 < 0.0000000000000002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter     Estimate
# Exponential.tailup   parsill      0.00322
# Exponential.tailup     range    529.99717
# Exponential.taildown   parsill      0.01647
# Exponential.taildown     range  11294.57447
# Exponential.Euclid   parsill      0.01397
# Exponential.Euclid     range 158323.16554
# locID   parsill      0.00201
# Nugget   parsill      0.02057
# 
# Residual standard error: 0.2371334
# Generalized R-squared: 0.2513855
# > varcomp(ssn1.glmssn4.G)
# VarComp Proportion
# 1    Covariates (R-sq) 0.25138548
# 2   Exponential.tailup 0.04286850
# 3 Exponential.taildown 0.21922058
# 4   Exponential.Euclid 0.18595572
# 5                locID 0.02678667
# 6               Nugget 0.27378306

#### ssn1.glmssn5.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
start.time <- Sys.time()
print(start.time)
ssn1.glmssn5.G <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           PASILTRCA + POWNRCA_FED, 
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
ssn1.resid5.G <- residuals(ssn1.glmssn5.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn5.G, file = 'ssn1_glmssn5_G.Rdata')

AIC(ssn1.glmssn5.G)
# -267.5148
#There does not appear to be a benefit to removing the population variable.

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + PASILTRCA + POWNRCA_FED, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.76878 -0.13597  0.01129  0.15310  0.77415 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)    0.87680893  0.11874820   7.384 < 0.0000000000000002 ***
#   sum_1095_days -0.00005369  0.00000606  -8.860 < 0.0000000000000002 ***
#   PALITHERODRCA  0.00247818  0.00038756   6.394 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.00528970  0.00150253   3.521              0.00046 ***
#   PASILTRCA      0.00652469  0.00204348   3.193              0.00147 ** 
#   POWNRCA_FED   -0.00169298  0.00030692  -5.516 < 0.0000000000000002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00574723
# Exponential.tailup     range    473.72263077
# Exponential.taildown   parsill      0.01644180
# Exponential.taildown     range  11326.95610124
# Exponential.Euclid   parsill      0.01411386
# Exponential.Euclid     range 167268.47136089
# locID   parsill      0.00000225
# Nugget   parsill      0.02057144
# 
# Residual standard error: 0.2384881
# Generalized R-squared: 0.2462756
# > varcomp(ssn1.glmssn5.G)
# VarComp    Proportion
# 1    Covariates (R-sq) 0.24627557668
# 2   Exponential.tailup 0.07616186937
# 3 Exponential.taildown 0.21788561517
# 4   Exponential.Euclid 0.18703584362
# 5                locID 0.00002978427
# 6               Nugget 0.27261131089

#### ssn1.glmssn1.S ####
#Now let's see if the autocovariance functions affect the AIC
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1.S <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR + STRMPWR + LAT_RAW, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1.S)
varcomp(ssn1.glmssn1.S)
#Plot the residuals
ssn1.resid1.S <- residuals(ssn1.glmssn1.S)
par(mfrow = c(1,2))
hist(ssn1.resid1.S)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid1.S)

resids1.S.df <- getSSNdata.frame(ssn1.resid1.S)
plot(resids1.S.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn1.S, file = 'ssn1_glmssn1_S.Rdata')

AIC(ssn1.glmssn1.S)
#-262.8767

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR + STRMPWR + LAT_RAW, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Spherical.tailup", "Spherical.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.725211 -0.130364  0.006165  0.151989  0.792216 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.9981292730  1.4574090835   0.685              0.49364    
# sum_1095_days -0.0000537524  0.0000061249  -8.776 < 0.0000000000000002 ***
#   PALITHERODRCA  0.0025025748  0.0003938098   6.355 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.0040887516  0.0016403509   2.493              0.01289 *  
#   XSLOPE_MAP     0.0011789909  0.0025779013   0.457              0.64755    
# PASILTRCA      0.0055523222  0.0021809832   2.546              0.01110 *  
#   MIN_Z         -0.0000303168  0.0000390525  -0.776              0.43781    
# upDist        -0.0000002621  0.0000001788  -1.466              0.14316    
# APOPRCA2010    0.0002292405  0.0000931712   2.460              0.01410 *  
#   PASUSCEP5_DE  -0.0011220226  0.0027471353  -0.408              0.68307    
# POWNRCA_FED   -0.0014504391  0.0004480669  -3.237              0.00126 ** 
#   POWNRCA_PRI    0.0001988436  0.0004082062   0.487              0.62632    
# PAOWNRCA_AGR   0.0041264565  0.0077011564   0.536              0.59224    
# STRMPWR       -0.0000626079  0.0000396136  -1.580              0.11441    
# LAT_RAW       -0.0012275076  0.0333004390  -0.037              0.97060    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Spherical.tailup   parsill      0.0063521
# Spherical.tailup     range    321.2121182
# Spherical.taildown   parsill      0.0148872
# Spherical.taildown     range   9952.7987863
# Exponential.Euclid   parsill      0.0114714
# Exponential.Euclid     range 132084.4250195
# locID   parsill      0.0000554
# Nugget   parsill      0.0206667
# 
# Residual standard error: 0.2311555
# Generalized R-squared: 0.2640405
# > varcomp(ssn1.glmssn1.S)
# VarComp  Proportion
# 1  Covariates (R-sq) 0.264040526
# 2   Spherical.tailup 0.087490826
# 3 Spherical.taildown 0.205048872
# 4 Exponential.Euclid 0.158002230
# 5              locID 0.000763471
# 6             Nugget 0.284654075

#### ssn1.glmssn2.S ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn2.S <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                        ssn1,
                        EstMeth = "ML",
                        CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                        addfunccol = "afvArea",
                        family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn2.S)
varcomp(ssn1.glmssn2.S)
#Plot the residuals
ssn1.resid2.S <- residuals(ssn1.glmssn2.S)
par(mfrow = c(1,2))
hist(ssn1.resid2.S)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid2.S)

resids2.S.df <- getSSNdata.frame(ssn1.resid2.S)
plot(resids2.S.df[,"_fit_"],resids2.S.df[,"_resid_"])

save(ssn1.glmssn2.S, file = 'ssn1_glmssn2_S.Rdata')

AIC(ssn1.glmssn2.S)
#-272.7678

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + PASILTRCA + APOPRCA2010 + POWNRCA_FED, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Spherical.tailup", 
#                                             "Spherical.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.76678 -0.13633  0.01097  0.15277  0.77002 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    0.869362547  0.117601114   7.392 < 0.0000000000000002 ***
#   sum_1095_days -0.000053574  0.000006045  -8.863 < 0.0000000000000002 ***
#   PALITHERODRCA  0.002496616  0.000386475   6.460 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.004517868  0.001524345   2.964              0.00313 ** 
#   PASILTRCA      0.006553825  0.002037618   3.216              0.00135 ** 
#   APOPRCA2010    0.000227681  0.000089560   2.542              0.01121 *  
#   POWNRCA_FED   -0.001618030  0.000306630  -5.277 < 0.0000000000000002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Spherical.tailup   parsill      0.00713458
# Spherical.tailup     range    686.18254654
# Spherical.taildown   parsill      0.01459922
# Spherical.taildown     range   9978.80070637
# Exponential.Euclid   parsill      0.01380916
# Exponential.Euclid     range 157825.45440611
# locID   parsill      0.00000931
# Nugget   parsill      0.02055161
# 
# Residual standard error: 0.2368626
# Generalized R-squared: 0.2520574
# > varcomp(ssn1.glmssn2.S)
# VarComp   Proportion
# 1  Covariates (R-sq) 0.2520574199
# 2   Spherical.tailup 0.0951138505
# 3 Spherical.taildown 0.1946278590
# 4 Exponential.Euclid 0.1840952780
# 5              locID 0.0001241361
# 6             Nugget 0.2739814565

#### ssn1.glmssn1.ELE ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#For the fixed autocovariance functions in the first step they used linear with sill tail down to maximize the autocorrelation between 
#flow un-connected sites and relied on the exponential tail up to minimize autocorrelation between unconnected sites.
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1.ELE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR + STRMPWR + LAT_RAW, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID","Exponential.tailup", "LinearSill.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1.ELE)
varcomp(ssn1.glmssn1.ELE)
#Plot the residuals
ssn1.resid1.ELE <- residuals(ssn1.glmssn1.ELE)
par(mfrow = c(1,2))
hist(ssn1.resid1.ELE)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid1.ELE)

resids.df <- getSSNdata.frame(ssn1.resid1.ELE)
plot(resids.df[,"_fit_"],resids.df[,"_resid_"])

save(ssn1.glmssn1.ELE, file = 'ssn1_glmssn1_ELE.Rdata')

AIC(ssn1.glmssn1.ELE)
#-262.3387

#AIC is not improved by use of a different autocovariance model selection. Noted as well, the significant variables
#are consistently the same between each of the model compositions.

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR + STRMPWR + LAT_RAW, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "LinearSill.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.724771 -0.129768  0.004191  0.150793  0.793053 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.9894855249  1.4478411185   0.683              0.49455    
# sum_1095_days -0.0000540864  0.0000061291  -8.825 < 0.0000000000000002 ***
#   PALITHERODRCA  0.0025053108  0.0003931128   6.373 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.0040321547  0.0016410912   2.457              0.01423 *  
#   XSLOPE_MAP     0.0011874494  0.0025755882   0.461              0.64490    
# PASILTRCA      0.0054683287  0.0021807159   2.508              0.01236 *  
#   MIN_Z         -0.0000311348  0.0000390724  -0.797              0.42578    
# upDist        -0.0000002745  0.0000001789  -1.535              0.12531    
# APOPRCA2010    0.0002276656  0.0000931940   2.443              0.01479 *  
#   PASUSCEP5_DE  -0.0011238847  0.0027394602  -0.410              0.68173    
# POWNRCA_FED   -0.0014221902  0.0004482975  -3.172              0.00157 ** 
#   POWNRCA_PRI    0.0002002506  0.0004079753   0.491              0.62368    
# PAOWNRCA_AGR   0.0042517540  0.0077069204   0.552              0.58133    
# STRMPWR       -0.0000609293  0.0000396698  -1.536              0.12497    
# LAT_RAW       -0.0008907747  0.0330859976  -0.027              0.97853    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0069953
# Exponential.tailup     range    559.3789945
# LinearSill.taildown   parsill      0.0141920
# LinearSill.taildown     range   6894.5681549
# Exponential.Euclid   parsill      0.0115075
# Exponential.Euclid     range 128817.6079221
# locID   parsill      0.0000299
# Nugget   parsill      0.0207015
# 
# Residual standard error: 0.2311411
# Generalized R-squared: 0.2633196
# > varcomp(ssn1.glmssn1.ELE)
# VarComp   Proportion
# 1   Covariates (R-sq) 0.2633196182
# 2  Exponential.tailup 0.0964563339
# 3 LinearSill.taildown 0.1956899400
# 4  Exponential.Euclid 0.1586735023
# 5               locID 0.0004126171
# 6              Nugget 0.2854479886


#### ssn1.glmssn6.G ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn6.G <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + bin_PALITHERODRCA + sqrt_PADISRSA_1YR + 
                         log10_PASILTRCA + log10_APOPRCA2010 + log10_POWNRCA_FED, 
                       ssn1,
                       EstMeth = "ML",
                       CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                       addfunccol = "afvArea",
                       family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn6.G)

summary(ssn1.glmssn1.ELE)
varcomp(ssn1.glmssn1.ELE)
#Plot the residuals
ssn1.resid1.ELE <- residuals(ssn1.glmssn1.ELE)
par(mfrow = c(1,2))
hist(ssn1.resid1.ELE)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid1.ELE)

resids.df <- getSSNdata.frame(ssn1.resid1.ELE)
plot(resids.df[,"_fit_"],resids.df[,"_resid_"])

save(ssn1.glmssn6.G, file = 'ssn1_glmssn6_G.Rdata')

AIC(ssn1.glmssn6.G)
#[1] -256.9108

#The transformation of the explanatory variables does not improve AIC or RMSPE (based on ML)

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ log10_sum_1095_days + bin_PALITHERODRCA + 
#            sqrt_PADISRSA_1YR + log10_PASILTRCA + log10_APOPRCA2010 + 
#            log10_POWNRCA_FED, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.725649 -0.134722  0.005367  0.162232  0.771079 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)          3.361267   0.568091   5.917 < 0.0000000000000002 ***
#   log10_sum_1095_days -0.910710   0.101990  -8.929 < 0.0000000000000002 ***
#   bin_PALITHERODRCA    0.154180   0.027764   5.553 < 0.0000000000000002 ***
#   sqrt_PADISRSA_1YR    0.012671   0.008282   1.530              0.12640    
# log10_PASILTRCA      0.592717   0.229464   2.583              0.00998 ** 
#   log10_APOPRCA2010    0.051196   0.014505   3.529              0.00044 ***
#   log10_POWNRCA_FED   -0.066784   0.015690  -4.256              0.00002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00560538
# Exponential.tailup     range    511.51611462
# Exponential.taildown   parsill      0.01600455
# Exponential.taildown     range  10886.22083255
# Exponential.Euclid   parsill      0.01721961
# Exponential.Euclid     range 124893.93346665
# locID   parsill      0.00000260
# Nugget   parsill      0.02039891
# 
# Residual standard error: 0.2433743
# Generalized R-squared: 0.2218016

#### ssn1.glmssn1.NULL ####
#Check the aspatial model
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1.NULL <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                         ssn1,
                         EstMeth = "ML",
                         CorModels = c("locID"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1.NULL)
varcomp(ssn1.glmssn1.NULL)
#Plot the residuals
ssn1.resid1.NULL <- residuals(ssn1.glmssn1.NULL)
par(mfrow = c(1,2))
hist(ssn1.resid1.NULL)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid1.NULL)

resids2.df <- getSSNdata.frame(ssn1.resid1.NULL)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

save(ssn1.glmssn1.NULL, file = 'ssn1_glmssn1_NULL.Rdata')

AIC(ssn1.glmssn1.NULL)
#-206.4454

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + PASILTRCA + APOPRCA2010 + POWNRCA_FED, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.8003182 -0.1311341  0.0005135  0.1506414  0.7969014 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    1.018177136  0.078115465  13.034 < 0.0000000000000002 ***
#   sum_1095_days -0.000067292  0.000004802 -14.014 < 0.0000000000000002 ***
#   PALITHERODRCA  0.003113405  0.000300904  10.347 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.005175383  0.001496490   3.458              0.00057 ***
#   PASILTRCA      0.003998552  0.001320541   3.028              0.00254 ** 
#   APOPRCA2010    0.000247019  0.000089454   2.761              0.00589 ** 
#   POWNRCA_FED   -0.001424182  0.000249365  -5.711 < 0.0000000000000002 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter Estimate
# locID   parsill   0.0279
# Nugget   parsill   0.0223
# 
# Residual standard error: 0.2241273
# Generalized R-squared: 0.4444287
# > varcomp(ssn1.glmssn1.NULL)
# VarComp Proportion
# 1 Covariates (R-sq)  0.4444287
# 2             locID  0.3086908
# 3            Nugget  0.2468805

#### ssn1.glmssn2.NULL ####
#Try a null model
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
ssn1.resid2.NULL <- residuals(ssn1.glmssn2.NULL)
par(mfrow = c(1,2))
hist(ssn1.resid2.NULL)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid2.NULL)

resids2.df <- getSSNdata.frame(ssn1.resid2.NULL)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])

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

#### ssn1.glmssn.EEE ####
1st.start.time <- Sys.time()
print(1st.start.time)
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.EEE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                         ssn1,
                         EstMeth = "REML",
                         CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.ELG ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ELG <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "LinearSill.taildown", "Gaussian.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.ELE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ELE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "LinearSill.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.EME ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.EME <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Mariah.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.ESE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ESE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Spherical.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.SSE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.SSE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

save(ssn1.glmssn.SSE, file = 'ssn1_glmssn_SSE.Rdata')

#### ssn1.glmssn.LLE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.LLE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","LinearSill.tailup", "LinearSill.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.MME ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.MME <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Mariah.tailup", "Mariah.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.MEE ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.MEE <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Mariah.tailup", "Exponential.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.SSS ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.SSS <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Spherical.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.ESG ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.ESG <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                            PASILTRCA + APOPRCA2010 + POWNRCA_FED, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Exponential.tailup", "Spherical.taildown", "Gaussian.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.SSE.2 ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.SSE.2 <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PADISRSA_1YR, 
                          ssn1,
                          EstMeth = "REML",
                          CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                          addfunccol = "afvArea",
                          family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

#### ssn1.glmssn.SSE.3 ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.SSE.3 <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                              XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                              APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                              PAOWNRCA_AGR + STRMPWR + LAT_RAW, 
                            ssn1,
                            EstMeth = "REML",
                            CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                            addfunccol = "afvArea",
                            family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

###################################################
### check the residuals
###################################################
ssn1.resid1 <- residuals(ssn1.glmssn.SSE)
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
cv.out <- CrossValidationSSN(ssn1.glmssn.SSE)
par(mfrow = c(1, 2))
plot(ssn1.glmssn.SSE$sampinfo$z,
     cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(ssn1)[, "FSS_26Aug14"]),
      cv.out[, "cv.se"], pch = 19,
      xlab = "Observed Data", ylab = "LOOCV Prediction SE")

###################################################
### cross validation stats
###################################################
##ML##
CrossValidationStatsSSN(ssn1.glmssn4.G)
#           bias    std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001191478 -0.00201606 0.1883938 0.1913374 0.9864037 0.8237548 0.9131545 0.9438059

CrossValidationStatsSSN(ssn1.glmssn5.G)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001095451 -0.001863036 0.1890676 0.1919526 0.9866342 0.8224777 0.9118774 0.9438059

CrossValidationStatsSSN(ssn1.glmssn6.G)
#           bias     std.bias   RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001430853 -0.002417528 0.18973 0.1924555 0.9877347 0.8339719 0.9118774 0.9386973

CrossValidationStatsSSN(ssn1.glmssn1.S)
#           bias     std.bias   RMSPE      RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001538842 -0.002422118 0.19059 0.191696 0.9978387 0.8250319 0.9106003 0.9438059

CrossValidationStatsSSN(ssn1.glmssn2.S)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90   cov.95
# 1 -0.001123455 -0.001900468 0.1883655 0.1910589 0.9880666 0.8224777 0.9157088 0.945083

CrossValidationStatsSSN(ssn1.glmssn1.ELE)
#           bias    std.bias    RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001388988 -0.00217256 0.190798 0.1917153 0.9988858 0.8263091 0.9106003 0.9463602


CrossValidationStatsSSN(ssn1.glmssn1.NULL)
#          bias      std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90   cov.95
# 1 -0.00042555 -0.0007098064 0.2048384 0.2056859 0.9988453 0.8186462 0.9118774 0.945083


##REML##
CrossValidationStatsSSN(ssn1.glmssn.EEE)
#           bias    std.bias    RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001219661 -0.00206018 0.188277 0.1918349 0.9811075 0.8250319 0.9131545 0.9438059

CrossValidationStatsSSN(ssn1.glmssn.ELE)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001048184 -0.001769277 0.1884244 0.1915877 0.9833104 0.8237548 0.9157088 0.9463602

CrossValidationStatsSSN(ssn1.glmssn.ELG)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80   cov.90    cov.95
# 1 -0.000918074 -0.001548822 0.1891589 0.1921395 0.9838551 0.8275862 0.916986 0.9463602

CrossValidationStatsSSN(ssn1.glmssn.EME)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90   cov.95
# 1 -0.001248669 -0.002107834 0.1883868 0.1920776 0.9803915 0.8275862 0.9118774 0.945083

CrossValidationStatsSSN(ssn1.glmssn.ESE)
#           bias    std.bias     RMSPE       RAV  std.MSPE    cov.80   cov.90   cov.95
# 1 -0.001173164 -0.00198104 0.1882641 0.1916325 0.9820012 0.8250319 0.916986 0.945083

CrossValidationStatsSSN(ssn1.glmssn.LLE)
#          bias    std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.00105426 -0.00178155 0.1885075 0.1916772 0.9828379 0.8237548 0.9157088 0.9463602

CrossValidationStatsSSN(ssn1.glmssn.MEE)
#         bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.0012244 -0.002067542 0.1883018 0.1918595 0.9811502 0.8263091 0.9144317 0.9438059

CrossValidationStatsSSN(ssn1.glmssn.MME)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90   cov.95
# 1 -0.001248691 -0.002107895 0.1883866 0.1920746 0.9804132 0.8275862 0.9118774 0.945083

CrossValidationStatsSSN(ssn1.glmssn.SSE)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80   cov.90   cov.95
# 1 -0.001182433 -0.001996497 0.1882473 0.1916278 0.9819218 0.8250319 0.916986 0.945083

CrossValidationStatsSSN(ssn1.glmssn.SSS)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80   cov.90   cov.95
# 1 -0.001169532 -0.001975525 0.1883464 0.1916555 0.9823161 0.8237548 0.916986 0.945083

CrossValidationStatsSSN(ssn1.glmssn.ESG)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001055536 -0.001779893 0.1888768 0.1920662 0.9829735 0.8250319 0.9182631 0.9463602

GR2(ssn1.glmssn.EEE) # 0.2447016
GR2(ssn1.glmssn.ELE) # 0.2450364
GR2(ssn1.glmssn.ELG) # 0.2645521
GR2(ssn1.glmssn.EME) # 0.2405389
GR2(ssn1.glmssn.ESE) # 0.2454929
GR2(ssn1.glmssn.LLE) # 0.2451181
GR2(ssn1.glmssn.MEE) # 0.244801
GR2(ssn1.glmssn.MME) # 0.240545
GR2(ssn1.glmssn.SSE) # 0.2457414
GR2(ssn1.glmssn.SSS) # 0.2502698
GR2(ssn1.glmssn.ESG) # 0.2638648

CrossValidationStatsSSN(ssn1.glmssn.SSE.2)
#           bias     std.bias     RMSPE      RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001963672 -0.003251992 0.1915419 0.196302 0.9750998 0.8326948 0.9067688 0.9412516
 CrossValidationStatsSSN(ssn1.glmssn.SSE.3)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001606171 -0.002531696 0.1902606 0.1929672 0.9842688 0.8275862 0.9157088 0.9438059

#### RMSE for the selected model ####
ssn1.glmssn.SSE.preds <- as.data.frame(getPreds(ssn1.glmssn.SSE))
ssn1.glmssn.SSE.preds$cv.pred.untran <- 10^(ssn1.glmssn.SSE.preds$cv.pred)
ssn1.glmssn.SSE.preds <- merge(ssn1.glmssn.SSE.preds, obs.vars[,c('pid','STATION_KEY','FSS_26Aug14')], by = 'pid')

library(hydroGOF)
rmse(ssn1.glmssn.SSE.preds$cv.pred.untran, ssn1.glmssn.SSE.preds$FSS_26Aug14)
#[1] 11.6405