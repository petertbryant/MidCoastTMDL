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

save(ssn1.glmssn1.G, file = 'ssn1_glmssn1_G.Rdata')

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
# -262.4

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ LAT_RAW + sum_1095_days + 
#            XSLOPE_MAP + MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + 
#            PALITHERODRCA + PACLAYRCA + PASILTRCA + PASUSCEP5_DE + DAROADX + 
#            DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7589 -0.1181  0.0225  0.1695  0.8735 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)   -0.73488051  1.76588839   -0.42              0.67741    
# LAT_RAW        0.03237487  0.04002806    0.81              0.41888    
# sum_1095_days -0.00005075  0.00000619   -8.20 < 0.0000000000000002 ***
#   XSLOPE_MAP    -0.00091961  0.00260610   -0.35              0.72428    
# MIN_Z         -0.00002336  0.00003883   -0.60              0.54766    
# STRMPWR       -0.00003906  0.00003907   -1.00              0.31772    
# PDISRSA_1YR    0.00170276  0.00059831    2.85              0.00455 ** 
#   POWNRCA_PRI    0.00103754  0.00028567    3.63              0.00030 ***
#   PALITHERODRCA  0.00246861  0.00039475    6.25 < 0.0000000000000002 ***
#   PACLAYRCA      0.00462522  0.00277946    1.66              0.09651 .  
# PASILTRCA      0.00473061  0.00222767    2.12              0.03403 *  
#   PASUSCEP5_DE   0.00153080  0.00290062    0.53              0.59783    
# DAROADX        0.02211653  0.02759652    0.80              0.42313    
# DAPOPRCA2010   1.36400505  0.36908082    3.70              0.00023 ***
#   PAOWNRCA_AGR   0.00688937  0.00730871    0.94              0.34617    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00481569
# Exponential.tailup     range    609.28230894
# Exponential.taildown   parsill      0.01624457
# Exponential.taildown     range  12018.89967427
# Exponential.Euclid   parsill      0.01441616
# Exponential.Euclid     range 172748.01809782
# locID   parsill      0.00000167
# Nugget   parsill      0.02076763
# 
# Residual standard error: 0.2372
# Generalized R-squared: 0.2582

#Scaling the data left p-values the same. Residual standard error is halved. GR2 is the same. AIC way lower.
#Below, see that untransformed RMSE is improved to 3.6
#Effect = No diff for variable selection.
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
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.39063    0.06499   6.011  < 2e-16 ***
#   LAT_RAW        0.07282    0.09004   0.809  0.41891    
# sum_1095_days -0.30144    0.03674  -8.205  < 2e-16 ***
#   XSLOPE_MAP    -0.01542    0.04372  -0.353  0.72441    
# MIN_Z         -0.03382    0.05617  -0.602  0.54734    
# STRMPWR       -0.10247    0.10249  -1.000  0.31772    
# PDISRSA_1YR    0.08261    0.02903   2.846  0.00455 ** 
#   POWNRCA_PRI    0.05642    0.01554   3.632  0.00030 ***
#   PALITHERODRCA  0.13426    0.02147   6.254  < 2e-16 ***
#   PACLAYRCA      0.08025    0.04825   1.663  0.09664 .  
# PASILTRCA      0.09275    0.04367   2.124  0.03401 *  
#   PASUSCEP5_DE   0.02404    0.04556   0.528  0.59794    
# DAROADX        0.10404    0.12990   0.801  0.42342    
# DAPOPRCA2010   0.18214    0.04928   3.696  0.00023 ***
#   PAOWNRCA_AGR   0.07511    0.07969   0.943  0.34619    
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
#AIC - scaled = -1216.322

# > varcomp(ssn1.glmssn3.G)
# VarComp Proportion
# 1    Covariates (R-sq)   0.258208
# 2   Exponential.tailup   0.063511
# 3 Exponential.taildown   0.214240
# 4   Exponential.Euclid   0.190126
# 5                locID   0.000022
# 6               Nugget   0.273892

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
# -263.8

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + XSLOPE_MAP + 
#            MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
#            PACLAYRCA + PASILTRCA + PASUSCEP5_DE + DAROADX + DAPOPRCA2010 + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7283 -0.1172  0.0304  0.1706  0.8440 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.6904609  0.1408171    4.90 < 0.0000000000000002 ***
#   sum_1095_days -0.0000509  0.0000062   -8.22 < 0.0000000000000002 ***
#   XSLOPE_MAP    -0.0009469  0.0026068   -0.36              0.71653    
# MIN_Z         -0.0000223  0.0000389   -0.57              0.56659    
# STRMPWR       -0.0000394  0.0000391   -1.01              0.31313    
# PDISRSA_1YR    0.0016548  0.0005966    2.77              0.00568 ** 
#   POWNRCA_PRI    0.0010549  0.0002852    3.70              0.00023 ***
#   PALITHERODRCA  0.0024502  0.0003956    6.19 < 0.0000000000000002 ***
#   PACLAYRCA      0.0042438  0.0027439    1.55              0.12236    
# PASILTRCA      0.0051609  0.0021614    2.39              0.01719 *  
#   PASUSCEP5_DE   0.0013843  0.0029028    0.48              0.63358    
# DAROADX        0.0187943  0.0273267    0.69              0.49181    
# DAPOPRCA2010   1.3698467  0.3690368    3.71              0.00022 ***
#   PAOWNRCA_AGR   0.0069295  0.0073087    0.95              0.34337    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00489110
# Exponential.tailup     range    610.70878354
# Exponential.taildown   parsill      0.01603179
# Exponential.taildown     range  11917.49655454
# Exponential.Euclid   parsill      0.01487036
# Exponential.Euclid     range 169368.48524018
# locID   parsill      0.00000202
# Nugget   parsill      0.02077131
# 
# Residual standard error: 0.2378
# Generalized R-squared: 0.2559

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
#265.6

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + MIN_Z + 
#            STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + 
#            PASILTRCA + PASUSCEP5_DE + DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -0.727 -0.117  0.030  0.170  0.845 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.6911887  0.1406626    4.91 < 0.0000000000000002 ***
#   sum_1095_days -0.0000509  0.0000062   -8.21 < 0.0000000000000002 ***
#   MIN_Z         -0.0000256  0.0000378   -0.68              0.49850    
# STRMPWR       -0.0000430  0.0000378   -1.14              0.25544    
# PDISRSA_1YR    0.0016827  0.0005911    2.85              0.00453 ** 
#   POWNRCA_PRI    0.0010662  0.0002837    3.76              0.00018 ***
#   PALITHERODRCA  0.0024596  0.0003949    6.23 < 0.0000000000000002 ***
#   PACLAYRCA      0.0042358  0.0027440    1.54              0.12308    
# PASILTRCA      0.0050763  0.0021495    2.36              0.01844 *  
#   PASUSCEP5_DE   0.0013505  0.0029028    0.47              0.64188    
# DAROADX        0.0166264  0.0267176    0.62              0.53393    
# DAPOPRCA2010   1.3599476  0.3682000    3.69              0.00024 ***
#   PAOWNRCA_AGR   0.0071517  0.0072855    0.98              0.32659    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00497507
# Exponential.tailup     range    591.05508722
# Exponential.taildown   parsill      0.01592501
# Exponential.taildown     range  12033.06910221
# Exponential.Euclid   parsill      0.01479863
# Exponential.Euclid     range 166798.45810273
# locID   parsill      0.00000545
# Nugget   parsill      0.02077449
# 
# Residual standard error: 0.2377
# Generalized R-squared: 0.2557

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
ssn1.resid5.G <- residuals(ssn1.glmssn6.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn6.G, file = 'ssn1_glmssn6_G.Rdata')

AIC(ssn1.glmssn6.G)
#-267.4

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + MIN_Z + 
#            STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + 
#            PASILTRCA + DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7298 -0.1183  0.0266  0.1732  0.8454 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)    0.71678921  0.13036263    5.50 < 0.0000000000000002 ***
#   sum_1095_days -0.00005074  0.00000619   -8.19 < 0.0000000000000002 ***
#   MIN_Z         -0.00002855  0.00003716   -0.77              0.44262    
# STRMPWR       -0.00004148  0.00003765   -1.10              0.27089    
# PDISRSA_1YR    0.00167641  0.00059116    2.84              0.00469 ** 
#   POWNRCA_PRI    0.00105191  0.00028208    3.73              0.00021 ***
#   PALITHERODRCA  0.00242616  0.00038936    6.23 < 0.0000000000000002 ***
#   PACLAYRCA      0.00378665  0.00257603    1.47              0.14198    
# PASILTRCA      0.00492874  0.00213302    2.31              0.02111 *  
#   DAROADX        0.01626029  0.02672956    0.61              0.54315    
# DAPOPRCA2010   1.34415563  0.36697702    3.66              0.00027 ***
#   PAOWNRCA_AGR   0.00731847  0.00728007    1.01              0.31508    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill      0.0049590
# Exponential.tailup     range    593.3572497
# Exponential.taildown   parsill      0.0158219
# Exponential.taildown     range  12016.2922163
# Exponential.Euclid   parsill      0.0147994
# Exponential.Euclid     range 160188.7775105
# locID   parsill      0.0000151
# Nugget   parsill      0.0207850
# 
# Residual standard error: 0.2374
# Generalized R-squared: 0.2546

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
ssn1.resid5.G <- residuals(ssn1.glmssn7.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn7.G, file = 'ssn1_glmssn7_G.Rdata')

AIC(ssn1.glmssn7.G)
#-269.1
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + MIN_Z + 
#            STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + 
#            PASILTRCA + DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, 
#          family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                             "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#          EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7315 -0.1212  0.0259  0.1710  0.8430 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)    0.72178311  0.12911791    5.59 < 0.0000000000000002 ***
#   sum_1095_days -0.00005093  0.00000619   -8.23 < 0.0000000000000002 ***
#   MIN_Z         -0.00002843  0.00003715   -0.77              0.44426    
# STRMPWR       -0.00004191  0.00003766   -1.11              0.26622    
# PDISRSA_1YR    0.00166267  0.00059108    2.81              0.00503 ** 
#   POWNRCA_PRI    0.00106617  0.00028129    3.79              0.00016 ***
#   PALITHERODRCA  0.00240914  0.00038782    6.21 < 0.0000000000000002 ***
#   PACLAYRCA      0.00383287  0.00256938    1.49              0.13617    
# PASILTRCA      0.00489176  0.00212562    2.30              0.02164 *  
#   DAPOPRCA2010   1.36760461  0.36546110    3.74              0.00020 ***
#   PAOWNRCA_AGR   0.00724796  0.00728283    1.00              0.31994    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00497231
# Exponential.tailup     range    591.88431593
# Exponential.taildown   parsill      0.01581833
# Exponential.taildown     range  12035.30540401
# Exponential.Euclid   parsill      0.01418686
# Exponential.Euclid     range 150404.51234959
# locID   parsill      0.00000516
# Nugget   parsill      0.02079759
# 
# Residual standard error: 0.2362
# Generalized R-squared: 0.2548

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
ssn1.resid5.G <- residuals(ssn1.glmssn8.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn8.G, file = 'ssn1_glmssn8_G.Rdata')

AIC(ssn1.glmssn8.G)
#[1] -268.9

#Call:
# glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#          PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#          DAROADX + DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, 
#        family = "Gaussian", CorModels = c("locID", "Exponential.tailup", 
#                                           "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea", 
#        EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7639 -0.1144  0.0307  0.1764  0.8524 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.7091133  0.1318704    5.38 < 0.0000000000000002 ***
#   sum_1095_days -0.0000510  0.0000062   -8.22 < 0.0000000000000002 ***
#   STRMPWR       -0.0000425  0.0000376   -1.13              0.25907    
# PDISRSA_1YR    0.0017088  0.0005884    2.90              0.00379 ** 
#   POWNRCA_PRI    0.0010682  0.0002818    3.79              0.00016 ***
#   PALITHERODRCA  0.0024207  0.0003917    6.18 < 0.0000000000000002 ***
#   PACLAYRCA      0.0038113  0.0025880    1.47              0.14124    
# PASILTRCA      0.0049241  0.0021531    2.29              0.02247 *  
#   DAROADX        0.0162076  0.0268271    0.60              0.54592    
# DAPOPRCA2010   1.3448488  0.3670809    3.66              0.00027 ***
#   PAOWNRCA_AGR   0.0072776  0.0072764    1.00              0.31754    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00533803
# Exponential.tailup     range    583.05660475
# Exponential.taildown   parsill      0.01531400
# Exponential.taildown     range  12199.10723400
# Exponential.Euclid   parsill      0.01626782
# Exponential.Euclid     range 165296.28805979
# locID   parsill      0.00000342
# Nugget   parsill      0.02072151
# 
# Residual standard error: 0.2401
# Generalized R-squared: 0.2505

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
ssn1.resid5.G <- residuals(ssn1.glmssn9.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn9.G, file = 'ssn1_glmssn9_G.Rdata')

AIC(ssn1.glmssn9.G)
#-270.5

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAPOPRCA2010 + PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7655 -0.1170  0.0284  0.1737  0.8498 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)    0.7140431  0.1305337    5.47 < 0.0000000000000002 ***
#   sum_1095_days -0.0000512  0.0000062   -8.26 < 0.0000000000000002 ***
#   STRMPWR       -0.0000429  0.0000376   -1.14              0.25468    
# PDISRSA_1YR    0.0016952  0.0005883    2.88              0.00407 ** 
#   POWNRCA_PRI    0.0010824  0.0002810    3.85              0.00013 ***
#   PALITHERODRCA  0.0024039  0.0003902    6.16 < 0.0000000000000002 ***
#   PACLAYRCA      0.0038554  0.0025816    1.49              0.13572    
# PASILTRCA      0.0048912  0.0021463    2.28              0.02294 *  
#   DAPOPRCA2010   1.3681898  0.3655822    3.74              0.00020 ***
#   PAOWNRCA_AGR   0.0072050  0.0072790    0.99              0.32256    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00533967
# Exponential.tailup     range    582.05250647
# Exponential.taildown   parsill      0.01531053
# Exponential.taildown     range  12217.72625209
# Exponential.Euclid   parsill      0.01559629
# Exponential.Euclid     range 155652.63866512
# locID   parsill      0.00000549
# Nugget   parsill      0.02073447
# 
# Residual standard error: 0.2387
# Generalized R-squared: 0.2507

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
ssn1.resid5.G <- residuals(ssn1.glmssn10.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn10.G, file = 'ssn1_glmssn10_G.Rdata')

AIC(ssn1.glmssn10.G)
#-269.9

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + STRMPWR + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + 
#            DAROADX + DAPOPRCA2010, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "Exponential.taildown", 
#                        "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7631 -0.1130  0.0302  0.1776  0.8563 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)    0.68268077  0.12934197    5.28 < 0.0000000000000002 ***
#   sum_1095_days -0.00005097  0.00000619   -8.23 < 0.0000000000000002 ***
#   STRMPWR       -0.00004145  0.00003761   -1.10              0.27080    
# PDISRSA_1YR    0.00183011  0.00057799    3.17              0.00160 ** 
#   POWNRCA_PRI    0.00105172  0.00028132    3.74              0.00020 ***
#   PALITHERODRCA  0.00244455  0.00039017    6.27 < 0.0000000000000002 ***
#   PACLAYRCA      0.00422665  0.00254901    1.66              0.09769 .  
# PASILTRCA      0.00525671  0.00212527    2.47              0.01360 *  
#   DAROADX        0.01571540  0.02680360    0.59              0.55783    
# DAPOPRCA2010   1.39818367  0.36381756    3.84              0.00013 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00534365
# Exponential.tailup     range    589.68853789
# Exponential.taildown   parsill      0.01552918
# Exponential.taildown     range  12178.46031453
# Exponential.Euclid   parsill      0.01589428
# Exponential.Euclid     range 170082.81327477
# locID   parsill      0.00000352
# Nugget   parsill      0.02073308
# 
# Residual standard error: 0.2398
# Generalized R-squared: 0.2513

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
ssn1.resid5.G <- residuals(ssn1.glmssn11.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn11.G, file = 'ssn1_glmssn11_G.Rdata')

AIC(ssn1.glmssn11.G)
# -270.5

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
                           DAROADX + DAPOPRCA2010, 
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
ssn1.resid5.G <- residuals(ssn1.glmssn12.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn12.G, file = 'ssn1_glmssn12_G.Rdata')

AIC(ssn1.glmssn12.G)
#[1] -269.9

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
ssn1.resid5.G <- residuals(ssn1.glmssn13.G)
par(mfrow = c(1,2))
hist(ssn1.resid5.G)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid5.G)

resids2.df <- getSSNdata.frame(ssn1.resid5.G)
plot(resids2.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn13.G, file = 'ssn1_glmssn13_G.Rdata')

AIC(ssn1.glmssn13.G)
#-272.3

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + 
#            POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + DAPOPRCA2010, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7706 -0.1178  0.0275  0.1736  0.8538 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)    0.68380422  0.12867409    5.31 < 0.0000000000000002 ***
#   sum_1095_days -0.00005155  0.00000619   -8.33 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.00182787  0.00057794    3.16              0.00162 ** 
#   POWNRCA_PRI    0.00104037  0.00027956    3.72              0.00021 ***
#   PALITHERODRCA  0.00244835  0.00038855    6.30 < 0.0000000000000002 ***
#   PACLAYRCA      0.00438817  0.00254440    1.72              0.08499 .  
# PASILTRCA      0.00520216  0.00212507    2.45              0.01459 *  
#   DAPOPRCA2010   1.40521055  0.36205840    3.88              0.00011 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00560968
# Exponential.tailup     range    568.82582213
# Exponential.taildown   parsill      0.01514371
# Exponential.taildown     range  11948.86240292
# Exponential.Euclid   parsill      0.01587337
# Exponential.Euclid     range 162835.04808843
# locID   parsill      0.00000254
# Nugget   parsill      0.02075749
# 
# Residual standard error: 0.2396
# Generalized R-squared: 0.2495

#### ssn1.glmssn14.G ####
#Per Ver Hoef and Peterson 2010 using a two step procedure to first select model variables then select autocovariance functions and using
#ML instead of REML for AIC comparisons  (from Appendix B: "Note that REML was not used for parameter estimation
#in the first phase of model selection because a side effect of REML is that information criteria,
#such as AIC, can only be used if the explanatory variable set remains fixed (Verbeke and
#                                                                             Molenberghs 2000)")
#Now without STRMPWR 
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
#-271.4

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
#   Min      1Q  Median      3Q     Max 
# -0.7910 -0.1183  0.0209  0.1693  0.8261 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)    0.76764491  0.11600158    6.62 < 0.0000000000000002 ***
#   sum_1095_days -0.00005326  0.00000611   -8.72 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.00181398  0.00057977    3.13              0.00182 ** 
#   POWNRCA_PRI    0.00107803  0.00027948    3.86              0.00012 ***
#   PALITHERODRCA  0.00250477  0.00038740    6.47 < 0.0000000000000002 ***
#   PASILTRCA      0.00607549  0.00205231    2.96              0.00317 ** 
#   DAPOPRCA2010   1.44038103  0.36254696    3.97              0.00008 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00573265
# Exponential.tailup     range    596.62225459
# Exponential.taildown   parsill      0.01523420
# Exponential.taildown     range  11578.62652266
# Exponential.Euclid   parsill      0.01421794
# Exponential.Euclid     range 141888.83856282
# locID   parsill      0.00000641
# Nugget   parsill      0.02075213
# 
# Residual standard error: 0.2365
# Generalized R-squared: 0.249

#scaled
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
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.44297    0.03642  12.161  < 2e-16 ***
#   sum_1095_days -0.31639    0.03629  -8.719  < 2e-16 ***
#   PDISRSA_1YR    0.08801    0.02813   3.129  0.00182 ** 
#   POWNRCA_PRI    0.05862    0.01520   3.857  0.00012 ***
#   PALITHERODRCA  0.13621    0.02107   6.466  < 2e-16 ***
#   PASILTRCA      0.11911    0.04023   2.961  0.00316 ** 
#   DAPOPRCA2010   0.19232    0.04841   3.973    8e-05 ***
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
# > AIC(ssn1.glmssn14.G)
# [1] -1225.341

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

#############################
#### Ecological Approach ####
#The above uses a mathematical approach to model selection. The below uses a set of 
#variables that include potential interactions based on our knowledge of how landscape
#affect sediment dynamics. 
em <- read.csv('ecological_models.csv')
formulas <- as.character(em$Formula)
formulas <- formulas[formulas != ""]
formulas <- paste('log10_FSS_26Aug14 ~',formulas)
ssn.ecomod.list <- list()
for(i in 39:length(formulas)) {
  fit <- glmssn(as.formula(formulas[i]),
                            ssn1,
                            EstMeth = "ML",
                            CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                            addfunccol = "afvArea",
                            family = "Gaussian")
  save(fit, file = paste("model",i,".Rdata",sep=''))
  ssn.ecomod.list <- c(ssn.ecomod.list, list(fit))
}
ecomod.compare <- InfoCritCompare(ssn.ecomod.list)
write.csv(ecomod.compare, 'ecological_models_comparison.csv')
save(ssn.ecomod.list, file = 'allModelsList.Rdata')

#started at 12:58
print(Sys.time)

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
#We can use InfoCritCompare to evaluate model performance instead of calling
#CrossValidationStatsSSN for each individual model 
##ML##
variable.compare <- InfoCritCompare(list(ssn1.glmssn1.G,ssn1.glmssn3.G,ssn1.glmssn4.G,ssn1.glmssn5.G,
                                         ssn1.glmssn6.G,ssn1.glmssn7.G,ssn1.glmssn8.G,ssn1.glmssn9.G,
                                         ssn1.glmssn10.G,ssn1.glmssn11.G,ssn1.glmssn12.G,ssn1.glmssn13.G,
                                         ssn1.glmssn14.G))

##REML##
compare.models <- InfoCritCompare(list(ssn1.glmssn.ELG,ssn1.glmssn.ELE,ssn1.glmssn.EME,ssn1.glmssn.ESE,
                                       ssn1.glmssn.SSE,ssn1.glmssn.LLE,ssn1.glmssn.MME,ssn1.glmssn.MEE,
                                       ssn1.glmssn.SSS,ssn1.glmssn.ESG))
# 1  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 2  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 3  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 4  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 5  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 6  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 7  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 8  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 9  log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# 10 log10_FSS_26Aug14 ~ sum_1095_days + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + DAPOPRCA2010      REML
# Variance_Components neg2LogL    AIC       bias  std.bias  RMSPE    RAV std.MSPE
# 1     Exponential.tailup + LinearSill.taildown + Gaussian.Euclid + locID + Nugget   -223.5 -207.5 -0.0005968 -0.001020 0.1902 0.1920   0.9900
# 2  Exponential.tailup + LinearSill.taildown + Exponential.Euclid + locID + Nugget   -223.5 -207.5 -0.0007075 -0.001211 0.1897 0.1918   0.9885
# 3      Exponential.tailup + Mariah.taildown + Exponential.Euclid + locID + Nugget   -221.3 -205.3 -0.0009215 -0.001569 0.1896 0.1922   0.9864
# 4   Exponential.tailup + Spherical.taildown + Exponential.Euclid + locID + Nugget   -223.8 -207.8 -0.0008255 -0.001408 0.1894 0.1917   0.9878
# 5     Spherical.tailup + Spherical.taildown + Exponential.Euclid + locID + Nugget   -224.0 -208.0 -0.0008063 -0.001375 0.1894 0.1916   0.9878
# 6   LinearSill.tailup + LinearSill.taildown + Exponential.Euclid + locID + Nugget   -223.4 -207.4 -0.0006832 -0.001165 0.1895 0.1918   0.9852
# 7           Mariah.tailup + Mariah.taildown + Exponential.Euclid + locID + Nugget   -221.3 -205.3 -0.0009217 -0.001569 0.1896 0.1922   0.9864
# 8      Mariah.tailup + Exponential.taildown + Exponential.Euclid + locID + Nugget   -222.8 -206.8 -0.0008914 -0.001519 0.1895 0.1919   0.9870
# 9       Spherical.tailup + Spherical.taildown + Spherical.Euclid + locID + Nugget   -224.4 -208.4 -0.0007432 -0.001266 0.1893 0.1916   0.9877
# 10     Exponential.tailup + Spherical.taildown + Gaussian.Euclid + locID + Nugget   -223.9 -207.9 -0.0007204 -0.001227 0.1899 0.1919   0.9888
# cov.80 cov.90 cov.95
# 1  0.8301 0.9093 0.9438
# 2  0.8314 0.9106 0.9451
# 3  0.8276 0.9106 0.9413
# 4  0.8276 0.9093 0.9451
# 5  0.8289 0.9093 0.9438
# 6  0.8289 0.9093 0.9464
# 7  0.8276 0.9106 0.9413
# 8  0.8263 0.9106 0.9451
# 9  0.8289 0.9093 0.9438
# 10 0.8289 0.9132 0.9451

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

#### Untransformed RMSE for the selected model ####
ssn1.glmssn.SSE.preds <- as.data.frame(getPreds(ssn1.glmssn.SSE))
ssn1.glmssn.SSE.preds$cv.pred.untran <- 10^(ssn1.glmssn.SSE.preds$cv.pred)
ssn1.glmssn.SSE.preds <- merge(ssn1.glmssn.SSE.preds, obs.vars[,c('pid','STATION_KEY','FSS_26Aug14')], by = 'pid')

library(hydroGOF)
rmse(ssn1.glmssn.SSE.preds$cv.pred.untran, ssn1.glmssn.SSE.preds$FSS_26Aug14)
#[1] 11.6405


ssn1.glmssn3.G.preds <- as.data.frame(getPreds(ssn1.glmssn3.G))
ssn1.glmssn3.G.preds$cv.pred.untran <- 10^(ssn1.glmssn3.G.preds$cv.pred)
ssn1.glmssn3.G.preds <- merge(ssn1.glmssn3.G.preds, obs.vars[,c('pid','STATION_KEY','FSS_26Aug14')], by = 'pid')

library(hydroGOF)
ssn1.glmssn3.G.preds$cv.pred.untran.unscale <- ssn1.glmssn3.G.preds$cv.pred.untran*(max(ssn1.glmssn3.G.preds$FSS_26Aug14)-min(ssn1.glmssn3.G.preds$FSS_26Aug14)) + min(ssn1.glmssn3.G.preds$FSS_26Aug14)
ssn1.glmssn3.G.preds$FSS_26Aug14.unscale <- ssn1.glmssn3.G.preds$FSS_26Aug14*(max(ssn1.glmssn3.G.preds$FSS_26Aug14)-min(ssn1.glmssn3.G.preds$FSS_26Aug14)) + min(ssn1.glmssn3.G.preds$FSS_26Aug14)
rmse(ssn1.glmssn3.G.preds$cv.pred.untran.unscale, ssn1.glmssn3.G.preds$FSS_26Aug14.unscale)
#[1] 3.610497


