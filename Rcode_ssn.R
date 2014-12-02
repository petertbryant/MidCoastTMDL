library(SSN)
library(plyr)

options(scipen = 100)

# ssn1 <- importSSN('bugs.ssn')
# obs <- getSSNdata.frame(ssn1, Name = "Obs")
# obs <- rename(obs, c(
#               "STATION_KE" = "STATION_KEY",
#               "APOPRCA201" = "APOPRCA2010",
#               "sum_60_day" = "sum_60_days",
#               "sum_180_da" = "sum_180_days",
#               "sum_365_da" = "sum_365_days",
#               "sum_1095_d" = "sum_1095_days",
#               "PPT_1981_2" = "PPT_1981_2010",
#               "FSS_26Aug1" = "FSS_26Aug14",
#               "PLITHERODR" = "PLITHERODRCA",
#               "PLITHERODR.1" = "PLITHERODRSA",
#               "PSILTRCA" = "PSILTRCA",
#               "PDISRCA_1Y" = "PDISRCA_1YR",
#               "PDISRCA_3Y" = "PDISRCA_3YR",
#               "PDISRSA_1Y" = "PDISRSA_1YR",
#               "POWNRCA_PR" = "POWNRCA_PRI",
#               "POWNRCA_FE" = "POWNRCA_FED",
#               "POWNRSA_FE" = "POWNRSA_FED",
#               "PALITHEROD" = "PALITHERODRCA",
#               "PALITHEROD.1" = "PALITHERODRSA",
#               "PASILT_CLA" = "PASILT_CLAYRCA",
#               "PADISRSA_1" = "PADISRSA_1YR",
#               "PAOWNRCA_F" = "PAOWNRCA_FED",
#               "PAOWNRCA_A" = "PAOWNRCA_AGR",
#               "PAOWNRSA_F" = "PAOWNRSA_FED",
#               "PAOWNRSA_A" = "PAOWNRSA_AGR",
#               "log10_FSS_" = "log10_FSS_26Aug14",
#               "log10_sum_" = "log10_sum_1095_days",
#               "sqrt_PADIS" = "sqrt_PADISRSA_1YR",
#               "bin_PALITH" = "bin_PALITHERODRCA",
#               "log10_XSLO" = "log10_XSLOPE_MAP",
#               "log10_PASI" = "log10_PASILTRCA",
#               "log10_MIN_" = "log10_MIN_Z"))
# ssn1 <- putSSNdata.frame(obs, ssn1, Name = 'Obs')

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
#Because variance is greater than mean, we are transforming the response variable.
ssn1.glmssn1 <- glmssn(log10_FSS_26Aug14  ~ LAT_RAW + upDist + APOPRCA2010 + sum_1095_days + XSLOPE_MAP + 
                         MIN_Z + STRMPWR + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
                         PACLAYRCA + PASILTRCA + PALFCVRSA + PLFUNRSA, 
                       ssn1,
                       EstMeth = "ML",
                        CorModels = c("locID","Exponential.tailup", "LinearSill.taildown",
                                      "Gaussian.Euclid"), addfunccol = "afvArea",family = "Gaussian")
summary(ssn1.glmssn1)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ LAT_RAW + upDist + APOPRCA2010 + 
#            sum_1095_days + XSLOPE_MAP + MIN_Z + STRMPWR + PDISRSA_1YR + 
#            POWNRCA_PRI + PALITHERODRCA + PACLAYRCA + PASILTRCA + PALFCVRSA + 
#            PLFUNRSA, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                            "Exponential.tailup", "LinearSill.taildown", "Gaussian.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.750119 -0.123449  0.008115  0.155868  0.820372 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.2086582196  1.4101493672   0.148              0.88241    
# LAT_RAW        0.0141312440  0.0319749647   0.442              0.65865    
# upDist        -0.0000005617  0.0000003219  -1.745              0.08142 .  
# APOPRCA2010    0.0001836552  0.0000895467   2.051              0.04061 *  
#   sum_1095_days -0.0000535969  0.0000060515  -8.857 < 0.0000000000000002 ***
#   XSLOPE_MAP     0.0053603935  0.0027792543   1.929              0.05413 .  
# MIN_Z         -0.0000426462  0.0000379442  -1.124              0.26140    
# STRMPWR       -0.0000514359  0.0000390860  -1.316              0.18858    
# PDISRSA_1YR    0.0013372260  0.0006176840   2.165              0.03070 *  
#   POWNRCA_PRI    0.0008701609  0.0002779782   3.130              0.00181 ** 
#   PALITHERODRCA  0.0022739105  0.0003701083   6.144 < 0.0000000000000002 ***
#   PACLAYRCA      0.0036029579  0.0025135479   1.433              0.15215    
# PASILTRCA      0.0035600369  0.0020918941   1.702              0.08919 .  
# PALFCVRSA     -0.0057580383  0.0029607458  -1.945              0.05216 .  
# PLFUNRSA       0.0030494114  0.0010166632   2.999              0.00279 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter     Estimate
# Exponential.tailup   parsill     0.005991
# Exponential.tailup     range   559.702403
# LinearSill.taildown   parsill     0.016398
# LinearSill.taildown     range  5367.612905
# Gaussian.Euclid   parsill     0.008972
# Gaussian.Euclid     range 88276.459124
# locID   parsill     0.000127
# Nugget   parsill     0.020713
# 
# Residual standard error: 0.2284765
# Generalized R-squared: 0.2891661

varcomp(ssn1.glmssn1)
# VarComp  Proportion
# 1   Covariates (R-sq) 0.289166069
# 2  Exponential.tailup 0.081576192
# 3 LinearSill.taildown 0.223293811
# 4     Gaussian.Euclid 0.122177005
# 5               locID 0.001731488
# 6              Nugget 0.282055434

AIC(ssn1.glmssn1)
# [1] -272.0046

ssn1.glmssn2 <- glmssn(log10_FSS_26Aug14  ~ upDist + APOPRCA2010 + sum_1095_days + XSLOPE_MAP + 
                         PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
                         PASILTRCA + PALFCVRSA + PLFUNRSA, 
                       ssn1,
                       EstMeth = "ML",
                       CorModels = c("locID","Exponential.tailup", "LinearSill.taildown",
                                     "Gaussian.Euclid"), addfunccol = "afvArea",family = "Gaussian")
summary(ssn1.glmssn2)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ upDist + APOPRCA2010 + sum_1095_days + 
#            XSLOPE_MAP + PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
#            PASILTRCA + PALFCVRSA + PLFUNRSA, ssn.object = ssn1, family = "Gaussian", 
#          CorModels = c("locID", "Exponential.tailup", "LinearSill.taildown", 
#                        "Gaussian.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.80646 -0.11997  0.01288  0.15677  0.80692 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.8831581133  0.1501943105   5.880 < 0.0000000000000002 ***
#   upDist        -0.0000005661  0.0000003260  -1.736              0.08288 .  
# APOPRCA2010    0.0001830359  0.0000884738   2.069              0.03890 *  
#   sum_1095_days -0.0000560325  0.0000059974  -9.343 < 0.0000000000000002 ***
#   XSLOPE_MAP     0.0038026053  0.0026111623   1.456              0.14572    
# PDISRSA_1YR    0.0013130808  0.0006151520   2.135              0.03311 *  
#   POWNRCA_PRI    0.0009001382  0.0002764444   3.256              0.00118 ** 
#   PALITHERODRCA  0.0023322042  0.0003745540   6.227 < 0.0000000000000002 ***
#   PASILTRCA      0.0045988366  0.0019759172   2.327              0.02020 *  
#   PALFCVRSA     -0.0056106975  0.0029700479  -1.889              0.05925 .  
# PLFUNRSA       0.0029888468  0.0010157479   2.943              0.00335 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill     0.00823648
# Exponential.tailup     range   632.12549633
# LinearSill.taildown   parsill     0.01449263
# LinearSill.taildown     range  6371.08418432
# Gaussian.Euclid   parsill     0.01025917
# Gaussian.Euclid     range 86261.53029903
# locID   parsill     0.00000452
# Nugget   parsill     0.02065366
# 
# Residual standard error: 0.2316171
# Generalized R-squared: 0.2776071

varcomp(ssn1.glmssn2)
# VarComp    Proportion
# 1   Covariates (R-sq) 0.27760713069
# 2  Exponential.tailup 0.11091082107
# 3 LinearSill.taildown 0.19515498858
# 4     Gaussian.Euclid 0.13814796966
# 5               locID 0.00006086399
# 6              Nugget 0.27811822601

AIC(ssn1.glmssn2)
# [1] -274.6108

ssn1.glmssn3 <- glmssn(log10_FSS_26Aug14 ~ upDist + APOPRCA2010 + sum_1095_days + 
                       PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + PALFCVRSA + 
                       PLFUNRSA, 
                       ssn1,
                       EstMeth = "ML",
                       CorModels = c("locID","Exponential.tailup", "LinearSill.taildown",
                                     "Gaussian.Euclid"), addfunccol = "afvArea",family = "Gaussian")
summary(ssn1.glmssn3)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ upDist + APOPRCA2010 + sum_1095_days + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + PALFCVRSA + 
#            PLFUNRSA, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                            "Exponential.tailup", "LinearSill.taildown", "Gaussian.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.79710 -0.12511  0.01423  0.15217  0.81504 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.8863318634  0.1498957873   5.913 < 0.0000000000000002 ***
#   upDist        -0.0000005672  0.0000003245  -1.748              0.08087 .  
# APOPRCA2010    0.0001772103  0.0000885365   2.002              0.04568 *  
#   sum_1095_days -0.0000560459  0.0000059955  -9.348 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.0012931126  0.0006167507   2.097              0.03635 *  
#   POWNRCA_PRI    0.0008897663  0.0002765297   3.218              0.00135 ** 
#   PALITHERODRCA  0.0023314779  0.0003738505   6.236 < 0.0000000000000002 ***
#   PASILTRCA      0.0049399989  0.0019569127   2.524              0.01179 *  
#   PALFCVRSA     -0.0048844500  0.0029157446  -1.675              0.09430 .  
# PLFUNRSA       0.0026183547  0.0009792552   2.674              0.00766 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Exponential.tailup   parsill     0.00653471
# Exponential.tailup     range   640.13805800
# LinearSill.taildown   parsill     0.01623594
# LinearSill.taildown     range  5260.78979390
# Gaussian.Euclid   parsill     0.01014906
# Gaussian.Euclid     range 84940.65583375
# locID   parsill     0.00000484
# Nugget   parsill     0.02072066
# 
# Residual standard error: 0.2316143
# Generalized R-squared: 0.2756966

varcomp(ssn1.glmssn3)
# VarComp    Proportion
# 1   Covariates (R-sq) 0.27569658729
# 2  Exponential.tailup 0.08822989565
# 3 LinearSill.taildown 0.21921337793
# 4     Gaussian.Euclid 0.13702989361
# 5               locID 0.00006531561
# 6              Nugget 0.27976492991

AIC(ssn1.glmssn3)
# [1] -274.7873

start.time <- Sys.time()
print(start.time)
ssn1.glmssn4 <- glmssn(log10_FSS_26Aug14  ~ APOPRCA2010 + sum_1095_days +  
                         PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + 
                         PASILTRCA + PLFUNRSA, 
                       ssn1,
                       EstMeth = "ML",
                       CorModels = c("locID","Exponential.tailup", "LinearSill.taildown",
                                     "Gaussian.Euclid"), addfunccol = "afvArea",family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn4)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ APOPRCA2010 + sum_1095_days + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + PLFUNRSA, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Exponential.tailup", "LinearSill.taildown", "Gaussian.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.81446 -0.12924  0.01079  0.15462  0.80814 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    0.686449317  0.116551979   5.890 < 0.0000000000000002 ***
#   APOPRCA2010    0.000193020  0.000088019   2.193              0.02861 *  
#   sum_1095_days -0.000055667  0.000006039  -9.218 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.001246433  0.000614505   2.028              0.04287 *  
#   POWNRCA_PRI    0.000929898  0.000276511   3.363              0.00081 ***
#   PALITHERODRCA  0.002403055  0.000366187   6.562 < 0.0000000000000002 ***
#   PASILTRCA      0.005425862  0.001997000   2.717              0.00673 ** 
#   PLFUNRSA       0.003032563  0.000918040   3.303              0.00100 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter    Estimate
# Exponential.tailup   parsill     0.00125
# Exponential.tailup     range   762.46313
# LinearSill.taildown   parsill     0.01720
# LinearSill.taildown     range  5436.43099
# Gaussian.Euclid   parsill     0.01430
# Gaussian.Euclid     range 87765.54993
# locID   parsill     0.00408
# Nugget   parsill     0.02095
# 
# Residual standard error: 0.2403815
# Generalized R-squared: 0.2585819

varcomp(ssn1.glmssn4)
# VarComp Proportion
# 1   Covariates (R-sq) 0.25858188
# 2  Exponential.tailup 0.01603197
# 3 LinearSill.taildown 0.22073464
# 4     Gaussian.Euclid 0.18342053
# 5               locID 0.05237384
# 6              Nugget 0.26885714

AIC(ssn1.glmssn4)
# [1] -272.9861

start.time <- Sys.time()
print(start.time)
ssn1.glmssn5 <- glmssn(log10_FSS_26Aug14 ~ upDist + APOPRCA2010 + sum_1095_days + 
                         PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + PALFCVRSA + 
                         PLFUNRSA, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Spherical.tailup", "Spherical.taildown",
                                     "Exponential.Euclid"), addfunccol = "afvArea",family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn5)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ upDist + APOPRCA2010 + sum_1095_days + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + PALFCVRSA + 
#            PLFUNRSA, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                            "Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.79199 -0.11870  0.01587  0.15676  0.81371 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    0.8745573233  0.1625993898   5.379 < 0.0000000000000002 ***
#   upDist        -0.0000006155  0.0000003522  -1.748              0.08091 .  
# APOPRCA2010    0.0001808462  0.0000906998   1.994              0.04651 *  
#   sum_1095_days -0.0000552654  0.0000061580  -8.975 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.0011784078  0.0006268304   1.880              0.06049 .  
# POWNRCA_PRI    0.0009775245  0.0002823703   3.462              0.00057 ***
#   PALITHERODRCA  0.0022268745  0.0004074429   5.465 < 0.0000000000000002 ***
#   PASILTRCA      0.0051253229  0.0021046553   2.435              0.01511 *  
#   PALFCVRSA     -0.0044959819  0.0029612671  -1.518              0.12936    
# PLFUNRSA       0.0024949482  0.0009855494   2.532              0.01155 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Spherical.tailup   parsill      0.0080113
# Spherical.tailup     range    774.8934761
# Spherical.taildown   parsill      0.0131319
# Spherical.taildown     range   9789.4793013
# Exponential.Euclid   parsill      0.0168734
# Exponential.Euclid     range 174916.3278718
# locID   parsill      0.0000131
# Nugget   parsill      0.0206540
# 
# Residual standard error: 0.2422473
# Generalized R-squared: 0.253741

CrossValidationStatsSSN(ssn1.glmssn5)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80   cov.90    cov.95
# 1 -0.001481071 -0.002504653 0.1886093 0.1917013 0.9833145 0.8288633 0.908046 0.9374202


start.time <- Sys.time()
print(start.time)
ssn1.glmssn6 <- glmssn(log10_FSS_26Aug14 ~ APOPRCA2010 + sum_1095_days + 
                         PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA +  
                         PLFUNRSA, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Spherical.tailup", "Spherical.taildown",
                                     "Exponential.Euclid"), addfunccol = "afvArea",family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn6)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ APOPRCA2010 + sum_1095_days + 
#            PDISRSA_1YR + POWNRCA_PRI + PALITHERODRCA + PASILTRCA + PLFUNRSA, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.80693 -0.12223  0.01584  0.15554  0.80235 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    0.676250301  0.134395682   5.032 < 0.0000000000000002 ***
#   APOPRCA2010    0.000195354  0.000090349   2.162              0.03091 *  
#   sum_1095_days -0.000055050  0.000006175  -8.915 < 0.0000000000000002 ***
#   PDISRSA_1YR    0.001152338  0.000625132   1.843              0.06566 .  
# POWNRCA_PRI    0.001031541  0.000281731   3.661              0.00027 ***
#   PALITHERODRCA  0.002290563  0.000397955   5.756 < 0.0000000000000002 ***
#   PASILTRCA      0.005726514  0.002124669   2.695              0.00719 ** 
#   PLFUNRSA       0.002907993  0.000927324   3.136              0.00178 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter         Estimate
# Spherical.tailup   parsill      0.007718290
# Spherical.tailup     range    735.585624394
# Spherical.taildown   parsill      0.013033033
# Spherical.taildown     range  10010.936409111
# Exponential.Euclid   parsill      0.022251313
# Exponential.Euclid     range 214720.881564404
# locID   parsill      0.000000763
# Nugget   parsill      0.020841600
# 
# Residual standard error: 0.2526757
# Generalized R-squared: 0.2421446

CrossValidationStatsSSN(ssn1.glmssn6)
#           bias     std.bias     RMSPE       RAV  std.MSPE    cov.80    cov.90    cov.95
# 1 -0.001634635 -0.002763872 0.1890935 0.1916397 0.9862844 0.8263091 0.9042146 0.9399745

#These models aren't an improvement over the model selection process that used only the susceptibility metric

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
