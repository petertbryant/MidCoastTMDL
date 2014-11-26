#### SubSSN Test ####
library(SSN)
library(stringr)

options('scipen' = 100)

ssn1 <- importSSN("C:/users/pbryant/desktop/midcoasttmdl-gis/revisedssn/LSN05/lsn.ssn", o.write=FALSE)
obs<- getSSNdata.frame(ssn1, Name = "Obs")
obs.complete <- read.csv("ssn_RF_data.csv")
obs.complete$SVN <- str_trim(obs.complete$SVN)
obs.fss2 <- read.csv('fss2_s2_data.csv')
obs.fss2 <- within(obs.fss2, rm(X))

vars <- c("STATION_KEY", "SITE_NAME", "SVN", "DATE","YEAR",names(obs.fss2))

obs.vars <- obs.complete[,vars]

obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", "afvArea",
                         "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", "HU_12_NAME", "HU_08", "LONG_RAW", "NHDHigh",
                         "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", "HU_12")],
                  obs.vars, 
                  by = "SVN",
                  all.x = TRUE)

obs.vars$log10_FSS_26Aug14 <- log10(obs.vars$FSS_26Aug14)
obs.vars$log10_sum_1095_days <- log10(obs.vars$sum_1095_day)
obs.vars$bin_PALITHERODRCA <- ifelse(obs.vars$PALITHERODRCA < 90,0,1)
obs.vars$sqrt_PADISRSA_1YR <- sqrt(obs.vars$PADISRSA_1YR)
obs.vars[which(obs.vars$XSLOPE_MAP <= 0),'XSLOPE_MAP'] <- 0.0001
obs.vars$log10_XSLOPE_MAP <- log10(obs.vars$XSLOPE_MAP)
obs.vars$log10_PASILTRCA <- log10(obs.vars$PASILTRCA)
obs.vars$log10_MIN_Z <- log10(obs.vars$MIN_Z)
obs.vars$log10_STRMPWR <- log10(obs.vars$STRMPWR + (1 - min(obs.vars$STRMPWR)))
obs.vars$sqrt_upDist <- sqrt(obs.vars$upDist)
obs.vars$log10_APOPRCA2010 <- log10(obs.vars$APOPRCA2010 + 1)
obs.vars$log10_PASUSCEP5_DE <- log10(obs.vars$PASUSCEP5_DE + 1)
obs.vars$log10_POWNRCA_FED <- log10(obs.vars$POWNRCA_FED + 1)
obs.vars$log10_POWNRCA_PRI <- log10(obs.vars$POWNRCA_PRI + 1)
obs.vars$bin_PAOWNRCA_AGR <- ifelse(obs.vars$PAOWNRCA_AGR > 0,1,0)

obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
row.names(obs.vars) <- obs.vars$pid
ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

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

#create data frame to store model scores
comp <- data.frame(model = character(), family = character(),GR2 = numeric(), AIC = numeric())
model.name <- paste('ssn1.glmssn',i,sep='')
comp <- rbind(comp, data.frame(model = model.name, family = model$args$family, GR2 = GR2(model), AIC = AIC(model)))

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

summary(ssn1.glmssn1.G)
varcomp(ssn1.glmssn1.G)
#Plot the residuals
ssn1.resid1.G <- residuals(ssn1.glmssn1.G)
par(mfrow = c(1,2))
hist(ssn1.resid1.G, breaks)
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
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -18.546  -4.888  -1.970   1.804  53.999 
# 
# Coefficients:
#   Estimate   Std. Error t value             Pr(>|t|)    
# (Intercept)    5.623025313  4.709658845   1.194              0.23287    
# sum_1095_days -0.001026898  0.000222037  -4.625 < 0.0000000000000002 ***
#   PALITHERODRCA  0.081223292  0.014791623   5.491 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.172814282  0.058168280   2.971              0.00306 ** 
#   XSLOPE_MAP    -0.074951697  0.087416124  -0.857              0.39148    
# PASILTRCA      0.215183845  0.079556273   2.705              0.00699 ** 
#   MIN_Z          0.000106806  0.001404543   0.076              0.93940    
# upDist        -0.000013656  0.000007113  -1.920              0.05523 .  
# APOPRCA2010    0.009817071  0.003285066   2.988              0.00289 ** 
#   PASUSCEP5_DE   0.012348931  0.100958514   0.122              0.90268    
# POWNRCA_FED   -0.031292494  0.016590162  -1.886              0.05964 .  
# POWNRCA_PRI    0.002480982  0.014925941   0.166              0.86803    
# PAOWNRCA_AGR   0.247389657  0.272615793   0.907              0.36444    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter      Estimate
# Exponential.tailup   parsill      0.000593
# Exponential.tailup     range 921433.225445
# Exponential.taildown   parsill     21.810591
# Exponential.taildown     range   7230.255432
# Exponential.Euclid   parsill     24.537307
# Exponential.Euclid     range  83565.776858
# locID   parsill      0.012102
# Nugget   parsill     25.746163
# 
# Residual standard error: 8.49157
# Generalized R-squared: 0.1706755
# > varcomp(ssn1.glmssn1.G)
# VarComp     Proportion
# 1    Covariates (R-sq) 0.170675500643
# 2   Exponential.tailup 0.000006818036
# 3 Exponential.taildown 0.250851075688
# 4   Exponential.Euclid 0.282211974044
# 5                locID 0.000139188399
# 6               Nugget 0.296115443190

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

#### ssn1.glmssn.S ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn.S <- glmssn(log10_FSS_26Aug14  ~ sum_1095_days + PALITHERODRCA + PADISRSA_1YR + 
                           XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
                           APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
                           PAOWNRCA_AGR, 
                         ssn1,
                         EstMeth = "REML",
                         CorModels = c("locID","Spherical.tailup", "Spherical.taildown"),
                         addfunccol = "afvArea",
                         family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn.S)
varcomp(ssn1.glmssn.S)
#Plot the residuals
ssn1.resid.S <- residuals(ssn1.glmssn.S)
par(mfrow = c(1,2))
hist(ssn1.resid.S)
hist(ssn1, "log10_FSS_26Aug14")

qqnorm(ssn1.resid.S)

resids.S.df <- getSSNdata.frame(ssn1.resid.S)
plot(resids.S.df[,"_fit_"],resids2.df[,"_resid_"])
#qqplot is linear and residual plot shows no trend. log transformation of response appears appropriate.

save(ssn1.glmssn.S, file = 'ssn1_glmssn_S.Rdata')

AIC(ssn1.glmssn.S)
#-80.82557

# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ sum_1095_days + PALITHERODRCA + 
#            PADISRSA_1YR + XSLOPE_MAP + PASILTRCA + MIN_Z + upDist + 
#            APOPRCA2010 + PASUSCEP5_DE + POWNRCA_FED + POWNRCA_PRI + 
#            PAOWNRCA_AGR, ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                                "Spherical.tailup", "Spherical.taildown"), addfunccol = "afvArea", 
#          EstMeth = "REML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.751179 -0.137770 -0.008497  0.146170  0.809293 
# 
# Coefficients:
#   Estimate    Std. Error t value             Pr(>|t|)    
# (Intercept)    1.0708502404  0.1047391674  10.224 < 0.0000000000000002 ***
#   sum_1095_days -0.0000634176  0.0000055522 -11.422 < 0.0000000000000002 ***
#   PALITHERODRCA  0.0030479876  0.0003528380   8.638 < 0.0000000000000002 ***
#   PADISRSA_1YR   0.0035425519  0.0016509183   2.146               0.0322 *  
#   XSLOPE_MAP    -0.0001570873  0.0025544645  -0.061               0.9510    
# PASILTRCA      0.0036637375  0.0015909580   2.303               0.0215 *  
#   MIN_Z         -0.0000814072  0.0000388250  -2.097               0.0363 *  
#   upDist        -0.0000001813  0.0000001368  -1.325               0.1857    
# APOPRCA2010    0.0002317467  0.0000953237   2.431               0.0153 *  
#   PASUSCEP5_DE  -0.0016464008  0.0026316778  -0.626               0.5318    
# POWNRCA_FED   -0.0015223025  0.0003874449  -3.929              0.00009 ***
#   POWNRCA_PRI    0.0000656226  0.0003835786   0.171               0.8642    
# PAOWNRCA_AGR   0.0066393822  0.0079166850   0.839               0.4019    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Spherical.tailup   parsill     0.000000502
# Spherical.tailup     range 16812.282854215
# Spherical.taildown   parsill     0.024413530
# Spherical.taildown     range 11157.613316094
# locID   parsill     0.004142475
# Nugget   parsill     0.021679210
# 
# Residual standard error: 0.2241333
# Generalized R-squared: 0.388666
# > varcomp(ssn1.glmssn.S)
# VarComp    Proportion
# 1  Covariates (R-sq) 0.38866603832
# 2   Spherical.tailup 0.00000611499
# 3 Spherical.taildown 0.29709578431
# 4              locID 0.05041105489
# 5             Nugget 0.26382100750

#### ssn1.glmssn1 ####
start.time <- Sys.time()
print(start.time)
ssn1.glmssn1 <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + bin_PALITHERODRCA + sqrt_PADISRSA_1YR + 
                         log10_XSLOPE_MAP + log10_PASILTRCA + log10_MIN_Z + sqrt_upDist + 
                         log10_APOPRCA2010 + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI + 
                         bin_PAOWNRCA_AGR, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                       addfunccol = "afvArea",
                       family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn1)

#### ssn1.glmssn1 ####
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ log10_sum_1095_days + bin_PALITHERODRCA + 
#            sqrt_PADISRSA_1YR + log10_XSLOPE_MAP + log10_PASILTRCA + 
#            log10_MIN_Z + sqrt_upDist + log10_APOPRCA2010 + log10_PASUSCEP5_DE + 
#            log10_POWNRCA_FED + log10_POWNRCA_PRI + bin_PAOWNRCA_AGR, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.70001 -0.11048  0.01134  0.16472  0.81896 
# 
# Coefficients:
#   Estimate  Std. Error t value             Pr(>|t|)    
# (Intercept)          3.61730762  0.57063708   6.339 < 0.0000000000000002 ***
#   log10_sum_1095_days -0.82411383  0.10142135  -8.126 < 0.0000000000000002 ***
#   bin_PALITHERODRCA    0.14191390  0.02720501   5.216 < 0.0000000000000002 ***
#   sqrt_PADISRSA_1YR    0.00108812  0.00849765   0.128              0.89814    
# log10_XSLOPE_MAP    -0.02072140  0.01113362  -1.861              0.06310 .  
# log10_PASILTRCA      0.37365774  0.22981999   1.626              0.10439    
# log10_MIN_Z         -0.08644950  0.02610882  -3.311              0.00097 ***
#   sqrt_upDist         -0.00002882  0.00013168  -0.219              0.82678    
# log10_APOPRCA2010    0.01418311  0.01679748   0.844              0.39873    
# log10_PASUSCEP5_DE  -0.10902844  0.03680660  -2.962              0.00315 ** 
#   log10_POWNRCA_FED   -0.04878066  0.01638483  -2.977              0.00300 ** 
#   log10_POWNRCA_PRI    0.03369402  0.01553072   2.170              0.03035 *  
#   bin_PAOWNRCA_AGR     0.03539024  0.03249651   1.089              0.27647    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter        Estimate
# Exponential.tailup   parsill      0.00387923
# Exponential.tailup     range    677.48029416
# Exponential.taildown   parsill      0.01827873
# Exponential.taildown     range  11987.05595284
# Exponential.Euclid   parsill      0.01278944
# Exponential.Euclid     range 186073.58948814
# locID   parsill      0.00000207
# Nugget   parsill      0.02065162
# 
# Residual standard error: 0.2357988
# Generalized R-squared: 0.2741876
# > varcomp(ssn1.glmssn1)
# VarComp    Proportion
# 1    Covariates (R-sq) 0.27418758523
# 2   Exponential.tailup 0.05063917907
# 3 Exponential.taildown 0.23860913515
# 4   Exponential.Euclid 0.16695237186
# 5                locID 0.00002700334
# 6               Nugget 0.26958472535
# 
# AIC(ssn1.glmssn1)
# [1] -223.7547
# 
##Plot the residuals
# ssn1.resid1 <- residuals(ssn1.glmssn1)
# par(mfrow = c(1,2))
# hist(ssn1.resid1, breaks = 100)
# hist(ssn1, "FSS_26Aug14")

#### ssn1.glmssn2 ####
#Using the variables indicated as significant in the above ssn1.glmssn1
#and then also changing the correlation models to spherical for tail up and tail down
start.time <- Sys.time()
print(start.time)
ssn1.glmssn2 <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + bin_PALITHERODRCA + log10_MIN_Z + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("locID","Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"),
                       addfunccol = "afvArea",
                       family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn2)
#### ssn1.glmssn2 summary ####
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ log10_sum_1095_days + bin_PALITHERODRCA + 
#            log10_MIN_Z + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID", 
#                                                                "Spherical.tailup", "Spherical.taildown", "Exponential.Euclid"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.720637 -0.128808 -0.004385  0.161257  0.863685 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)          4.41150    0.38525  11.451 < 0.0000000000000002 ***
#   log10_sum_1095_days -0.85014    0.10098  -8.419 < 0.0000000000000002 ***
#   bin_PALITHERODRCA    0.13415    0.02714   4.944 < 0.0000000000000002 ***
#   log10_MIN_Z         -0.11772    0.02252  -5.228 < 0.0000000000000002 ***
#   log10_PASUSCEP5_DE  -0.13761    0.03521  -3.908              0.00010 ***
#   log10_POWNRCA_FED   -0.04997    0.01646  -3.036              0.00248 ** 
#   log10_POWNRCA_PRI    0.04097    0.01538   2.663              0.00790 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter       Estimate
# Spherical.tailup   parsill     0.00586751
# Spherical.tailup     range   725.69654261
# Spherical.taildown   parsill     0.01512397
# Spherical.taildown     range 10087.32625469
# Exponential.Euclid   parsill     0.01201286
# Exponential.Euclid     range 91319.26275710
# locID   parsill     0.00000798
# Nugget   parsill     0.02067752
# 
# Residual standard error: 0.2317107
# Generalized R-squared: 0.2566409
# 
# VarComp   Proportion
# 1  Covariates (R-sq) 0.2566408575
# 2   Spherical.tailup 0.0812382754
# 3 Spherical.taildown 0.2093979328
# 4 Exponential.Euclid 0.1663233027
# 5              locID 0.0001104314
# 6             Nugget 0.2862892002

# "Exponential.tailup" (default), "LinearSill.tailup", "Spherical.tailup", "Mariah.tailup" "Epanech.tailup"; 
# "Exponential.taildown" (default), "LinearSill.taildown", "Spherical.taildown", "Mariah.taildown", "Epanech.taildown"; 
# "Spherical.Euclid", "Gaussian.Euclid", "Exponential.Euclid" (default), "Cauchy.Euclid"
# 
# AIC(ssn1.glmssn2)
# [1] -256.1144
# 
# ssn1.resid2 <- residuals(ssn1.glmssn2)
# par(mfrow = c(1,2))
# hist(ssn1.resid2, breaks = 100)
# hist(ssn1, "FSS_26Aug14")

#### ssn1.glmssn3 ####
#Use the significant variables but use the same autocovariance functions except for the non-important ones
start.time <- Sys.time()
print(start.time)
ssn1.glmssn3 <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + bin_PALITHERODRCA + log10_MIN_Z + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI, 
                       ssn1,
                       EstMeth = "REML",
                       CorModels = c("Exponential.taildown", "Exponential.Euclid"),
                       addfunccol = "afvArea",
                       family = "Gaussian")
end.time <- Sys.time()
print(end.time)
print(end.time - start.time)

summary(ssn1.glmssn3)
#### ssn1.glmssn3 summary ####
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ log10_sum_1095_days + bin_PALITHERODRCA + 
#            log10_MIN_Z + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("Exponential.taildown", 
#                                                                "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.717776 -0.127442 -0.003731  0.160735  0.863771 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)          4.37933    0.38421  11.398 < 0.0000000000000002 ***
#   log10_sum_1095_days -0.84022    0.10075  -8.339 < 0.0000000000000002 ***
#   bin_PALITHERODRCA    0.13612    0.02712   5.019 < 0.0000000000000002 ***
#   log10_MIN_Z         -0.12031    0.02238  -5.376 < 0.0000000000000002 ***
#   log10_PASUSCEP5_DE  -0.13763    0.03522  -3.907              0.00010 ***
#   log10_POWNRCA_FED   -0.05025    0.01637  -3.069              0.00222 ** 
#   log10_POWNRCA_PRI    0.04033    0.01527   2.642              0.00842 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter   Estimate
# Exponential.taildown   parsill     0.0209
# Exponential.taildown     range  9125.5378
# Exponential.Euclid   parsill     0.0117
# Exponential.Euclid     range 92582.3434
# Nugget   parsill     0.0210
# 
# Residual standard error: 0.2314826
# Generalized R-squared: 0.2578336
# varcomp(ssn1.glmssn3)
# VarComp Proportion
# 1    Covariates (R-sq)  0.2578336
# 2 Exponential.taildown  0.2893477
# 3   Exponential.Euclid  0.1621260
# 4               Nugget  0.2906927
# 
# AIC(ssn1.glmssn3)
# [1] -260.9818
# 
# ssn1.resid3 <- residuals(ssn1.glmssn3)
# par(mfrow = c(1,2))
# hist(ssn1.resid13, breaks = 100)
# hist(ssn1, "FSS_26Aug14")

#### ssn1.glmssn.NULL ####
# start.time <- Sys.time()
# print(start.time)
# ssn1.glmssn.NULL <- glmssn(log10_FSS_26Aug14  ~ log10_sum_1095_days + bin_PALITHERODRCA + log10_MIN_Z + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI, 
#                        ssn1,
#                        EstMeth = "REML",
#                        CorModels = c("locID"),
#                        addfunccol = "afvArea",
#                        family = "Gaussian")
# end.time <- Sys.time()
# print(end.time)
# print(end.time - start.time)
# 
# summary(ssn1.glmssn.NULL)
# Call:
#   glmssn(formula = log10_FSS_26Aug14 ~ log10_sum_1095_days + bin_PALITHERODRCA + 
#            log10_MIN_Z + log10_PASUSCEP5_DE + log10_POWNRCA_FED + log10_POWNRCA_PRI, 
#          ssn.object = ssn1, family = "Gaussian", CorModels = c("locID"), 
#          addfunccol = "afvArea", EstMeth = "REML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.767239 -0.143024 -0.005353  0.140465  0.873587 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)          4.93759    0.30525  16.175 < 0.0000000000000002 ***
#   log10_sum_1095_days -0.99022    0.07644 -12.954 < 0.0000000000000002 ***
#   bin_PALITHERODRCA    0.17795    0.02082   8.548 < 0.0000000000000002 ***
#   log10_MIN_Z         -0.13131    0.01792  -7.328 < 0.0000000000000002 ***
#   log10_PASUSCEP5_DE  -0.11363    0.02666  -4.262              0.00002 ***
#   log10_POWNRCA_FED   -0.05137    0.01140  -4.507              0.00001 ***
#   log10_POWNRCA_PRI    0.04442    0.01414   3.142              0.00174 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Covariance Parameters:
#   Covariance.Model Parameter Estimate
# locID   parsill   0.0266
# Nugget   parsill   0.0218
# 
# Residual standard error: 0.2199623
# Generalized R-squared: 0.469224
# > varcomp(ssn1.glmssn.NULL)
# VarComp Proportion
# 1 Covariates (R-sq)  0.4692240
# 2             locID  0.2921727
# 3            Nugget  0.2386033
# AIC(ssn1.glmssn.NULL)
# [1] -211.868
# 
# ssn1.resid.NULL <- residuals(ssn1.glmssn.NULL)
# par(mfrow = c(1,2))
# hist(ssn1.resid.NULL, breaks = 100)
# hist(ssn1, "FSS_26Aug14")

###################################################
### check the residuals
###################################################
ssn1.resid1 <- residuals(ssn1.glmssn1)
names( getSSNdata.frame(ssn1.resid1) )
plot(ssn1.resid1)

ssn1.resid1.P <- residuals(ssn1.glmssn1.P)
par(mfrow = c(1,2))
hist(ssn1.resid1.P, breaks = 100)
hist(ssn1, "FSS_26Aug14")
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

cv.out <- CrossValidationSSN(ssn1.glmssn1.P)
par(mfrow = c(1, 2))
plot(ssn1.glmssn1.P$sampinfo$z,
     cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(ssn1)[, "FSS_26Aug14"]),
      cv.out[, "cv.se"], pch = 19,
      xlab = "Observed Data", ylab = "LOOCV Prediction SE")
###################################################
### cross validation stats
###################################################
CrossValidationStatsSSN(ssn1.glmssn1)

###################################################
### R2
###################################################
GR2(ssn1.glmssn1.P)
varcomp(ssn1.glmssn1.P)
varcomp(ssn1.glmssn2)

