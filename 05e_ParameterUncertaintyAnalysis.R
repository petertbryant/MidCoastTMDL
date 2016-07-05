library(SSN)
library(plyr)
library(RODBC)
library(ggplot2)
library(reshape2)

options(stringsAsFactors = FALSE)

#Get the obs data frame
obs_all <- getSSNdata.frame(ssn1, Name = "Obs")

# #Split the observations up to generate an estimate and prediction set
# lng_rng <- seq_along(obs_all$SVN)
# lng_est <- sample(lng_rng, 500)
# lng_prd <- sample(lng_rng[-1 * lng_est], 260)
# obs_est <- obs_all[lng_est, ]
# obs_prd <- obs_all[lng_prd, ]
# 
# #the obs data frame needs to have NA's populated for the removed data
# obs_all[lng_prd,'log10_BSTI'] <- NA
# 
# #Put back the est only
# ssn1 <- putSSNdata.frame(obs_all, ssn1, Name = 'Obs')
# ssn1 <- putSSNdata.frame(obs_prd, ssn1, Name = 'Preds_obs')

tmp <- glmssn(as.formula(obs_all[,c('log10_BSTI',names(bsti.s2))]),
              EstMeth = "ML",
              ssn1,
              CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
              addfunccol = "afvArea",
              family = "Gaussian")

df_ci <- confint.glmssn(tmp)
df_prd <- getSSNdata.frame(ssn1, Name = 'obs_prd')

#Covariance matrix of fixed effects
covb <- tmp$estimates$covb
#corb <- cov2cor(covb)
#Cholesky decomposition
L <- chol(covb)
#Number of variables
nvars <- dim(L)[1]
#Number of parameter sets to generate
nparsets <- 500
#Make matrix of paramter sets with uniform distribtion
r <- t(L) %*% matrix(runif(nvars*nparsets), nrow = nvars, ncol = nparsets)
r <- t(r)
#Pull in the components used to calculate betahat
X2 <- tmp$sampinfo$X
Vi <- tmp$estimates$Vi
zt <- tmp$sampinfo$z
betahat <- tmp$estimates$betahat

bi <- t(betahat)
covbi <- t(X2) %*% Vi %*% zt %*% bi
solve(interim)

df_pred <- data.frame(zobs = df_prd$log10_BSTI, 
                      b = runif(n = length(df_prd$log10_BSTI), 
                                min = df_ci[1, 1], 
                                max = df_ci[1, 2]),
                      beta1 = runif(n = length(df_prd$log10_BSTI), 
                                    min = df_ci[2, 1], 
                                    max = df_ci[2, 2]),
                      beta2 = runif(n = length(df_prd$log10_BSTI), 
                                    min = df_ci[3, 1], 
                                    max = df_ci[3, 2]),
                      x = df_prd$x,
                      y = df_prd$y)

InfoCritCompare2(list(tmp))
predict.glmssn(tmp, predpointsID = )