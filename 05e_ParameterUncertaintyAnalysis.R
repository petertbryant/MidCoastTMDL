library(SSN)
library(plyr)
library(RODBC)
library(ggplot2)
library(reshape2)

options(stringsAsFactors = FALSE)

# #
# nbatches <- 1
# nreps <- 10000
# for (k in 1:nbatches) { 
#   pb <- txtProgressBar(min = 0, max = nreps, style = 3)
#   predtable <- matrix(rep(0,nreps*nrow(dhists)),ncol=nreps) # create empty matrix for preds
#   for (j in 1:nreps) {
#    #i <- sample(1:11,1,prob=weights) # select a model
#     fitmod <- ssn1_glmssn1 # get that model
#     randparm <- mvrnorm(n=1,mu=fitmod$estimates$betahat,Sigma=fitmod$estimates$covb) # get parameter ests
#     mm <- model.matrix(terms(fitmod),dhists)
#     randeff <- ranef(fitmod)$huc12[match(as.character(dhists$huc12),row.names(ranef(fitmod)$huc12[1])),1] # use known random effect value where available
#     # For other hucs, generate random value and apply it to each record
#     randhuc12 <- rnorm(length(unique(dhists$huc12)),0,as.numeric(VarCorr(fitmod)[1])^.5)
#     names(randhuc12) <- as.character(unique(dhists$huc12))
#     randeff[is.na(randeff)] <- randhuc12[match(as.character(dhists$huc12),names(randhuc12))]
#     predtable[,j] <- round(plogis(0,location=((mm%*%randparm)+randeff),scale=1,lower.tail=F),3)	
#     setTxtProgressBar(pb,j)
#   }
# }


#Get the obs data frame
# obs_all <- getSSNdata.frame(ssn1, Name = "Obs")

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

tmp <- glmssn(log10_BSTI ~ sum_1095_days + XSLOPE_MAP + MIN_Z + KFACT_MARCA + POP_DARCA + OWN_FED_PRCA + TYPEF_PRCA + ROADLEN_DRSA + DIS_3YR_PRSA + HDWTR,
              EstMeth = "REML",
              ssn1,
              CorModels = c("locID",'Exponential.Euclid','Exponential.taildown'),
              addfunccol = "afvArea",
              family = "Gaussian")



df_ci <- confint.glmssn(tmp, level = .9)
df_ci <- cbind(df_ci, tmp$estimates$betahat)
df_ci <- as.data.frame(df_ci)
df_ci$parms <- rownames(df_ci)
df_ci <- plyr::rename(df_ci, c("5 %" = "lci", "95 %" = "uci", "V3" = "est"))
df_ci_tmp <- df_ci[!df_ci$parms %in% c('(Intercept)', 'HDWTR100'),]
for (i in 1:nrow(df_ci)) {
  df_ci_tmp <- df_ci[i,]
  g <- ggplot(data = df_ci_tmp, aes(x = parms, y = est)) + #geom_point(aes(y = est)) + 
    geom_pointrange(aes(ymin = lci, ymax = uci)) + #facet_wrap( ~ parms) +
    geom_hline(yintercept = 0, colour = 'red', lwd = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
}

df_prd <- getSSNdata.frame(ssn1, Name = 'obs_prd')

# #Covariance matrix of fixed effects
# covb <- tmp$estimates$covb
# #corb <- cov2cor(covb)
# #Cholesky decomposition
# L <- chol(covb)
# #Number of variables
# nvars <- dim(L)[1]
# #Number of parameter sets to generate
# nparsets <- 500
# #Make matrix of paramter sets with uniform distribtion
# r <- t(L) %*% matrix(runif(nvars*nparsets), nrow = nvars, ncol = nparsets)
# r <- t(r)
# #Pull in the components used to calculate betahat
# X2 <- tmp$sampinfo$X
# Vi <- tmp$estimates$Vi
# zt <- tmp$sampinfo$z
# betahat <- tmp$estimates$betahat
# 
# bi <- t(betahat)
# covbi <- t(X2) %*% Vi %*% zt %*% bi
# solve(interim)
# 
# df_pred <- data.frame(zobs = df_prd$log10_BSTI, 
#                       b = runif(n = length(df_prd$log10_BSTI), 
#                                 min = df_ci[1, 1], 
#                                 max = df_ci[1, 2]),
#                       beta1 = runif(n = length(df_prd$log10_BSTI), 
#                                     min = df_ci[2, 1], 
#                                     max = df_ci[2, 2]),
#                       beta2 = runif(n = length(df_prd$log10_BSTI), 
#                                     min = df_ci[3, 1], 
#                                     max = df_ci[3, 2]),
#                       x = df_prd$x,
#                       y = df_prd$y)
# 
# InfoCritCompare2(list(tmp))
# predict.glmssn(tmp, predpointsID = )
# 
# vars <- c('sum_1095_days', 'XSLOPE_MAP', 'MIN_Z', 'KFACT_MARCA', 'POP_DARCA', 
#           'OWN_FED_PRCA', 'TYPEF_PRCA', 'ROADLEN_DRSA', 'DIS_3YR_PRSA')
# cvars <- c()
# for (i in 1:length(vars)) {
#   tmp.cvars <- pcor[,vars[i]][apply(pcor, MARGIN = 1, FUN = function(x) x > 0.4)[,vars[i]]]
#   cvars <- c(cvars, tmp.cvars)
# }