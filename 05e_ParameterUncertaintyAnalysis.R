library(SSN)
library(plyr)
library(RODBC)
library(ggplot2)
library(reshape2)
library(MASS)

options(stringsAsFactors = FALSE)

#Get fit object from 06_scenarios.R
fit <- models[[7]]
#Generate confidence intervals using assumptions of indepent normality
df_ci <- confint.glmssn(fit, level = .9)
df_ci <- cbind(df_ci, fit$estimates$betahat)
df_ci <- as.data.frame(df_ci)
df_ci$parms <- rownames(df_ci)
df_ci <- plyr::rename(df_ci, c("5 %" = "lci", "95 %" = "uci", "V3" = "est"))
df_ci_tmp <- df_ci[!df_ci$parms %in% c('(Intercept)', 'HDWTR1.21463119808789'),]
for (i in 1:nrow(df_ci)) {
  df_ci_tmp <- df_ci[i,]
  g <- ggplot(data = df_ci, aes(x = parms, y = est)) + #geom_point(aes(y = est)) + 
    geom_pointrange(aes(ymin = lci, ymax = uci)) + #facet_wrap( ~ parms) +
    geom_hline(yintercept = 0, colour = 'red', lwd = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
}

ssn_preds <- predict(fit, interval = "prediction", predpointsID = "preds")
#Brute force monte carlo simulations varying parameter estimates using
#multivariate normal distributions
obs <- getSSNdata.frame(fit)
set.seed(11)
randparm <- mvrnorm(n=500, mu = fit$estimates$betahat, Sigma = fit$estimates$covb)
randparm
bhats <- fit$estimates$betahat
dimnames(bhats)[[1]][7] <- 'HDWTR1'
dimnames(randparm)[[2]][7] <- 'HDWTR1'
for (i in 1:nrow(randparm)) {
  varied.preds <- predict.vary(betahat = as.data.frame(t(bhats)), ss = obs, 
                               r_vec = randparm[i,])
  varied.preds <- plyr::rename(varied.preds, c('BSTI_prd' = paste0("BSTI_prd_", i)))
  if (i == 1) {
    predtable <- varied.preds
  } else {
    predtable <- merge(predtable, varied.preds)
  }
}

#Check out the results
boxplot(10^t(predtable[1:10,])[2:501,])
