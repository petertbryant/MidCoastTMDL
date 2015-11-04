library(SSN)
library(stringr)
library(MASS)
library(plyr)
library(reshape2)

options('scipen' = 2)

load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_vi_median_20150722_1035.RData")
load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_20150722_1035.RData")
#### Variable selection ####
# Values drop off and then level out. Arbitrarily going with 50% of the variables.
# grab all variable names with median values > 1.004880e-04 = 50% of the data
# This 50% of the data reflects 50% of the original list of variables prior to scaling
# Scaling had the effect of dropping variables that were all 0s anyway.
fss2.s2.col <- fss2.s1.vi.median[1:ceiling(nrow(fss2.s1.vi.median) / 2), ][, 1]
#fss2.s2.col <- c("FSS_26Aug14",(fss2.s1.vi.median[,'var_name']))
fss2.s2 <- fss2.s1[, colnames(fss2.s1) %in% fss2.s2.col]

obs.complete <- read.csv("ssn_RF_data.csv")
obs.complete$SVN <- str_trim(obs.complete$SVN)
ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN05/lsn.ssn", 
                  predpts = "preds_up", o.write = TRUE)

obs<- getSSNdata.frame(ssn1, Name = "Obs")
vars <- c("STATION_KEY", "SVN", "DATE","YEAR",'FSS_26Aug14')

load(file = "matrix_svd.Rdata")

obs.vars <- obs.complete[,vars]

obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", 
                         "afvArea", "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", 
                         "HU_12_NAME", "HU_08", "LONG_RAW", "NHDHigh",
                         "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", 
                         "HU_12")],
                  obs.vars, 
                  by = "SVN",
                  all.x = TRUE)

#Normalize by maximum range
melted <- melt(fss2.s2[,names(fss2.s2[,-grep('_P',names(fss2.s2))])])
min.max <- ddply(melted, .(variable), 
                 summarize, 
                 min_val = min(value), 
                 max_val = max(value))
fss2.s2[,names(fss2.s2[,-grep('_P',names(fss2.s2))])] <- as.data.frame(lapply(
  fss2.s2[,names(fss2.s2[,-grep('_P',names(fss2.s2))])], function(x) {((x-min(x))/
      (max(x)-min(x)))*100}))

fss2.s2.cov <- cov(fss2.s2)
fss2.s2.cov.svd <- svd(fss2.s2.cov)
# save(fss2.s2.cov.svd, file = 'matrix_svd.Rdata')
# save(fss2.s2.cov, file = 'matrix_cov.Rdata')

d <- fss2.s2.cov.svd$d

# d/d[1]
# 
# plot(d,log="y")

v <- fss2.s2.cov.svd$v

# v[50,]
junk <- data.frame(s1 = (d[1]/d[1] * v[1,])^2, 
                   s2 = (d[2]/d[1] * v[2,])^2,
#                    s3 = (d[3]/d[1] * v[3,])^2,
#                    s4 = (d[4]/d[1] * v[4,])^2,
                   var = names(fss2.s2))
arrange(junk, desc(s1))[,c(1,3)]
arrange(junk, desc(s2))[,c(2,3)]
arrange(junk, desc(s3))[,c(3,5)]
# d[2] * v[2,]


S4 <- diag(d[1:12])

v4 <- t(v[1:12,])

sp4 <- v4 %*% S4

fss2.sp4 <- as.matrix(fss2.s2) %*% sp4

obs.vars <- cbind(obs.vars, fss2.sp4)
#obs.vars <- cbind(obs.vars, fss2.s2)

obs.vars$log10_FSS_26Aug14 <- log10(obs.vars$FSS_26Aug14)

# melted <- melt(obs.vars[,c('log10_FSS_26Aug14',as.character(1:4))])
# min.max <- ddply(melted, .(variable), 
#                  summarize, 
#                  min_val = min(value), 
#                  max_val = max(value))
# obs.vars[,c(as.character(1:2),'log10_FSS_26Aug14')] <- as.data.frame(lapply(
#   obs.vars[,c(as.character(1:2),'log10_FSS_26Aug14')], function(x) {(x-min(x))/
#       (max(x)-min(x))}))

obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
row.names(obs.vars) <- obs.vars$pid

names(obs.vars)[names(obs.vars) %in% 1:12] <- paste("V", names(obs.vars)
                                                   [names(obs.vars) %in% 1:12], 
                                                   sep = "")

ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

fit4 <- glmssn(as.formula(obs.vars[c('log10_FSS_26Aug14', paste("V", 1:12, 
                                                                sep = ""))]),
               ssn.object = ssn1,
               EstMeth = "ML",
               CorModels = c('Exponential.Euclid','Exponential.taildown'))
summary(fit4)
varcomp(fit4)
InfoCritCompare(list(fit4))

fit1 <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14',
                                      "V1", "V2")]),
           ssn.object = ssn1,
           EstMeth = "ML",
           CorModels = c('Exponential.Euclid','Exponential.taildown'))
summary(fit1)
varcomp(fit1)
#Plot the residuals
fit1_resids <- residuals(fit1)
par(mfrow = c(1,2))
hist(fit1_resids)
hist(ssn1, "FSS_26Aug14")

qqnorm(fit1_resids)

resids.df <- getSSNdata.frame(fit1_resids)
plot(resids.df[,"_fit_"],resids.df[,"_resid_"])

fit2 <- glmssn(as.formula(obs.vars[,c('log10_FSS_26Aug14',
                                      "V1")]),
               ssn.object = ssn1,
               EstMeth = "ML",
               CorModels = c('Exponential.Euclid','Exponential.taildown'))
summary(fit2)

fit3 <- glmssn(log10_FSS_26Aug14 ~ 1,
               ssn.object = ssn1,
               EstMeth = "ML",
               CorModels = c('Exponential.Euclid','Exponential.taildown'))
summary(fit3)

InfoCritCompare(list(fit1, fit2, fit3))
