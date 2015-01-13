# -----------------------------------------------------------
# SSN
library(SSN)
library(stringr)
library(MASS)
library(plyr)
library(reshape2)
ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN05/lsn.ssn", predpts = "preds", o.write = TRUE)
#ssn1 <- importSSN('C:/users/pbryant/desktop/midcoasttmdl-gis/revisedssn/lsn05/lsn.ssn', o.write = TRUE)
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
#obs.vars <- obs.vars[!is.na(obs.vars$STRMPWR),]

#don't run when going through. only for data exploration
#obs.vars <- arrange(obs.vars, STATION_KEY, desc(YEAR))
#obs.vars.sub <- obs.vars[!duplicated(obs.vars$STATION_KEY),]

#### Response TRANSFORMATION####
names(obs.vars)
####FSS_26Aug14####
# hist(obs.vars.sub$FSS_26Aug14)
# plot(density(obs.vars.sub$FSS_26Aug14))
# hist(log10(obs.vars.sub$FSS_26Aug14))
# plot(density(log10(obs.vars.sub$FSS_26Aug14)))
# shapiro.test(log10(obs.vars.sub$FSS_26Aug14))
# ks.test(log10(obs.vars.sub$FSS_26Aug14), "pnorm")
#log should be just fine
obs.vars$log10_FSS_26Aug14 <- log10(obs.vars$FSS_26Aug14)

#### Explanatory Transformations ####
#sum_1095_days####
# qqnorm(log(obs.vars.sub$sum_1095_days), pch = 16)
# qqline(log(obs.vars.sub$sum_1095_days), col = 'green', lty = 2)
# hist(obs.vars.sub$sum_1095_days)
# hist(log(obs.vars.sub$sum_1095_days)) #looks good
# shapiro.test(obs.vars.sub$sum_1095_days)
# shapiro.test(log10(obs.vars.sub$sum_1095_days)) #p-value = 0.00274
# ks.test(log(obs.vars.sub$sum_1095_days), "pnorm") #p-value < 2.2e-16
# plot(density(((log10(obs.vars.sub$sum_1095_days))^(0.02) - 1)/0.02))
# #or we could try the boxcox
# boxcox(FSS_26Aug14 ~ sum_1095_days, data = obs.vars.sub, lambda = seq(-0.1, 0.1, .01))
# ks.test((obs.vars.sub$sum_1095_days^(0.02) - 1)/0.02, 'pnorm')
# plot(density(((obs.vars.sub$sum_1095_days)^(0.02) - 1)/0.02))
#log is about the same as boxcox. so we go with log
#obs.vars$log10_sum_1095_days <- log10(obs.vars$sum_1095_day)

#PALITHERODRCA####
# hist(obs.vars$PALITHERODRCA, 100)
# hist(1/(obs.vars$PALITHERODRCA))
# hist(asin(sqrt(obs.vars$PALITHERODRCA)))
# shapiro.test(log10(1/(obs.vars$PALITHERODRCA+1)))
# 
# boxcox(FSS_26Aug14 ~ PALITHERODRCA, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# hist((obs.vars.sub$PALITHERODRCA^(0.06) - 1)/0.06)
# ks.test((obs.vars.sub$PALITHERODRCA^(0.06) - 1)/0.06, 'pnorm')
# shapiro.test((obs.vars.sub$PALITHERODRCA^(0.06) - 1)/0.06)
# #This is very much either on or off. And looking at the partial dependence plot
# #there seems to be a break at 90% on the effect with FSS. We will use that as our
# #cut to convert this variable to binary
#obs.vars$bin_PALITHERODRCA <- ifelse(obs.vars$PALITHERODRCA < 90,0,1)

#PADISRSA_1YR####
# hist(obs.vars$PADISRSA_1YR)
# hist((sqrt(obs.vars$PADISRSA_1YR)))
# hist((log(obs.vars$PADISRSA_1YR+1)))
# x <- sqrt(obs.vars.sub$PADISRSA_1YR)
# # Test for normality
# qqnorm(x)
# hist(x, breaks = 10)
# plot(density(x))
# levene.test(x, obs.vars.sub$FSS_26Aug14)
# # Test Statistic = 1.4788, p-value = 0.02038 # p > 0.05 looks good.
# kurtosis(x, na.rm=TRUE)
# # [1] -0.3640666 
# skewness(x, na.rm=TRUE, type = 2)
# # [1] -0.8385455 for n = 100 needs to be between -0.391 and 0.391 see http://www.amstat.org/publications/jse/v19n2/doane.pdf
# shapiro.test(x)
# store <- data.frame('pvalue' = rep(NA,1000))
# for (i in 1:1000) {
#   y <- shapiro.test(sample(x,20)) #p-value = 9.691e-14
#   store$pvalue[i] = y$p.value
# }
# length(store[(store$pvalue < 0.05),])/1000*100
# #Let's go with sqrt() for this one
# boxcox(FSS_26Aug14 ~ PADISRSA_1YR, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# hist((obs.vars.sub$PADISRSA_1YR^(0.13) - 1)/0.13)
# plot(density(((obs.vars.sub$PADISRSA_1YR^(0.13) - 1)/0.13)))
# ks.test((obs.vars.sub$PADISRSA_1YR^(0.13) - 1)/0.13, 'pnorm')
# shapiro.test((obs.vars.sub$PADISRSA_1YR^(0.13) - 1)/0.13)
#stick with sqrt
#obs.vars$sqrt_PADISRSA_1YR <- sqrt(obs.vars$PADISRSA_1YR)

#XSLOPE_MAP####
# hist(obs.vars$XSLOPE_MAP)
# hist(log10(obs.vars$XSLOPE_MAP+3))
# plot(density(log10(obs.vars$XSLOPE_MAP+3)))
# shapiro.test(sqrt(obs.vars$XSLOPE_MAP))
# 
# boxcox(FSS_26Aug14 ~ XSLOPE_MAP, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# hist(((obs.vars.sub$XSLOPE_MAP+3)^(0.1) - 1)/0.1)
# plot(density(((obs.vars.sub$XSLOPE_MAP+3)^(0.1) - 1)/0.1))
# ks.test((obs.vars.sub$XSLOPE_MAP^(0.1) - 1)/0.1, 'pnorm')
# shapiro.test((obs.vars.sub$XSLOPE_MAP^(0.1) - 1)/0.1)
# #log is about the same as boxcox. so we go with log
# #but to take the log we need to adjust the negative slopes
# #given that negative and zero slopes are likely due to mapping errors
# #we will convert them to very small values just above 0.
#obs.vars[which(obs.vars$XSLOPE_MAP <= 0),'XSLOPE_MAP'] <- 0.0001
#obs.vars$log10_XSLOPE_MAP <- log10(obs.vars$XSLOPE_MAP)

#PASILTRCA####
# hist(obs.vars$PASILTRCA)
# hist(log10(obs.vars$PASILTRCA))
# plot(density(log10(obs.vars$PASILTRCA)))
# shapiro.test(log10(obs.vars$PASILTRCA))
# 
# boxcox(FSS_26Aug14 ~ PASILTRCA, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# hist((obs.vars.sub$PASILTRCA^(0.09) - 1)/0.09)
# plot(density((obs.vars.sub$PASILTRCA^(0.09) - 1)/0.09))
# ks.test((obs.vars.sub$PASILTRCA^(0.09) - 1)/0.09, 'pnorm')
# shapiro.test((obs.vars.sub$PASILTRCA^(0.09) - 1)/0.09)
#log is about the same as boxcox. so we go with log
#obs.vars$log10_PASILTRCA <- log10(obs.vars$PASILTRCA)

#MIN_Z####
# hist(obs.vars$MIN_Z, 10)
# hist(log10(obs.vars$MIN_Z), 10)
# plot(density(log10(obs.vars$MIN_Z)))
# shapiro.test(log10(obs.vars$MIN_Z))
# 
# boxcox(FSS_26Aug14 ~ MIN_Z, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# hist((obs.vars.sub$MIN_Z^(0.08) - 1)/0.08)
# plot(density((obs.vars.sub$MIN_Z^(0.08) - 1)/0.08))
# ks.test((obs.vars.sub$MIN_Z^(0.08) - 1)/0.08, 'pnorm')
# shapiro.test((obs.vars.sub$MIN_Z^(0.08) - 1)/0.08)
#log is about the same as boxcox. so we go with log
#obs.vars$log10_MIN_Z <- log10(obs.vars$MIN_Z)

#STRMPWR####
# hist(obs.vars$STRMPWR)
# plot(density(obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)]))
# hist(log10(obs.vars.sub$STRMPWR))
# boxcox(FSS_26Aug14 ~ STRMPWR, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# plot(density(log10(obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)] + (1 - obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)]))))
# shapiro.test(log10(obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)] + (1 - obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)]))
# ks.test(log10(obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)] + (1 - obs.vars.sub$STRMPWR[!is.na(obs.vars.sub$STRMPWR)])), "pnorm")
#obs.vars[which(obs.vars$STRMPWR <= 0),'STRMPWR'] <- 0.0001
#obs.vars$log10_STRMPWR <- log10(obs.vars$STRMPWR + (1 - min(obs.vars$STRMPWR)))

#"LAT_RAW"####
# hist(obs.vars.sub$LAT_RAW)
# plot(density(obs.vars.sub$LAT_RAW))
# shapiro.test(obs.vars.sub$LAT_RAW)
# hist((obs.vars$LAT_RAW))
# boxcox(FSS_26Aug14 ~ LAT_RAW, data = obs.vars.sub, lambda = seq(0, 0.3, .01))
# plot(density(log10(obs.vars.sub$LAT_RAW[!is.na(obs.vars.sub$LAT_RAW)] + (1 - obs.vars.sub$LAT_RAW[!is.na(obs.vars.sub$LAT_RAW)]))))
# shapiro.test(log10(obs.vars.sub$LAT_RAW[!is.na(obs.vars.sub$LAT_RAW)] + (1 - obs.vars.sub$LAT_RAW[!is.na(obs.vars.sub$LAT_RAW)]))
# ks.test(log10(obs.vars.sub$LAT_RAW[!is.na(obs.vars.sub$LAT_RAW)] + (1 - obs.vars.sub$LAT_RAW[!is.na(obs.vars.sub$LAT_RAW)])), "pnorm")
#Don't gain much with transformation. Leave untransformed.

#upDist####
# hist(obs.vars$upDist)
# plot(density(obs.vars.sub$upDist))
# shapiro.test(sqrt(obs.vars.sub$upDist))
# hist(sqrt(obs.vars.sub$upDist))
#obs.vars$sqrt_upDist <- sqrt(obs.vars$upDist)

#APOPRCA2010####
# hist(obs.vars$APOPRCA2010)
# plot(density(obs.vars.sub$APOPRCA2010))
# shapiro.test(log10(obs.vars.sub$APOPRCA2010 + 1))
# hist(log10(obs.vars.sub$APOPRCA2010 + 1))
#obs.vars$log10_APOPRCA2010 <- log10(obs.vars$APOPRCA2010 + 1)

#PASUSCEP5_DE####
# hist(obs.vars$PASUSCEP5_DE)
# plot(density(obs.vars.sub$PASUSCEP5_DE))
# shapiro.test(log10(obs.vars.sub$PASUSCEP5_DE + 1))
# hist(log10(obs.vars.sub$PASUSCEP5_DE + 1))
#obs.vars$log10_PASUSCEP5_DE <- log10(obs.vars$PASUSCEP5_DE + 1)

#POWNRCA_FED####
#hist(obs.vars$POWNRCA_FED)
#hist(((log10(obs.vars$POWNRCA_FED + 1))))
#Including transformation in case we want to use it not sure how much gain we get
#obs.vars$log10_POWNRCA_FED <- log10(obs.vars$POWNRCA_FED + 1)

#POWNRCA_PRI####
# hist(obs.vars$POWNRCA_PRI)
# hist(log10(obs.vars$POWNRCA_PRI + 1))
#Including transformation in case we want to use it not sure how much gain we get
#obs.vars$log10_POWNRCA_PRI <- log10(obs.vars$POWNRCA_PRI + 1)

#PAOWNRCA_AGR####
# hist(obs.vars$PAOWNRCA_AGR)
# hist(log10(obs.vars$PAOWNRCA_AGR + 1))
# summary(obs.vars$PAOWNRCA_AGR)
#obs.vars$bin_PAOWNRCA_AGR <- ifelse(obs.vars$PAOWNRCA_AGR > 0,1,0)

#### Variable Scaling #### 
melted <- melt(obs.vars[,c(names(obs.vars)[24:39])])
min.max <- ddply(melted, .(variable), summarize, min_val = min(value), max_val = max(value))
obs.vars[,c(names(obs.fss2),'log10_FSS_26Aug14')] <- as.data.frame(lapply(obs.vars[,c(names(obs.fss2),'log10_FSS_26Aug14')],
                                                   function(x) {(x-min(x))/(max(x)-min(x))}))

#### Put the data back into the ssn ####
#Now that we have the transformed variables we put them back in the 
#SSN object
obs.vars <- obs.vars[match(obs$pid, obs.vars$pid),]
row.names(obs.vars) <- obs.vars$pid
ssn1 <- putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

#save the ssn object to the github folder
#writeSSN(ssn1, filename = 'bugs.ssn')
