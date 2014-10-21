# -----------------------------------------------------------
# SSN
library(SSN)
library(stringr)
ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/lsn.ssn", o.write=FALSE)
obs<- getSSNdata.frame(ssn1, Name = "Obs")
obs.complete <- read.csv("ssn_RF_data.csv")
obs.complete$SVN <- str_trim(obs.complete$SVN)
obs.fss2 <- read.csv('fss2_s2_data.csv')
obs.fss2 <- within(obs.fss2, rm(X))

vars <- c("STATION_KEY", "SITE_NAME", "SVN", "YEAR",names(obs.fss2))

obs.vars <- obs.complete[,vars]

obs.vars <- merge(obs[,c("SVN","rid", "ratio", "locID", "netID", "pid", "upDist",  "afvArea",
                         "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", "HU_12_NAME", "HU_08", "LONG_RAW", "LAT_RAW", "NHDHigh",
                         "NHDh_Reach", "NHDP21_Rea", "NHDP12_COM", "HU_10", "HU_12")],
                  obs.vars, 
                  by = "SVN",
                  all.x = TRUE)
obs.vars <- obs.vars[obs$pid,]

row.names(obs.vars) <- obs.vars$pid

putSSNdata.frame(obs.vars, ssn1, Name = 'Obs')

obs.vars <- arrange(obs.vars, STATION_KEY, desc(YEAR))
obs.vars.sub <- obs.vars[!duplicated(obs.vars$STATION_KEY),]

model1 <- c("FSS_26Aug14","sum_1095_days","PADISRSA_1YR","PALITHERODRCA","XSLOPE_MAP","PASILTRCA","MIN_Z","POWNRCA_PRI")
model2 <- c("FSS_26Aug14","sum_1095_days","PADISRSA_1YR","PALITHERODRCA","XSLOPE_MAP","PASILTRCA","MIN_Z","PAOWNRSA_FED")
model3 <- c("FSS_26Aug14","sum_1095_days","PADISRSA_1YR","PALITHERODRCA","XSLOPE_MAP","PASILTRCA","MIN_Z","PAOWNRCA_FED")
model4 <- c("FSS_26Aug14","sum_1095_days","PADISRSA_1YR","PALITHERODRCA","XSLOPE_MAP","PASILTRCA","MIN_Z","PAOWNRCA_AGR")
model5 <- c("FSS_26Aug14","sum_1095_days","PADISRSA_1YR","PALITHERODRCA","XSLOPE_MAP","PASILTRCA","MIN_Z","PAOWNRSA_AGR")

qqnorm(log(obs.vars.sub$sum_1095_days), pch = 16)
qqline(log(obs.vars.sub$sum_1095_days), col = 'green', lty = 2)
hist(obs.vars.sub$sum_1095_days)
hist(log(obs.vars.sub$sum_1095_days)) #looks good
shapiro.test(obs.vars.sub$sum_1095_days)
shapiro.test(log10(obs.vars.sub$sum_1095_days)) #p-value = 0.00274
ks.test(log(obs.vars.sub$sum_1095_days), "pnorm") #p-value < 2.2e-16
#let's go with the log for this one

hist(obs.vars$PADISRSA_1YR)
hist((sqrt(obs.vars$PADISRSA_1YR)))
hist((log(obs.vars$PADISRSA_1YR+1)))

x <- sqrt(obs.vars.sub$PADISRSA_1YR)
# Test for normality
qqnorm(x)
hist(x, breaks = 10)
plot(density(x))
levene.test(x, obs.vars.sub$FSS_26Aug14)
# Test Statistic = 1.4788, p-value = 0.02038 # p > 0.05 looks good.
kurtosis(x, na.rm=TRUE)
# [1] -0.3640666 
skewness(x, na.rm=TRUE, type = 2)
# [1] -0.8385455 for n = 100 needs to be between -0.391 and 0.391 see http://www.amstat.org/publications/jse/v19n2/doane.pdf
shapiro.test(x)
store <- data.frame('pvalue' = rep(NA,1000))
for (i in 1:1000) {
  y <- shapiro.test(sample(x,20)) #p-value = 9.691e-14
  store$pvalue[i] = y$p.value
}
length(store[(store$pvalue < 0.05),])/1000*100
#Let's go with sqrt() for this one

hist(obs.vars.sub$PALITHERODRCA, 100)
hist(1/(obs.vars$PALITHERODRCA))
hist(asin(sqrt(obs.vars$PALITHERODRCA)))
shapiro.test(log10(1/(obs.vars$PALITHERODRCA+1)))
obs.vars.sub$BPALITHERODRCA <- ifelse(obs.vars.sub$PALITHERODRCA < )

hist(obs.vars$XSLOPE_MAP)
hist(log10(obs.vars$XSLOPE_MAP+1))
shapiro.test(sqrt(obs.vars$XSLOPE_MAP))

hist(obs.vars$PASILTRCA)
hist(log10(obs.vars$PASILTRCA))
shapiro.test(log10(obs.vars$PASILTRCA))

hist(obs.vars$MIN_Z, 10)
hist(log10(obs.vars$MIN_Z), 10)
shapiro.test(sqrt(obs.vars$MIN_Z))

hist(obs.vars$POWNRCA_PRI, 100)
shapiro.test(obs.vars$POWNRCA_PRI)

hist(obs.vars$PAOWNRSA_FED)
hist(obs.vars$PAOWNRCA_FED)
hist(obs.vars$PAOWNRCA_AGR)
hist(obs.vars$PAOWNRSA_AGR)
