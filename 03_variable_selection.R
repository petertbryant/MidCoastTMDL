library(randomForest)
library(reshape)
library(plyr)

options(stringsAsFactors = FALSE)

#vars <- read.csv("VarNames_RF_v2.csv")
# load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_vi_median_STRMPWR_20151105_1429.RData")
# load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_STRMPWR_20151105_1429.RData")
load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_vi_median_20151216_0913.RData")
load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_20151216_0913.RData")

#### Variable selection ####
# Values drop off and then level out. Arbitrarily going with 50% of the variables.
# grab all variable names with median values > 1.004880e-04 = 50% of the data
# This 50% of the data reflects 50% of the original list of variables prior to scaling
# Scaling had the effect of dropping variables that were all 0s anyway.
fss2.s2.col <- fss2.s1.vi.median[1:ceiling(nrow(fss2.s1.vi.median) / 2), ][, 1]
fss2.s1.vi.median <- fss2.s1.vi.median[1:ceiling(nrow(fss2.s1.vi.median) / 2), ]
#fss2.s2.col <- c("FSS_26Aug14",(fss2.s1.vi.median[,'var_name']))
fss2.s2 <- fss2.s1[, colnames(fss2.s1) %in% c(fss2.s2.col,'HDWTR')]

# fss2.s2.col <- vars[vars$var %in% names(fss2.s2),]
# #fss2.s2.col <- fss2.s2.col[fss2.s2.col$var != 'FSS_26Aug14',]
# fss2.s2.col <- merge(fss2.s2.col, fss2.s1.vi.median[, c('var_name','median')], 
#                      by.x = 'var', by.y = 'var_name', all.x = TRUE)
# fss2.s2.col <- arrange(fss2.s2.col, desc(median))

#By category
correlation_threshold <- 0.4
all_keep <- c()
# for (j in 1:length(unique(fss2.s2.col$Category))) {
#   pcor <- cor(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 
#                                       unique(fss2.s2.col$Category)[j], 'var']])
#   pnames <- attr(pcor, "dimnames")[[1]]
#   pkeep <- pnames
#   for (i in length(pnames):1) {
#     if (any(round(abs(pcor[i, ][-i]), 2) >= correlation_threshold)) {
#       if (i != 1) {
#         pkeep <- pkeep[-i]
#         pcor <- pcor[-i, -i, drop=FALSE]
#       }
#     } 
#   }
#   all_keep <- c(all_keep, pkeep)
# }

# # #All together
pcor <- cor(fss2.s1[,setdiff(fss2.s1.vi.median$var_name,"DATE")])
pnames <- attr(pcor, "dimnames")[[1]]
pkeep <- pnames
for (i in length(pnames):1) {
  if (any(round(abs(pcor[i,][-i]),2) >= correlation_threshold)) {
    if (i != 1) {
      pkeep <- pkeep[-i]
      pcor <- pcor[-i,-i,drop=FALSE]
    }
  } 
}
all_keep <- pkeep

#Further remove variables to reduce the influence of correlation on raising variable importance
#fss2.s2 <- fss2.s2[,colnames(fss2.s2) %in% keeps.s2]
fss2.s2 <- fss2.s2[, c(all_keep,'HDWTR')]
colnames(fss2.s2)

#Test will all human influence variables - NOT YET RUN 10/30/2015 - THESE ARE ALL IN THERE ALREADY EXCEPT FOR OWN_FED WHICH IS HIGHLY CORRELATED WITH OWN_PRI
# fss2.s2 <- fss2.s1[, colnames(fss2.s1) %in% c('FSS_26Aug14', 'STRMPWR', 
#                                               'PPT_1981_2010', 'OWN_FED_PRCA',
#                                               'POP_DARCA', 'OWN_PRI_PRCA',
#                                               'OWN_AGR_PARCA', 'ROADLEN_DRSA',
#                                               'OWN_URB_PARCA')]

# remove any NAs
fss2.s2 <- data.frame(na.omit(fss2.s2))
# save(fss2.s2, file = 'fss2_s2_scaled.Rdata')
#write.csv(fss2.s2, 'fss2_s2_data_testing.csv')
#write.csv(fss2.s2, 'fss2_s2_data_scaled.csv')
