library(randomForest)
library(reshape)
library(plyr)

source('funCorrelationPlots.R')

options(stringsAsFactors = FALSE)

load("C:/users/pbryant/desktop/midcoasttmdl/
     fss2_s1_vi_median_20150722_1035.RData")
load("C:/users/pbryant/desktop/midcoasttmdl/fss2_s1_20150722_1035.RData")

#### Variable selection ####
# Values drop off and then level out. Arbitrarily going with 50% of the variables.
# grab all variable names with median values > 1.004880e-04 = 50% of the data
# This 50% of the data reflects 50% of the original list of variables prior to scaling
# Scaling had the effect of dropping variables that were all 0s anyway.
fss2.s2.col <- fss2.s1.vi.median[1:ceiling(nrow(fss2.s1.vi.median) / 2), ][, 1]
#fss2.s2.col <- c("FSS_26Aug14",(fss2.s1.vi.median[,'var_name']))
fss2.s2 <- fss2.s1[, colnames(fss2.s1) %in% fss2.s2.col]

fss2.s2.col <- vars[vars$var %in% names(fss2.s2),]
#fss2.s2.col <- fss2.s2.col[fss2.s2.col$var != 'FSS_26Aug14',]
fss2.s2.col <- merge(fss2.s2.col, fss2.s1.vi.median[, c('var_name','median')], 
                     by.x = 'var', by.y = 'var_name', all.x = TRUE)
fss2.s2.col <- arrange(fss2.s2.col, desc(median))

#By category
all_keep <- c()
for (j in 1:length(unique(fss2.s2.col$Category))) {
  pcor <- cor(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 
                                      unique(fss2.s2.col$Category)[j], 'var']])
  pnames <- attr(pcor, "dimnames")[[1]]
  pkeep <- pnames
  for (i in length(pnames):1) {
    if (any(round(abs(pcor[i, ][-i]), 2) >= 0.2)) {
      if (i != 1) {
        pkeep <- pkeep[-i]
        pcor <- pcor[-i, -i, drop=FALSE]
      }
    } 
  }
  all_keep <- c(all_keep, pkeep)
}

# # #All together
# pcor <- cor(fss2.s1[,setdiff(fss2.s1.vi.median$var_name,"DATE")])
# pnames <- attr(pcor, "dimnames")[[1]]
# pkeep <- pnames
# for (i in length(pnames):1) {
#   if (any(round(abs(pcor[i,][-i]),2) >= 0.3)) {
#     if (i != 1) {
#       pkeep <- pkeep[-i]
#       pcor <- pcor[-i,-i,drop=FALSE]
#     }
#   } 
# }

#### Correlation plots ####
#names(fss2.s2.col) <- fss2.s2.col

#Method: Look for r<=0.2. If lower importance variable has r>=0.2 with any higher importance variables select the highest 
#importance variable unless that variable is already conceptually represented.

# Precip  
#png('precip_cor.png')
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Precipitation', 'var']],
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
#dev.off()
# keep "sum_1095_days"

#Disturbance
pairs(fss2.s2[ ,fss2.s2.col[fss2.s2.col$Category == 'Disturbance', 'var']],
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
# everything is coorelated
# keep "DIS_1YR_PARSA"

# Lithology/soils
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Lithology and soils', 
                            'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
# almost everything is coorelated
# Keep "EROD_PARCA", "SILT_CLAY_PARCA"

# Ownership 
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Land use', 'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
# Keep"ROADLEN_DRSA","OWN_FED_PRCA","POP_DARCA","OWN_PRI_PRCA","OWN_AGR_PARCA"

#Susceptibility
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Landslide susceptibility', 
                            'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist)
#All correlated. 
#Keep "SUSCEP5_PARCA"

#Location
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Location', 'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist)
#Keep "MIN_Z"

#Stream attributes
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Stream attributes', 
                            'var']],
      lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist)
#Keep "STRMPWR","XSLOPE_MAP"

#Further remove variables to reduce the influence of correlation on raising variable importance
#fss2.s2 <- fss2.s2[,colnames(fss2.s2) %in% keeps.s2]
fss2.s2 <- fss2.s2[, colnames(fss2.s2) %in% c('FSS_26Aug14', all_keep)]
colnames(fss2.s2)

# remove any NAs
fss2.s2 <- data.frame(na.omit(fss2.s2))
#write.csv(fss2.s2, 'fss2_s2_data_testing.csv')
#write.csv(fss2.s2, 'fss2_s2_data_scaled.csv')
