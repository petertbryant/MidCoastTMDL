#This file needs to be run in 64 bit R

library(randomForest)
library(reshape)
library(plyr)

source('funCorrelationPlots.R')

options(stringsAsFactors = FALSE)

vars <- read.csv("VarNames_RF_v2.csv")
bugs <- read.csv("ssn_RF_data.csv")

# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data ####
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss2.s1 <- vars[vars$keep == 1, ]

fss2.s1 <- bugs[, colnames(bugs) %in% vars.fss2.s1$var]

fss2.s1$DATE <- as.POSIXct(fss2.s1$DATE)

# #The below steps were necessary during initial processing and are not necessary now that we are 
# #using the modified dataset
# #remove NAs in response variable
# fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]
# # remove any NAs - this has been taken care of prior to this but for good measure we'll leave it in
# #This removes the whole row where there is an NA. There are columns that are all NA due to scaling 
# fss2.s1 <- data.frame(na.omit(fss2.s1))

# #Normalize the data to the same scale - Completed on 12/17 at 1618/1620. Did not affect rank order
# Will not normalize to same scale at this point in the data processing since it affects the 
# value of the %inc MSE because the scale of the response is so small. 
# fss2.s1$LONG_RAW <- abs(fss2.s1$LONG_RAW)
# fss2.s1 <- as.data.frame(lapply(fss2.s1,function(x) {(x-min(x))/(max(x)-min(x))}))
# 
#This removes those variables where all the values are 0
fss2.s1 <- fss2.s1[, setdiff(names(fss2.s1), c("X2year_count_60_days",
                                               "X10year_count_60_days", 
                                               "X25year_count_60_days", 
                                               "X50year_count_60_days",
                                               "X100year_count_60_days", 
                                               "X10year_count_180_days", 
                                               "X25year_count_180_days",
                                               "X50year_count_180_days", 
                                               "X100year_count_180_days"))]
#Need to convert Inf values to 0
fss2.s1[, grep("X", names(fss2.s1))] <- as.data.frame(sapply(
  fss2.s1[, grep("X",names(fss2.s1))], function(x) {
    replace(x, is.infinite(x),0)
    }
  ))

# mtry value
mtry.fss2.s1 <- as.integer(((ncol(fss2.s1) - 1) / 3), 0)

# initialize the variable importance df to store the importance scores from each run
fss2.s1.vi <- data.frame(matrix(nrow = ncol(fss2.s1) - 1, ncol = 50))
fss2.s1.visd <- data.frame(matrix(nrow = ncol(fss2.s1) - 1, ncol = 50))

fss2.s1.col <- colnames(fss2.s1)
fss2.s1.col <- fss2.s1.col[!(fss2.s1.col == "FSS_26Aug14")]

#Run the randomForest
# WARNING - Takes about 30 min
beg <- Sys.time()
print(beg)
set.seed(100)
for (i in 1:50) {
  fss2.s1.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s1,
                             mtry = mtry.fss2.s1,
                             ntree = 2000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s1.vi[, i] <- fss2.s1.rf$importance[, 1]
  fss2.s1.visd[, i] <- fss2.s1.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s1.vi[, 51]<- fss2.s1.col
fss2.s1.vi[, 52]<-c(1:length(fss2.s1.col))
colnames(fss2.s1.vi)[51] <- "var_name"
colnames(fss2.s1.vi)[52] <- "var_index"
fss2.s1.visd[, 51]<- fss2.s1.col
fss2.s1.visd[, 52]<- c(1:length(fss2.s1.col))
colnames(fss2.s1.visd)[51] <- "var_name"
colnames(fss2.s1.visd)[52] <- "var_index"

# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file = paste0("fss2_s1_vi_", timestamp, ".RData"))
save(fss2.s1.visd, file = paste0("fss2_s1_visd_", timestamp, ".RData"))
timestamp
load("fss2_s1_vi_20150720_1627.RData")
load("fss2_s1_visd_20150720_1627.RData")

#### s1 boxplot ####
fss2.s1.vi.l <- melt(fss2.s1.vi, id = c("var_name", "var_index"))

#png('varImpALL.png', width = 960, height = 960)
bymedian <- with(fss2.s1.vi.l, reorder(var_index, value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        ylab = "Variable index", xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
#dev.off()

#R2
1 - sum((fss2.s1$FSS_26Aug14 - predict(fss2.s1.rf))^2) / 
  sum((fss2.s1$FSS_26Aug14 - mean(fss2.s1$FSS_26Aug14))^2)
#0.5491

fss2.s1.vi.median <- cast(fss2.s1.vi.l, var_name + var_index ~ ., 
                          value ='value', median)
colnames(fss2.s1.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s1.vi.median <- fss2.s1.vi.median[with(fss2.s1.vi.median, order(-median)), ]

## Save the precursors to s2 with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi.median, file=paste0("fss2_s1_vi_median_", timestamp, ".RData"))
save(fss2.s1, file=paste0("fss2_s1_", timestamp, ".RData"))
timestamp
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

# --- #
##### Random Forest Step 2 ####
colnames(fss2.s2)

# mtry and ntree values 
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(nrow = ncol(fss2.s2) - 1, ncol = 50))
fss2.s2.visd <- data.frame(matrix(nrow = ncol(fss2.s2) - 1, ncol = 50))

fss2.s2.col <- colnames(fss2.s2)
fss2.s2.col <- fss2.s2.col[!(fss2.s2.col == "FSS_26Aug14")]

beg <- Sys.time()
set.seed(100)
for (i in 1:50) {
  fss2.s2.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s2, 
                             ntree = 1000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s2.vi[, i] <- importance(fss2.s2.rf, type = 1, conditional = TRUE)
  fss2.s2.visd[, i] <- fss2.s2.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s2.vi[, 51]<- fss2.s2.col
fss2.s2.vi[, 52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.vi)[51] <- "var_name"
colnames(fss2.s2.vi)[52] <- "var_index"
fss2.s2.visd[, 51]<- fss2.s2.col
fss2.s2.visd[, 52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.visd)[51] <- "var_name"
colnames(fss2.s2.visd)[52] <- "var_index"

# --- #
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s2.vi, file = paste0("fss2_s2_vi_", timestamp, ".RData"))
save(fss2.s2.visd, file = paste0("fss2_s2_visd_", timestamp, ".RData"))
timestamp
load("fss2_s2_vi_20141027_2010.RData")
load("fss2_s2_visd_20141027_2010.RData")

#### s2 boxplot ####

fss2.s2.vi.l <- melt(fss2.s2.vi, id = c("var_name", "var_index"))

# bymedian <- sort(sapply(fss2.s2.rm.rf.vi, median))
# index.merge <- data.frame('variable' = names(bymedian), 'index' = 1:11)
# fm <- melt(fss2.s2.rm.rf.vi)
# fm <- merge(fm, index.merge, by = 'variable', all.x = TRUE)

fss2.vars <- fss2.s1.vi[fss2.s1.vi$var_name %in% setdiff(names(fss2.s2),
                                                         "FSS_26Aug14"),]
fss2.s2.vi.l <- melt(fss2.vars, id = c("var_name", "var_index"))
png('varImp_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s2.vi.l, reorder(var_name, value, median))
par(yaxt = "n",mar = c(5, 8, 4, 5), cex = 2)
boxplot(value ~ bymedian, data = fss2.s2.vi.l,
        xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", 
        horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:100, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
dev.off()

png('varImpALL_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s2.vi.l, reorder(var_name, value, median))
par(yaxt="n", mar=c(5, 8, 4, 5), cex=2)
boxplot(value ~ bymedian, data = fss2.s2.vi.l,
        xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", 
        horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:14, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
dev.off()

#R2
1 - sum((fss2.s2$FSS_26Aug14 - predict(fss2.s2.rf))^2) / 
  sum((fss2.s2$FSS_26Aug14 - mean(fss2.s2$FSS_26Aug14))^2)
#0.5658064
#This corroborates the use of these variables by showing that with a smaller subest of variables
#we achieve essentially the same R2. When we look at the bottom 2/3 of the variables we get a much smaller R2.

#### Median df creation ####

fss2.s2.vi.median <- cast(fss2.s2.vi.l,var_name + var_index ~ ., 
                          value ='value', median)
colnames(fss2.s2.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s2.vi.median <- fss2.s2.vi.median[with(fss2.s2.vi.median, 
                                            order(-median)), ]

#### Partial Dependence Plots ####
for (i in 1:length(fss2.s2.vi$var_name)) {
  filename <- paste("partialPlot_", fss2.s2.vi$var_name[i], ".png", sep = "")
  png(filename, width = 960, height = 960)
  partialPlot(fss2.s2.rf, 
              fss2.s2, 
              x.var = fss2.s2.vi$var_name[i], 
              ylab = 'Mean FSS', 
              ylim = c(9, 18))
  dev.off()  
}
