library(randomForest)
library(reshape)
library(plyr)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

#### Correlation Matrix Functions ####
# source
# http://stackoverflow.com/questions/15271103/how-to-modify-this-correlation-matrix-plot
#
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}

# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data ####
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.
vars[vars$var == 'STATION_KEY','fss2.rf_keep'] <- 1
vars.fss2.s1 <- vars[vars$fss2.rf_keep == 1,]

fss2.s1 <- bugs[,colnames(bugs) %in% vars.fss2.s1$var]

fss2.s1 <- within(fss2.s1, rm('STATION_KEY'))

# keep a copy with all the data
fss2.s1.na <- fss2.s1

#There are mostly NA's in these variables since we don't have lidar in most of the coast range
#and so aren't that useful. We'll remove them here.
fss2.s1 <- fss2.s1[,!(colnames(fss2.s1) %in% c("PSUSCEP4_LI","PSUSCEP5_LI",
                                               "PASUSCEP4_LI", "PASUSCEP5_LI"))]
colnames(fss2.s1)

# remove NAs in response variable
fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]
# remove any NAs - this has been taken care of prior to this but for good measure we'll leave it in
fss2.s1 <- data.frame(na.omit(fss2.s1))

# mtry value
mtry.fss2.s1 <- as.integer(((ncol(fss2.s1)-1) / 3),0)

# initialize the variable importance df to store the importance scores from each run
fss2.s1.vi <- data.frame(matrix(, nrow = ncol(fss2.s1)-1, ncol = 50))
fss2.s1.visd <- data.frame(matrix(, nrow = ncol(fss2.s1)-1, ncol = 50))

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
  fss2.s1.vi[,i] <- fss2.s1.rf$importance[,1]
  fss2.s1.visd[,i] <- fss2.s1.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s1.vi[,51]<- fss2.s1.col
fss2.s1.vi[,52]<-c(1:length(fss2.s1.col))
colnames(fss2.s1.vi)[51] <- "var_name"
colnames(fss2.s1.vi)[52] <- "var_index"
fss2.s1.visd[,51]<- fss2.s1.col
fss2.s1.visd[,52]<-c(1:length(fss2.s1.col))
colnames(fss2.s1.visd)[51] <- "var_name"
colnames(fss2.s1.visd)[52] <- "var_index"

# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file=paste0("fss2_s1_vi_",timestamp,".RData"))
save(fss2.s1.visd, file=paste0("fss2_s1_visd_",timestamp,".RData"))
timestamp
load("fss2_s1_vi_20141210_1515.RData")
load("fss2_s1_visd_20141210_1515.RData")

#### s1 boxplot ####
fss2.s1.vi.l <- melt(fss2.s1.vi, id=c("var_name","var_index"))

png('varImpALL.png', width = 960, height = 960)
bymedian <- with(fss2.s1.vi.l, reorder(var_index, value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        ylab = "Variable index", xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
dev.off()

#R2
1-sum((fss2.s1$FSS_26Aug14-predict(fss2.s1.rf))^2)/sum((fss2.s1$FSS_26Aug14-mean(fss2.s1$FSS_26Aug14))^2)
#0.5491

fss2.s1.vi.median <- cast(fss2.s1.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss2.s1.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s1.vi.median <- fss2.s1.vi.median[with(fss2.s1.vi.median, order(-median)), ]

## Save the precursors to s2 with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi.median, file=paste0("fss2_s1_vi_median_",timestamp,".RData"))
save(fss2.s1, file=paste0("fss2_s1_",timestamp,".RData"))
timestamp
load("fss2_s1_vi_median_20141210_1515.RData")
load("fss2_s1_20141210_1515.RData")

#### Variable selection ####
# Values drop off and then level out. Arbitrarily going with 50% of the variables.
# grab all variable names with median values > 0.4805920 = 50% of the data
fss2.s2.col <- fss2.s1.vi.median[fss2.s1.vi.median$median >= 0.4805920,][,1]
fss2.s2.col <- c("FSS_26Aug14",fss2.s2.col)
fss2.s2 <- fss2.s1[,colnames(fss2.s1) %in% fss2.s2.col]

#### Correlation plots ####
fss2.s2.col

#Method: Look for r<=0.2. If lower importance variable has r>=0.2 with any higher importance variables select the highest 
#importance variable unless that variable is already conceptually represented.

# Precip  "sum_1095_days" "PPT_1981_2010" "sum_365_days"  "sum_180_days"  "sum_60_days"
png('precip_cor.png')
pairs(fss2.s2[,fss2.s2.col[c(3,4,7,33,37)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()
# keep "sum_1095_days"

# Disturb [1] "PDISRSA_1YR"  "PADISRSA_1YR" "PDISRCA_1YR"  "PDISRCA_3YR"  "PADISRCA_1YR" "PDISRCA_10YR" "PADISRCA_3YR"
pairs(fss2.s2[,fss2.s2.col[c(9,10,13,21,41,45,47)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# everything is coorelated
# keep "PDISRSA_1YR",

# Lithology/soils
#[1] "PALITHERODRCA"  "PALITHERODRSA"  "PASILTRCA"      "PACLAYRCA"      "PASILT_CLAYRCA" "PASANDRCA"      "MAKFACTRCA"     "PSILTRCA"      
#[9] "PCLAYRCA"       "PLITHERODRSA"   "PLITHERODRCA"   "PSANDRCA"       "PSILT_CLAYRCA"  "MKFACTRCA"      "PALITHCOMPRCA"
pairs(fss2.s2[,fss2.s2.col[c(2,6,14,15,17,19,20,22,23,29,32,34,35,48,56)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# almost everything is coorelated
# Keep "PALITHERODRCA", "PASILTRCA", "PACLAYRCA"

# Ownership 
#  [1] "DAPOPRCA2010" "APOPRCA2010"  "POWNRCA_PRI"  "PAOWNRSA_PRI" "POPRCA2010"   "POWNRSA_PRI"  "PAOWNRCA_AGR" "POWNRCA_FED"  "PAOWNRSA_FED"
# [10] "PAOWNRCA_PRI" "POWNRSA_FED"  "PAOWNRCA_URB" "PAOWNRSA_AGR" "DAROADX"
pairs(fss2.s2[,fss2.s2.col[c(12,24,31,38,39,40,46,50,51,52,53,54,55,57)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# Keep"DAPOPRCA2010","POWNRCA_PRI","PAOWNRCA_AGR,"DAROADX"

#Susceptibility
# "PASUSCEP5_DE" "PASUSCEP4_DE" "PSUSCEP4_DE"  "PSUSCEP5_DE" 
pairs(fss2.s2[,fss2.s2.col[c(16,18,25,43)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
#All correlated. 
#Keep "PASUSCEP5_DE"

# Others/The rest
#[1] "STRMPWR"    "XSLOPE_MAP" "MIN_Z"      "LONG_RAW"   "LAT_RAW"    "upDist"     "afvArea"    "PATYPEF"    "ARCASQM"    "ARSASQM"   
#[11] "Q0001A"   
pairs(fss2.s2[,fss2.s2.col[-c(1,3,4,7,33,16,18,25,43,37,9,10,13,21,41,45,47,2,6,14,15,17,19,20,22,23,29,32,34,35,48,56,12,24,31,38,39,40,46,50,51,52,53,54,55,57)]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
# udist and Long coorelated
# keeep [1] "STRMPWR", "XSLOPE_MAP","MIN_Z","LAT_RAW"

keeps.s2 <- c("FSS_26Aug14",
              "sum_1095_days", 
              "PDISRSA_1YR",
              "PALITHERODRCA", "PASILTRCA", "PACLAYRCA",
              "DAPOPRCA2010","POWNRCA_PRI","PAOWNRCA_AGR","DAROADX",
              "PASUSCEP5_DE",
              "STRMPWR", "XSLOPE_MAP","MIN_Z","LAT_RAW")

pairs(fss2.s2[,fss2.s2.col[fss2.s2.col %in% keeps.s2[-1]]],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

#Further remove variables to reduce the influence of correlation on raising variable importance
fss2.s2 <- fss2.s2[,colnames(fss2.s2) %in% keeps.s2]
colnames(fss2.s2)


# remove any NAs
fss2.s2 <- data.frame(na.omit(fss2.s2))
#write.csv(fss2.s2, 'fss2_s2_data.csv')

# --- #
##### Random Forest Step 2 ####
colnames(fss2.s2)

# mtry and ntree values 
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 50))
fss2.s2.visd <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 50))

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
  fss2.s2.vi[,i] <- importance(fss2.s2.rf, type = 1, conditional = TRUE)
  fss2.s2.visd[,i] <- fss2.s2.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s2.vi[,51]<- fss2.s2.col
fss2.s2.vi[,52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.vi)[51] <- "var_name"
colnames(fss2.s2.vi)[52] <- "var_index"
fss2.s2.visd[,51]<- fss2.s2.col
fss2.s2.visd[,52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.visd)[51] <- "var_name"
colnames(fss2.s2.visd)[52] <- "var_index"

# --- #
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s2.vi, file=paste0("fss2_s2_vi_",timestamp,".RData"))
save(fss2.s2.visd, file=paste0("fss2_s2_visd_",timestamp,".RData"))
timestamp
load("fss2_s2_vi_20141027_2010.RData")
load("fss2_s2_visd_20141027_2010.RData")

#### s2 boxplot ####

fss2.s2.vi.l <- melt(fss2.s2.vi, id=c("var_name","var_index"))

# bymedian <- sort(sapply(fss2.s2.rm.rf.vi, median))
# index.merge <- data.frame('variable' = names(bymedian), 'index' = 1:11)
# fm <- melt(fss2.s2.rm.rf.vi)
# fm <- merge(fm, index.merge, by = 'variable', all.x = TRUE)

#png('varImpALL_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s2.vi.l, reorder(var_name, value, median))
par(yaxt="n",mar=c(5, 8, 4, 5))
boxplot(value ~ bymedian, data = fss2.s2.vi.l,
        xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:15, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
#dev.off()

#R2
1-sum((fss2.s2$FSS_26Aug14-predict(fss2.s2.rf))^2)/sum((fss2.s2$FSS_26Aug14-mean(fss2.s2$FSS_26Aug14))^2)
#0.5658064
#This corroborates the use of these variables by showing that with a smaller subest of variables
#we achieve essentially the same R2. When we look at the bottom 2/3 of the variables we get a much smaller R2.

#### Median df creation ####

fss2.s2.vi.median <- cast(fss2.s2.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss2.s2.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s2.vi.median <- fss2.s2.vi.median[with(fss2.s2.vi.median, order(-median)), ]

#### Partial Dependence Plots ####
for (i in 1:length(fss2.s2.vi$var_name)) {
  filename <- paste("partialPlot_", fss2.s2.vi$var_name[i], ".png",sep = "")
  png(filename, width = 960, height = 960)
  partialPlot(fss2.s2.rf, 
              fss2.s2, 
              x.var = fss2.s2.vi$var_name[i], 
              ylab = 'Mean FSS', 
              ylim = c(9,18))
  dev.off()  
}
