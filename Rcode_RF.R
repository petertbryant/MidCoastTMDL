library(party)
library(randomForest)
library(reshape)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

bugs <- arrange(bugs, STATION_KEY, desc(YEAR))
bugs <- bugs[!duplicated(bugs$STATION_KEY),]

# -----------------------------------------------------------
# Correlation Matrix Functions
# source
# http://stackoverflow.com/questions/15271103/how-to-modify-this-correlation-matrix-plot
# -----------------------------------------------------------
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
# FSS1 - Random forests including some of the physical habitat data
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss1.s1 <- vars[vars$fss1.rf_keep == 1,]
fss1.s1 <- bugs[,colnames(bugs) %in% c(vars.fss1.s1$var)]

# remove NAs in response variable
fss1.s1 <- fss1.s1[(!is.na(fss1.s1$FSS_26Aug14)),]

#remove NAs
#fss1.s1 <- data.frame(na.omit(fss1.s1))

# impute the NAs
set.seed(111)
fss1.s1.imputed <- rfImpute(FSS_26Aug14 ~ ., fss1.s1, ntree=2000, iter=3)


#Output for inclusion in chart for presentation
#write.csv(data.frame(variable = colnames(fss1.s1)),'fss1_s1_variable_categories.csv')
var.cat <- read.csv('fss1_s1_variable_categories.csv')
var.cat <- var.cat[!is.na(var.cat$Category),]
hbp <- barplot(summary(var.cat$Category), cex.names = 0.7, ylab = 'Variable Count', main = 'Count of Variable Categories')
text(x= hbp, y= summary(var.cat$Category)+3, labels=as.character(summary(var.cat$Category)), xpd=TRUE)

# mtry and ntree values 
mtry.fss1.s1 <- as.integer(((ncol(fss1.s1)-1) / 3),0)

colnames(fss1.s1.imputed)


# initialize the variable importance df
fss1.s1.vi <- data.frame(matrix(, nrow = ncol(fss1.s1)-1, ncol = 50))
fss1.s1.visd <- data.frame(matrix(, nrow = ncol(fss1.s1)-1, ncol = 50))

fss1.s1.col <- colnames(fss1.s1)
fss1.s1.col <- fss1.s1.col[!(fss1.s1.col == "FSS_26Aug14")]

# WARNING - Takes about 1 hour
beg <- Sys.time()
set.seed(42)
for (i in 1:50) {
  fss1.s1.rf <- randomForest(FSS_26Aug14 ~ ., 
                                data = fss1.s1.imputed, 
                                ntree = 2000, 
                                keep.forest = TRUE, 
                                importance = TRUE)
  fss1.s1.vi[,i] <- fss1.s1.rf$importance[,1]
  fss1.s1.visd[,i] <- fss1.s1.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss1.s1.vi[,51]<- fss1.s1.col
fss1.s1.vi[,52]<-c(1:length(fss1.s1.col))
colnames(fss1.s1.vi)[51] <- "var_name"
colnames(fss1.s1.vi)[52] <- "var_index"
fss1.s1.visd[,51]<- fss1.s1.col
fss1.s1.visd[,52]<-c(1:length(fss1.s1.col))
colnames(fss1.s1.visd)[51] <- "var_name"
colnames(fss1.s1.visd)[52] <- "var_index"



# ----------
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss1.s1.vi, file=paste0("fss1_s1_vi_",timestamp,".RData"))
save(fss1.s1.visd, file=paste0("fss1_s1_visd_",timestamp,".RData"))

#load("fss1_s1_v1_20141024_0542.RData")
#load("fss1_s1_v1sd_20141024_0542.RData")
# ----------

fss1.s1.vi.l <- melt(fss1.s1.vi, id=c("var_name","var_index"))
fss1.s1.visd.l <- melt(fss1.s1.visd, id=c("var_name","var_index"))


bymedian <- with(fss1.s1.vi.l, reorder(var_index, value, median))
boxplot(value ~ bymedian, data = fss1.s1.vi.l,
        ylab = "Variable index", xlab = "Importance", 
        varwidth = TRUE,
        col = "lightgray")


bymedian_vi <- with(fss1.s1.vi.l, reorder(var_index, -value, median))
bymedian_visd <- with(fss1.s1.visd.l, reorder(var_index, -value, median))


fss1.s1.vi.median <- cast(fss1.s1.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss1.s1.vi.median )[3] <- "median"

fss1.s1.visd.median <- cast(fss1.s1.visd.l,var_name + var_index ~ ., value ='value', median)
colnames(fss1.s1.visd.median )[3] <- "mediansd"

# sort the data so largest at the top
fss1.s1.vi.median <- fss1.s1.vi.median[with(fss1.s1.vi.median, order(-median)), ]
fss1.s1.visd.median <- fss1.s1.visd.median[with(fss1.s1.visd.median, order(-mediansd)), ]

# plot the vi
boxplot(value ~ bymedian_vi, data = fss1.s1.vi.l,
        xlab = "Variable index", ylab = "Importance", 
        varwidth = TRUE,
        col = "lightgray")

# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss2.s1 <- vars[vars$fss2.rf_keep == 1,]
fss2.s1 <- bugs[,colnames(bugs) %in% vars.fss2.s1$var]

# remove NAs in response variable
fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]

colnames(fss2.s1)

# mtry and ntree values 
mtry.fss2.s1 <- as.integer(((ncol(fss2.s1)-1) / 3),0)

# initialize the variable importance df
fss2.s1.vi <- data.frame(matrix(, nrow = ncol(fss2.s1)-1, ncol = 50))
fss2.s1.col <- colnames(fss2.s1)
fss2.s1.col <- fss2.s1.col[!(fss2.s1.col == "FSS_26Aug14")]

# WARNING - Takes a long time to run.
set.seed(42)
for (i in 1:50) {
  print(i)
  fss2.s1.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s1, controls = cforest_unbiased(ntree = 2000, mtry = mtry.fss2.s1))
  fss2.s1.vi[,i]<- varimp(fss2.s1.cf, conditional=FALSE)
}

test <- ctree(FSS_26Aug14 ~ ., data = fss2.s1)
plot(test)

# Add var names and index
fss2.s1.vi[,51]<- fss.col
fss2.s1.vi[,52]<-c(1:length(fss.col))
colnames(fss2.s1.vi)[51] <- "var_name"
colnames(fss2.s1.vi)[52] <- "var_index"

# ----------
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file=paste0("fss2_vi_s1_",timestamp,".RData"))

load("fss2_vi_s1_20141019_1451.RData")
# ----------

fss2.s1.vi.l <- melt(fss2.s1.vi, id=c("var_name","var_index"))

png('varImpALL.png', width = 960, height = 960)
bymedian <- with(fss2.s1.vi.l, reorder(var_index, value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        ylab = "Variable index", xlab = "Importance", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
dev.off()

fss2.s1.vi.median <- cast(fss2.s1.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss2.s1.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s1.vi.median <- fss2.s1.vi.median[with(fss2.s1.vi.median, order(-median)), ]

# Variable removal
# After the first 15 largest median values starts to flatten out so 
# we will take the top 30 to step 2.

# grab all variable names with median values > 0.5
fss2.s2.col <- fss2.s1.vi.median[fss2.s1.vi.median$median > 0.5,][,1]
fss2.s2.col <- c("FSS_26Aug14",fss2.s2.col)
fss2.s2 <- fss2.s1[,colnames(fss2.s1) %in% fss2.s2.col]

#since the variables are mostly station specific and not sample specific we should 
#subset the data frame
fss2.s2 <- bugs[,colnames(bugs) %in% fss2.s2.col]

#Further remove variables to reduce the influence of correlation on raising variable importance
fss2.s2.rm <- within(fss2.s2, rm('PALITHERODRSA','PPT_1981_2010','sum_365_days','sum_60_days','sum_180_days',
                                 'PDISRSA_1YR','PDISRCA_1YR','PDISRCA_3YR','PASILT_CLAYRCA','PASANDRCA',
                                 'PSILTRCA','PLITHERODRSA','PAOWNRSA_FED','POWNRCA_PRI','POWNRSA_FED',
                                 'PAOWNRCA_FED','PAOWNRCA_AGR','MAKFACTRCA','PLITHERODRCA'))

# remove any NAs
fss2.s2.rm <- data.frame(na.omit(fss2.s2.rm))

#write.csv(fss2.s2, 'fss2_s2_data.csv')

# -----------------------------------------------------------
# Correlation plots
colnames(fss2.s2)

# Precip
png('precip_cor.png')
pairs(fss2.s2[,c(4:6)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

# Disturb
pairs(fss2.s2[,c(14:16,27)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Lithology/soils
pairs(fss2.s2[,c(11:13,20:26)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Ownership
pairs(fss2.s2[,c(17:19,28:31)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Disturb and ownership
pairs(fss2.s2[,c(14:16,27,17:19,28:31)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

#Re-run the Random Forest on the subset of variables

#This is for cforest from party but since we can't get the partial dependencies we 
#aren't using this here now
# mtry and ntree values 
# mtry.fss2.s2.rm <- as.integer(((ncol(fss2.s2.rm)-1) / 3),0)
# set.seed(42)
# for (i in 1:50) {
#   print(i)
#   fss2.s2.rm.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s2.rm, controls = cforest_unbiased(ntree = 2000, mtry = mtry.fss2.s2.rm))
#   fss2.s2.rm.vi[,i]<- varimp(fss2.s2.rm.cf, conditional=FALSE)
# }
# 
# fss2.s2.rm.rf.vi <- data.frame(matrix(, ncol = 11, nrow = 50))
# names(fss2.s2.rm.rf.vi) <- names(fss2.s2.rm)[-6]

#Here we use randomForest so that we can use the object for the partial dependencies

library(randomForest)
# initialize the variable importance df
fss2.s2.rm.vi <- data.frame(matrix(, nrow = ncol(fss2.s2.rm)-1, ncol = 50))
set.seed(100)
for (i in 1:50) {
  fss2.s2.rm.rf <- randomForest(FSS_26Aug14 ~ ., 
                                data = fss2.s2.rm, 
                                ntree = 2000, 
                                keep.forest = TRUE, 
                                importance = TRUE)
  fss2.s2.rm.rf.vi[i,] <- importance(fss2.s2.rm.rf, conditional = TRUE)
}

View(fss2.s2.rm.rf.vi)
bymedian <- sort(sapply(fss2.s2.rm.rf.vi, median))
index.merge <- data.frame('variable' = names(bymedian), 'index' = 1:11)
fm <- melt(fss2.s2.rm.rf.vi)
fm <- merge(fm, index.merge, by = 'variable', all.x = TRUE)
png('varImp.png', width = 960, height = 960)
par(yaxt="n",mar=c(5, 8, 4, 5))
boxplot(value ~ index, data = fm,
        xlab = "Importance",
        varwidth = TRUE,
        col = "lightgray",
        horizontal = TRUE)
lablist.y<-names(bymedian)
axis(2, labels = FALSE)
text(y = 1:11, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
dev.off()


# ref.varimp <- importance(fss2.s2.rm.rf, conditional = TRUE)  
# fss2.s2.rm.rf$importance
# print(fss2.s2.rm.rf)
# plot(fss2.s2.rm.rf)
# varImpPlot(fss2.s2.rm.rf, main = "Variable Importance", pch = 19)

# STEP 2. Here we follow reccomendations by Strobl et al (2008) and Strobl et al (2009) 
# and calculate conditional variable importance. This reduces importance scores on 
# variables that get hight socres from step 1 becuase they highly correlated to 
# to other ones. Probably applies to some of the climate variables.
# We couldn't do this in step 1 becuase we had too many variables and NAs.
# Not enough computer memory. A smaller dataset is ideal.



# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 1))
fss2.s2.col <- fss2.s2.col[!(fss2.s2.col == "FSS_26Aug14")]

# mtry and ntree values 
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

set.seed(55)
fss2.s2.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s2, controls = cforest_unbiased(mtry = mtry.fss2.s2))
#fss2.s2.vi[,1] <- varimp(fss2.s2.cf, conditional=TRUE, threshold = 0.8)
test2 <- varimp(fss2.s2.cf, conditional=FALSE)

library(randomForest)
fss2.s2.rf <- randomForest(FSS_26Aug14 ~ ., 
                           data = fss2.s2, 
                           ntree = 2000, 
                           keep.forest = TRUE, 
                           importance = TRUE)

SedImpairedSites <- c(21792, 21842, 25297, 26816, 26818, 26822, 26964, 
                      29906, 30403, 33320, 33327, 33329, 33333, 33361,
                      33417, 33417, 33418, 34660, 34665, 34695, 37165)

mean.data <- data.frame('Imp' = c('Sed Impaired', 
                                       'Reference'),
                        'mean_FSS' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'FSS_26Aug14']),
                                            mean(bugs[bugs$REF == 'Y','FSS_26Aug14'])),
                        'mean_sum1095days' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'sum_1095_days']),
                                                    mean(bugs[bugs$REF == 'Y','sum_1095_days'])),
                        'mean_PALITHERODRCA' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'PALITHERODRCA']),
                                                 mean(bugs[bugs$REF == 'Y','PALITHERODRCA'])),
                        'mean_PADISRSA_1YR' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'PADISRSA_1YR']),
                                                 mean(bugs[bugs$REF == 'Y','PADISRSA_1YR'])),
                        'mean_PDISRSA_1YR' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'PDISRSA_1YR']),
                                                mean(bugs[bugs$REF == 'Y','PDISRSA_1YR'])),
                        'mean_XSLOPE_MAP' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'XSLOPE_MAP']),
                                               mean(bugs[bugs$REF == 'Y','XSLOPE_MAP'])),
                        'mean_MIN_Z' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'MIN_Z']),
                                              mean(bugs[bugs$REF == 'Y','MIN_Z'])),
                        'mean_STRMPWR' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'STRMPWR'],na.rm = TRUE),
                                         mean(bugs[bugs$REF == 'Y','STRMPWR'])),
                        'mean_PASILTRCA' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'PASILTRCA']),
                                           mean(bugs[bugs$REF == 'Y','PASILTRCA'])),
                        'mean_POWNRCA_PRI' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'POWNRCA_PRI']),
                                             mean(bugs[bugs$REF == 'Y','POWNRCA_PRI'])),
                        'mean_PAOWNRCA_FED' = c(mean(bugs[bugs$STATION_KEY %in% SedImpairedSites, 'PAOWNRCA_FED']),
                                               mean(bugs[bugs$REF == 'Y','PAOWNRCA_FED']))
                        )

png('pardep_sum1095days.png', width = 1200, height = 800)
nf <- layout(mat = matrix(c(1,2,0,3),2,2, byrow=TRUE),  width = c(0.5,3), height = c(2,0.5))
par(mar = c(2,4,2,2), cex = 1.5)
boxplot(fss2.s2$FSS_26Aug14, ylim = c(0,20), outline = TRUE, frame = FALSE, bxp = 0.5, axes = FALSE)
partialPlot(fss2.s2.rm.rf, fss2.s2, x.var = 'sum_1095_days', ylab = 'Mean FSS', ylim = c(0,20), xlim = c(3000,10000))
boxplot(fss2.s2$sum_1095_days, ylim = c(3000,10000), bxp = 0.5, outline = TRUE, frame = FALSE, axes = FALSE, horizontal = TRUE)
legend(7900,25,mean.data$Imp,pch = c(19,17))
dev.off()
#Interesting

partialPlot(fss2.s2.rm.rf, fss2.s2, x.var = 'PALITHERODRCA', ylab = 'Mean FSS', ylim = c(7,25), xlim = c(30,100))
points(mean.data$mean_PALITHERODRCA, mean.data$mean_FSS, pch = c(19, 17))
#Interesting

partialPlot(fss2.s2.rf, fss2.s2, x.var = 'PADISRSA_1YR', ylab = 'Mean FSS', ylim = c(7,25), xlim = c(0,15))
points(mean.data$mean_PADISRSA_1YR, mean.data$mean_FSS, pch = c(19, 17))
#not interesting

partialPlot(fss2.s2.rf, fss2.s2, x.var = 'PDISRSA_1YR', ylab = 'Mean FSS', ylim = c(7,25), xlim = c(0, 40))
points(mean.data$mean_PDISRSA_1YR, mean.data$mean_FSS, pch = c(19, 17))
#Subtle

partialPlot(fss2.s2.rf, fss2.s2, x.var = 'XSLOPE_MAP', ylab = 'Mean FSS', ylim = c(7,25), xlim = c(0,10))
points(mean.data$mean_XSLOPE_MAP, mean.data$mean_FSS, pch = c(19, 17))
#ref and impaired have about the same mean xslope_map

partialPlot(fss2.s2.rf, fss2.s2, x.var = "MIN_Z", ylab = 'Mean FSS', ylim = c(7, 25), xlim = c(0,600))
points(mz$mean_MIN_Z, mz$mean_FSS, pch = c(19, 17))
#kind of interesting

partialPlot(fss2.s2.rf, fss2.s2, x.var = "STRMPWR", ylab = 'Mean FSS', ylim = c(7,25), xlim = c(0,400))
points(mean.data$mean_STRMPWR, mean.data$mean_FSS, pch = c(19, 17))
#interesting

partialPlot(fss2.s2.rf, fss2.s2, x.var = "PASILTRCA", ylab = 'Mean FSS', ylim = c(7,25))
points(mean.data$mean_PASILTRCA, mean.data$mean_FSS, pch = c(19, 17))
#line is interesting but ref and impaired plot in straight line section

partialPlot(fss2.s2.rf, fss2.s2, x.var = "POWNRCA_PRI", ylab = 'Mean FSS', ylim = c(7,25))
points(mean.data$mean_POWNRCA_PRI, mean.data$mean_FSS, pch = c(19, 17))
#boring - straight line

partialPlot(fss2.s2.rf, fss2.s2, x.var = "PAOWNRCA_FED", ylab = 'Mean FSS', ylim = c(7,25))
points(mean.data$mean_PAOWNRCA_FED, mean.data$mean_FSS, pch = c(19, 17))
#boring - straight line

ref.varimp <- importance(fss2.s2.rf, conditional = TRUE)  
fss2.s2.rf$importance
print(fss2.s2.rf)
plot(fss2.s2.rf)

png('varImp.png', width = 1200, height = 960)
par(cex = 2)
varImpPlot(fss2.s2.rf, main = "Variable Importance", pch = 19)
dev.off()
