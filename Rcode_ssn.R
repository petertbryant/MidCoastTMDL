

library(party)
library(reshape)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

# -----------------------------------------------------------
# FSS2 - Random forests excluding the physical habitat data

# -----------------------------------------------------------
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
  fss2.s1.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s1, controls = cforest_unbiased(ntree = 2000, mtry = mtry.fss))
  fss2.s1.vi[,i]<- varimp(fss2.s1.cf, conditional=FALSE)
}

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

bymedian <- with(fss2.s1.vi.l, reorder(var_index, -value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        xlab = "Variable index", ylab = "Importance", 
        varwidth = TRUE,
        col = "lightgray")

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

# remove any NAs
fss2.s2 <- data.frame(na.omit(fss2.s2))

# -----------------------------------------------------------
# STEP 2. Here we follow reccomendations by Strobl et al (2008) and Strobl et al (2009) 
# and calculate conditional variable importance. This reduces importance scores on 
# variables that get hight socres from step 1 becuase they highly correlated to 
# to other ones. Probably applies to some of the climate variables.
# We couldn't do this in step 1 becuase we had too many variables and NAs.
# Not enough computer memory. A smaller dataset is ideal.

# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 50))
fss2.s2.col <- colnames(fss2.s2)
fss2.s2.col <- fss2.s2.col[!(fss2.s2.col == "FSS_26Aug14")]

# mtry and ntree values 
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

set.seed(84)
#for (i in 1:50) {
  #print(i)
 #fss2.s2.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s2,controls = cforest_unbiased(ntree = 2000, mtry = mtry.fss2.s2))
#fss2.s2.vi[,i]<- varimp(fss2.s2.cf, conditional=TRUE)
#}

fss2.s2.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s2,controls = cforest_unbiased(mtry = mtry.fss2.s2))
test <- varimp(fss2.s2.cf, conditional=TRUE)

# -----------------------------------------------------------

fss.ct <- ctree(FSS_26Aug14 ~ ., data = fss, controls = ctree_control(maxsurrogate = 3))
plot(fss.ct)


# -----------------------------------------------------------
library(randomForest)
         

#random forests modeling
fss.rf <- randomForest(FSS_26Aug14 ~ ., data=fss.na, importance=TRUE,  keep.forest=TRUE)
fss.rf  #---15 preds = 34.4  % variance  

save(fss.rf = "fss.na.rf.RData") 
print(ref.cal.rf) 

#which variables are most influential on the RF model?                  
ref.varimp <- importance(fss.rf, conditional = TRUE)  
fss.rf$importance
print(fss.rf)
plot(fss.rf)
varImpPlot(fss.rf)   
#Influential Predictors, in order of strength
# 1) Stream Power, 
# 2) Flow, 
# 3) Precip(areaW), 
# 4) Elev, % Erod (watershed), Precip (point), Ecoregion


#make predictions for Calibration sites
ref.val.2 <- ref.val[,c(7,8,11,13:18,22,28,46:49,51)] # all predictors
dim(ref.val.2)

ref.val.pred<-predict(ref.cal.rf, newdata=ref.val.2)

#RMSE
library(hydroGOF)
RF.rmse<-rmse(sim=ref.val.pred, obs=ref.val.2$fss.trans)  
RF.rmse  
#rmse = 0.288 ; rmse from original full ref dataset, reduced predictors = 0.281 ---> not a major loss of performance with smaller ref dataset
#in log10 scale--> untransform
10^RF.rmse      # RMSE = 1.94 FSS, or 1.9%


###### RF down to key predictors
ref.cal.2<-ref.cal[,c( 13, 17, 22, 48, 51)]	      #final predictors
colnames(ref.cal.2)
ref.cal.rf2 <- randomForest(fss.trans ~ ., data=ref.cal.2, importance=TRUE,  keep.forest=TRUE)
ref.cal.rf2   

#validation predictions
ref.val.2red <- ref.val[,c( 13, 17, 22, 48, 51)]
dim(ref.val.2red)

ref.val.pred.red<-predict(ref.cal.rf2, newdata=ref.val.2red)

#RMSE
RF_red.rmse<-rmse(sim=ref.val.pred.red, obs=ref.val.2red$fss.trans)  
RF_red.rmse  
#rmse = 0.301
#in log10 scale--> untransform
10^RF.rmse      # RMSE = 1.999,  or 2.0%

# -----------------------------------------------------------
# SSN
library(SSN)
bugs <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/lsn.ssn", o.write=FALSE)
obs<- getSSNdata.frame(bugs, Name = "Obs")
