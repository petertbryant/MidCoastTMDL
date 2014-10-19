


library(party)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

# -----------------------------------------------------------
# FSS - Random forests
vars.fss <- vars[vars$fss2.rf_keep == 1,]
fss <- bugs[,colnames(bugs) %in% vars.fss$var]

# remove NAs in response variable
fss <- fss[(!is.na(fss$FSS_26Aug14)),]
fss.na <- data.frame(na.omit(fss))

colnames(fss.na)

# factor chr vars
#fss.na$YEAR <- factor(fss.na$YEAR)

set.seed(42)
fss.cf <- cforest(FSS_26Aug14 ~ ., data = fss,controls = cforest_unbiased(ntree = 50, stump = TRUE))
fss.ct <- ctree(FSS_26Aug14 ~ ., data = fss, controls = ctree_control(maxsurrogate = 3))

plot(fss.ct)

df <- data.frame(sort(varimp(fss.cf, conditional=FALSE),decreasing = TRUE))

set.seed(99998)

plot(fss.cf)

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
