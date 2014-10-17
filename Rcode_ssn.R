


library(party)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

# -----------------------------------------------------------
# FSS - Random forests
vars.fss <- vars[vars$fss.rf_keep == 1,]
fss <- bugs[,colnames(bugs) %in% vars.fss$var]

# remove NAs in response variable
fss <- fss[(!is.na(fss$FSS_26Aug14)),]
fss.na <- na.omit(fss)

# factor chr vars
fss$ECO3_NAME <- factor(fss$ECO3_NAME)
fss$fishpres <- factor(fss$fishpres)
fss$YEAR <- factor(fss$YEAR)

set.seed(99998)

fss.cf <- cforest(FSS_26Aug14 ~ ., data = fss,controls = cforest_unbiased(ntree = 50))

fss.ct <- ctree(FSS_26Aug14 ~ ., data = fss, controls = ctree_control(maxsurrogate = 3))

plot(fss.ct)


set.seed(99998)
varimp(fss.cf)

plot(fss.cf)


# -----------------------------------------------------------
library(randomForest)
colnames(ref.cal) 
dim(ref.cal)          

#run with all predictors
#ref.cal.2<-ref.cal[,c(7,8,11,13:18,22,28,46:49,51)]
ref.cal.2<-ref.cal[,c(7,8,11,13,15,17,18,22,46:48,49,51)]        #reduce predictors--RF can't handle categorical with >32 categories

#"LONG"           "LAT"           "ECO3_NAME"      "Elev_FT"        "Map_Slope"      "Precip_mm"      "Temp_CX10"     
#"Strm_Power_NHD" "BASIN_NA_1"     "FSP_EROD_P"     "FSP_ER"        "MAFLOWU"        "SLOPE"          "AREAWTMAP"     
#"FSP_ER_40"       "fss.trans"  

colnames(ref.cal.2)
head(ref.cal.2)

#random forests modeling
ref.cal.rf <- randomForest(fss.trans ~ ., data=ref.cal.2, importance=TRUE,  keep.forest=TRUE)
ref.cal.rf  #---15 preds = 34.4  % variance  

save(ref.cal.rf,file = "ref.cal_ranfor.RData") 
print(ref.cal.rf) 

#which variables are most influential on the RF model?                  
ref.varimp <- importance(ref.cal.rf, conditional = TRUE)  
ref.cal.rf$importance
print(ref.cal.rf)
plot(ref.cal.rf)
varImpPlot(ref.cal.rf)   
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
