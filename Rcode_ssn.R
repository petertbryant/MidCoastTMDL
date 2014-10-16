


# library(raster)
# # import the edge shapefile
# edgepath <- "C:/WorkSpace/Biocriteria/WatershedCharaterization/SSN/LSN04/lsn.ssn/edges.shp"
# edge_shp <- shapefile(edgepath ,stringsAsFactors=FALSE)
# edgedf <- data.frame(edge_shp)
# rm(edge_shp, edgepath)

library(RODBC)
## READ DATA FROM ACCESS 2007
indb <- "C:/WorkSpace/Biocriteria/WatershedCharaterization/SSN/LSN04/Tables.mdb"
tablename <- "edge_table_final"
channel <-odbcConnectAccess2007(indb)
edgedf <- sqlFetch(channel, tablename)
close(channel)
rm(indb, tablename, channel)

library(SSN)
bugs <- importSSN("C:/WorkSpace/Biocriteria/WatershedCharaterization/SSN/LSN04/lsn.ssn", o.write=FALSE)
obs<- getSSNdata.frame(bugs, Name = "Obs")

colnames(edgedf)
colnames(obs)

#obs <- siteaccum(edgedf, obs, "RCASQM", "ARCASQM", "ratio", "rid", "rid")

dfx <- obs[,c("STATION_KE","rid", "ratio")]
dfy <- edgedf[,c("rid", "RCASQM", "ARCASQM")]
dfm <- merge(dfx, dfy, by.x="rid", by.y="rid", all.x=TRUE)
dfm$accum <- dfm$ARCASQM - (dfm$ratio * dfm$RCASQM)
dfm2 <- dfm[,c("STATION_KE","accum")]
obs2 <- merge(obs, dfm2, by = 'STATION_KE', all.x =TRUE)

# Function to calculate accumulated attributes at sites # THIS NEEDS WORK
siteaccum <- function(Edgedf, Sitesdf, EdgeVar, AEdgeVar, upDist, station, by.site, by.edge) {
  # get the cols
  dfx <- Sitesdf[,c(station, by.site, upDist)]
  dfy <- Edgedf[,c(by.edge, EdgeVar, AEdgeVar)]
  dfm <- merge(dfx, dfy, by.x=by.site, by.y=by.edge, all.x=TRUE)
  dfm$accum <- dfm[,AEdgeVar] - (dfm[,upDist] * dfm[,EdgeVar])
  #dfm2 <- dfm[,c(station,'accum')]
  #Sitesdf <- merge(Sitesdf, dfm2, by = station, all.x = TRUE)
  #return(Sitesdf)
  return(dfm$accum)
}

obs3 <- siteaccum(Edgedf = edgedf, 
                  Sitesdf = obs,
                  EdgeVar = "RCASQM", 
                  AEdgeVar = "ARCASQM", 
                  upDist = "ratio",
                  station = "STATION_KE",
                  by.site = "rid",
                  by.edge = "rid")



Aedgevars <- names(edgedf)[grep('^A',names(edgedf))]
edgevars <- gsub('^A','',Aedgevars)

obs2 <- obs
for (i in 1:length(edgevars)) {
  newcol <- siteaccum(Edgedf = edgedf, 
                      Sitesdf = obs,
                      EdgeVar = edgevars[i], 
                      AEdgeVar = Aedgevars[i], 
                      upDist = "ratio",
                      station = "STATION_KE",
                      by.site = "rid",
                      by.edge = "rid")
  obs2 <- cbind(obs2, newcol)
}

###################################
#######   random forests  ##########
#######################################
# --- EDIT THIS From another project.

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
