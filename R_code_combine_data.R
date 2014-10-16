library(plyr)
library(RODBC)
library(SSN)
library(foreign)
options(stringsAsFactors = FALSE)

# -----------------------------------------------------------
# READ DATA FROM ACCESS 2007
indb <- "//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/Tables.mdb"
tablename1 <- "ssn_edges_table_final"
tablename2 <- "ssn_sites_table_final"
tablename3 <- "stations_table"
tablename4 <- "FSS_by_SVN"
tablename5 <- "tbl_HASLIDAR_Station_watershed"
channel <-odbcConnectAccess2007(indb)
edgedf <- sqlFetch(channel, tablename1)
obs <- sqlFetch(channel, tablename2)
stations.df <- sqlFetch(channel, tablename3)
fss <- sqlFetch(channel, tablename4)
haslidar <- sqlFetch(channel, tablename5)
close(channel)
rm(indb, tablename1, tablename2, tablename3, tablename4, tablename5, channel)

#read in the SSN object in order to get the station data
#bugs <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/lsn.ssn", o.write=FALSE)
#obs<- getSSNdata.frame(bugs, Name = "Obs")
#rm(bugs)

# -----------------------------------------------------------
# Function to calculate accumulated attributes at sites
siteaccum <- function(Edgedf, Sitesdf, EdgeVar, AEdgeVar, upratio, station, by.site, by.edge) {
  # get the cols
  dfx <- Sitesdf[,c(station, by.site, upratio)]
  dfy <- Edgedf[,c(by.edge, EdgeVar, AEdgeVar)]
  dfm <- merge(dfx, dfy, by.x=by.site, by.y=by.edge, all.x=TRUE)
  accum <- dfm[,AEdgeVar] - (dfm[,upratio] * dfm[,EdgeVar])
  return(accum)
}

# -----------------------------------------------------------

colnames(obs)
colnames(edgedf)

# These are edgedf col that we want to delete from obs after the merge
edgedel <- c("OBJECTID", "arcid", "from_node", 
              "to_node","HydroID", "GridID","NextDownID",
              "DrainID", "Shape_Length","upDist", "RCA_PI")

# we add the accumulated cols because they are going 
# to be added again in the siteaccum function
edgedel <- c(edgedel,names(edgedf)[grep('^A',names(edgedf))])

# We want to keep these ones though since we don't accumulate them at the site
edgekeep <- c("AROADX", "AROADLENRCAM", "AROADLENRSAM", "ASPLASH")
edgedel <- edgedel[!(edgedel %in% edgekeep)]

# merge edge accumulations with the station and clean up
obs.a <- merge(obs, edgedf, by = 'rid', all.x = TRUE)
obs.a <- rename(obs.a, c('upDist.x' = 'upDist'))
obs.a <- within(obs.a, rm(OBJECTID.x,OBJECTID.y, upDist.y))
obs.a <- obs.a[, !colnames(obs.a) %in% edgedel]

colnames(obs.a)

# These are obs.a col that we want to exclude from the site accumulation function
# all other col will be accumulated
noaccum <- c("rid", "OBJECTID",  "POURID", "STATION_KEY", "SITE_NAME",
             "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", "HU_12_NAME",      
             "HU_08", "HU_10", "HU_12", "LONG_RAW",        
             "LAT_RAW", "NHDHigh", "NHDh_ReachCode", "NHDP21_ReachCode",
             "NHDP12_COMID", "RESOLUTION", "ECO3_NAME", "VERSION",         
             "ratio", "locID", "netID", "pid", "upDist", 
             "afvArea", "fishpres", "Shape_Length","fishpres", 
             "ROADX", "ROADLENRCAM","ROADLENRSAM",
             "AROADX", "AROADLENRCAM", "AROADLENRSAM",
             "SINUMAP","SPLASH", "ASPLASH")

# All the col left are the ones we want to accumulate, except the noaccum.
edgevars <- colnames(obs.a[, !colnames(obs.a) %in% noaccum])

# This adds the A to get the accumulated counterpart
Aedgevars <- paste("A",edgevars, sep="")

for (i in 1:length(edgevars)) {
  obs.a[,Aedgevars[i]] <- siteaccum(Edgedf = edgedf, 
                                   Sitesdf = obs.a,
                                   EdgeVar = edgevars[i], 
                                   AEdgeVar = Aedgevars[i], 
                                   upratio = "ratio",
                                   station = "STATION_KEY",
                                   by.site = "rid",
                                   by.edge = "rid")
}

colnames(obs.a)
rm(Aedgevars,edgevars, edgedel, edgekeep, noaccum, i)
# -----------------------------------------------------------

# precip data
precip <- read.csv('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Precip/R_output_Precip_samples_2014-09-08.csv')

# physical habitat data
phab.bugs <- read.csv('//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Station_Selection/Watershed_Char_phab_bugs_merge_SSN_FINAL.csv')

# -----------------------------------------------------------
# Pull in the slope data, fix the col names, and append it together

slope <- read.csv('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/SLOPES/Final/slopesmerge.csv')
ryan.slope <- read.csv('//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/SLOPES/Final/slopes_ryan.txt')
ryan.slope <- rename(ryan.slope, c('SLOPE_AVG_MAP' = 'XSLOPE_MAP', 'Z_Min' = 'MIN_Z', 'Z_Max' = 'MAX_Z'))

ryan.slope <- cbind(ryan.slope, data.frame("Source" = rep(NA, nrow(ryan.slope)), 
                                           "CONFIDENCE" = rep(NA, nrow(ryan.slope)), 
                                           "RchLenFin" = rep(NA, nrow(ryan.slope)),
                                           "XSLOPE" = rep(NA, nrow(ryan.slope))))

ryan.slope <- within(ryan.slope, rm(OBJECTID, TYPE, Shape_Length))
slope <- rename(slope, c("SLOPE_AVG" = 'XSLOPE', "SLOPE_AVG_MAP" = 'XSLOPE_MAP'))
slope.all <- rbind(slope, ryan.slope)
slope.all <- slope.all[!duplicated(slope.all$STATION_KEY),]
rm(slope, ryan.slope)

# -----------------------------------------------------------
# Delete?

# name.resolve <- data.frame('phab.bugs.names' = sort(names(phab.bugs)),
#                            'common.names' = sort(names(phab.bugs)))
# 
# precip.resolve <- data.frame('precip.names' = sort(names(precip)),
#                              'common.names' = sort(names(precip)))
# 
# name.resolve <- merge(name.resolve, precip.resolve, by = 'common.names', all = TRUE)
# 
# slope.resolve <- data.frame('slope.names' = sort(names(slope)),
#                             'common.names' = sort(names(slope)))
# 
# name.resolve <- merge(name.resolve, slope.resolve, by = 'common.names', all = TRUE)
# 
# all(precip$SVN %in% phab.bugs$SVN) #TRUE
# 
# all(sort(chars$PREDATOR_Nov05_score) == sort(phab.bugs[phab.bugs$SVN %in% chars$SVN,'PREDATOR_Nov05_score'])) #TRUE
# all(sort(precip$PREDATOR_Nov05_score) == sort(phab.bugs[phab.bugs$SVN %in% precip$SVN,'PREDATOR_Nov05_score'])) #TRUE
# 
# all(sort(precip$TS_May05) == sort(phab.bugs[phab.bugs$SVN %in% precip$SVN,'TS_May05'])) #TRUE

# -----------------------------------------------------------
# Put the dfs together without duplicating columns

precip2 <- precip[,names(precip)[!names(precip) %in% names(phab.bugs)]]
precip <- cbind(data.frame('SVN' = precip[,'SVN']),precip2)

comb <- merge(phab.bugs, precip, by = 'SVN', all = TRUE)
rm(precip, precip2)

comb2 <- comb[,names(comb)[!names(comb) %in% names(slope.all)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']),comb2)
rm(comb2)

comb <- merge(comb, slope.all, by = 'STATION_KEY', all = TRUE)
rm(slope.all)

comb2 <- comb[,names(comb)[!names(comb) %in% names(obs.a)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']),comb2)
rm(comb2)

comb <- merge(comb, obs.a, by = 'STATION_KEY', all.y = TRUE)
rm(obs.a)

comb2 <- comb[,names(comb)[!names(comb) %in% names(stations.df)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']), comb2)
rm(comb2)

comb <- merge(comb, stations.df, by = 'STATION_KEY', all.x = TRUE)
rm(stations.df)

colnames(comb)
# -----------------------------------------------------------
# Calculate stream power by first pulling in the flow data from nhdplus v21
nhd <- read.dbf('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/NHDplus_21/EROMExtension/EROM_MA0001.DBF')
nhd.flow <- nhd[,c('Comid','Q0001A')]

comb <- merge(comb, nhd.flow, by.x = "NHDP21_COMID", by.y = 'Comid', all.x = TRUE)
comb$STRMPWR <- comb$XSLOPE_MAP * comb$Q0001A

rm(nhd.flow, nhd)
# -----------------------------------------------------------
#Pull in updated FSS values
#names(comb)[grep('FSS',names(comb))] #FSS_May05
comb <- merge(comb, fss[,c('SVN','FSS_26Aug14')], by = 'SVN', all.x = TRUE)
comb <- within(comb, rm(FSS_May05))
rm(fss)

# -----------------------------------------------------------
# Set SUSCEP_LI data to NA if there is less than 100% LiDAR coverage for that watershed
has.lidar.id <- haslidar[haslidar$HASLIDAR =="No", ]
has.lidar.id <- unique(has.lidar.id$STATION_KE)

changecol <- c("SUSCEP1_LI","SUSCEP2_LI","SUSCEP3_LI","SUSCEP4_LI","SUSCEP5_LI",
               "ASUSCEP1_LI","ASUSCEP2_LI","ASUSCEP3_LI","ASUSCEP4_LI","ASUSCEP5_LI")
for (i in 1:length(changecol)) {
  comb[,changecol[i]] <- ifelse(comb$STATION_KEY %in% has.lidar.id,NA,comb[,changecol[i]])
}
rm(haslidar, has.lidar.id, changecol, i)

# -----------------------------------------------------------
# Delete?

na.list <- apply((apply(comb, 2, is.na)),2,table)
na.df <- data.frame('col' = names(na.list))
na.df$na.count <- NA
for (i in 1:length(na.list)) {
  if(!is.na((na.list[i][[1]][2]))) {
    na.df$na.count[i] <- na.list[i][[1]][2]
  }
}
na.df$na.count <- ifelse(is.na(na.df$na.count),0,na.df$na.count)
write.csv(na.df, 'VarNACount.csv')

# -----------------------------------------------------------
rf.vars <- read.csv('Var_RF.csv')
ncol(comb[,rf.vars[rf.vars$RF_Keep == 1,'col']])

colnames(comb)

