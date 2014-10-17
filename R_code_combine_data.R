library(plyr)
library(RODBC)
library(SSN)
library(foreign)
options(stringsAsFactors = FALSE)

# -----------------------------------------------------------
# FUNCTIONS

# Function to calculate accumulated attributes at sites
siteaccum <- function(Edgedf, Sitesdf, EdgeVar, AEdgeVar, upratio, station, by.site, by.edge) {
  # get the cols
  dfx <- Sitesdf[,c(station, by.site, upratio)]
  dfy <- Edgedf[,c(by.edge, EdgeVar, AEdgeVar)]
  dfm <- merge(dfx, dfy, by.x=by.site, by.y=by.edge, all.x=TRUE)
  accum <- dfm[,AEdgeVar] - (dfm[,upratio] * dfm[,EdgeVar])
  return(accum)
}

# where NAs exist it will pull value from the other col, returns a df
replacena <- function(df, x, y){
  df[, x] <- ifelse(is.na(df[, x]), df[, y],df[, x])
  df[, y] <- ifelse(is.na(df[, y]), df[, x],df[, y])
  return(df)
  }

propor<- function(df, numer, denom, newfield){
  df[,newfield] <- (df[,numer] / df[,denom]) * 100
  return(df)
}

# -----------------------------------------------------------
# Read in data from access tables

indb <- "//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/Tables.mdb"
tablename1 <- "ssn_edges_table_final"
tablename2 <- "ssn_sites_table_final"
tablename3 <- "FSS_by_SVN"
tablename4 <- "tbl_HASLIDAR_Station_watershed"
tablename5 <- "tb_PPT_annual_avg_by_STATION_KEY"
tablename6 <- "tbl_POPRCA2010_by_RID"
tablename7 <- "var_proportion_table"
channel <-odbcConnectAccess2007(indb)
edgedf <- sqlFetch(channel, tablename1)
obs <- sqlFetch(channel, tablename2)
#stations.df <- sqlFetch(channel, tablename3)
fss <- sqlFetch(channel, tablename3)
haslidar <- sqlFetch(channel, tablename4)
ppt <- sqlFetch(channel, tablename5)
pop <- sqlFetch(channel, tablename6)
pvar <- sqlFetch(channel, tablename7)
close(channel)
rm(indb, tablename1, tablename2, tablename3, tablename4, tablename5, tablename6, tablename7, channel)
# -----------------------------------------------------------
# Clean up the data and accumulate some of the variables

colnames(obs)
colnames(edgedf)

ppt <- within(ppt, rm(OBJECTID))
pop <- within(pop, rm(OBJECTID))

edgedf <- merge(edgedf, pop, by="rid", all.x=TRUE)

# remove all the NAs in the accumulated fields except fishpres
edgedf[!(names(edgedf) %in%"fishpres")][is.na(edgedf[!(names(edgedf) %in%"fishpres")])] <- 0


# These are edgedf col that we want to delete from obs after the merge
edge.rm <- c("OBJECTID", "arcid", "from_node", 
             "to_node","HydroID", "GridID","NextDownID",
             "DrainID", "Shape_Length", "RCA_PI")

# we add the accumulated cols because they are going 
# to be added again in the siteaccum function
edge.rm <- c(edge.rm,names(edgedf)[grep('^A',names(edgedf))])

# We want to keep these ones though since we don't accumulate them at the site
edgekeep <- c("AROADX", "AROADLENRCAM", "AROADLENRSAM", "ASPLASH")
edge.rm <- edge.rm[!(edge.rm %in% edgekeep)]

# merge edge accumulations with the station and clean up
obs.a <- merge(obs, edgedf, by = 'rid', all.x = TRUE)
obs.a <- rename(obs.a, c('upDist.x' = 'upDist'))
obs.a <- within(obs.a, rm(OBJECTID.x,OBJECTID.y, upDist.y))
obs.a <- obs.a[, !colnames(obs.a) %in% edge.rm]

colnames(obs.a)


# These are obs.a col that we want to exclude from the site accumulation function
# all other col will be accumulated
noaccum <- c("rid", "OBJECTID",  "POURID", "STATION_KEY", "SITE_NAME",
             "HU_6_NAME", "HU_8_NAME", "HU_10_NAME", "HU_12_NAME",      
             "HU_08", "HU_10", "HU_12", "LONG_RAW",        
             "LAT_RAW", "NHDHigh", "NHDh_ReachCode", "NHDP21_ReachCode",
             "NHDP12_COMID", "NHDP21_COMID", "RESOLUTION", "ECO3_NAME", "VERSION",         
             "ratio", "locID", "netID", "pid", "upDist", 
             "afvArea", "fishpres", "Shape_Length","fishpres", 
             "ROADX", "ROADLENRCAM","ROADLENRSAM",
             "AROADX", "AROADLENRCAM", "AROADLENRSAM",
             "SINUMAP","SPLASH", "ASPLASH")

# All the col left are the ones we want to accumulate, except the noaccum.
edgevars <- colnames(obs.a[, !colnames(obs.a) %in% noaccum])

# This adds the A to get the accumulated counterpart
Aedgevars <- paste("A",edgevars, sep="")

# Run the loop to accumulate all the fields
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
rm(Aedgevars,edgevars, edge.rm, edgekeep, noaccum, i)

# -----------------------------------------------------------
# Read in SVNs to remove per Shannon Hubler comments (see Dealing with Low Counts_SH_4 8 14_RM.xlsx)
svn.rm <- read.csv('//deqhq1/TMDL/TMDL_WR/MidCoast/Data/BenthicMacros/Raw_From_Shannon/SVNs_to_Remove_2014_08_04.csv')

# Read in precip data
precip <- read.csv('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Precip/R_output_Precip_samples_2014-09-08.csv')
rm.col <- c("EXCLUDE","NewSample","Bug_RefApril13", "TMDL",
            "MIDCOAST", "Samples","Ref_Samples", "REF",
            "FLAG", "FLAG_REASON","REF_shubler", "Bug_RefMay05",
            "SVN2KEY", "Date", "Month_Sampled","Day_Sampled", 
            "Year_Sampled", "HabitatSampled","SamplingAgency",
            "SamplingProtocol", "Field_QAQC",
            "Lab_QAQC", "Project_name", "TS_May05", "FSS_May05", 
            "PREDATOR_Nov05_model", "PREDATOR_Nov05_score","PREDATOR_Nov05_Condition",
            "PREDATOR_outlier","PREDATOR_Integrated_Rpt_2010", "Bug_Count_RIV",
            "PREDATOR_reference_model", "EcoT20_FURR",
            "EMAP_REF", "FURR_Score","HDI_FURR", "SITE_NAME", "Shannon_Original", 
            "Shannon.First.Download")
precip <- precip[, !colnames(precip) %in% rm.col]

# Read in physical habitat data
phab.bugs <- read.csv('//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Station_Selection/Watershed_Char_phab_bugs_merge_SSN_FINAL.csv')

rm(rm.col)
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

precip2 <- merge(precip, ppt, by = 'STATION_KEY', all.x = TRUE)
rm(ppt)

precip3 <- precip2[,names(precip2)[!names(precip2) %in% names(phab.bugs)]]
precip <- cbind(data.frame('SVN' = precip2[,'SVN']),precip3)

comb <- merge(phab.bugs, precip, by = 'SVN', all = TRUE)
rm(precip, precip2, precip3)

# remove the cols from comb that are in slope
comb2 <- comb[,names(comb)[!names(comb) %in% names(slope.all)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']),comb2)
rm(comb2)

comb <- merge(comb, slope.all, by = 'STATION_KEY', all = TRUE)
rm(slope.all)

comb2 <- comb[,names(comb)[!names(comb) %in% names(obs.a)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']),comb2)
rm(comb2)

comb <- merge(obs.a, comb, by = 'STATION_KEY', all.y = TRUE)
rm(obs.a)

colnames(comb)
# -----------------------------------------------------------
# Calculate stream power by first pulling in the flow data from nhdplus v21
nhd <- read.dbf('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/NHDplus_21/EROMExtension/EROM_MA0001.DBF')
nhd.flow <- nhd[,c('Comid','Q0001A')]

comb <- merge(comb, nhd.flow, by.x = "NHDP21_COMID", by.y = 'Comid', all.x = TRUE)
comb$STRMPWR <- comb$XSLOPE_MAP * comb$Q0001A

rm(nhd.flow, nhd)
# -----------------------------------------------------------
#Pull in updated FSS values and remove the SVNs that have low counts
#names(comb)[grep('FSS',names(comb))] #FSS_May05
comb <- merge(comb, fss[,c('SVN','FSS_26Aug14')], by = 'SVN', all.x = TRUE)
comb <- within(comb, rm(FSS_May05))

#comb <- comb[!(comb$SVN %in% svn.rm),] # These SVNs aren't in the file
rm(fss, svn.rm)
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
# Calculate Year specfic disturbance # FIX

var <- paste0(disvar[i],"_",comb$YEAR[i])

comb[,disvar[i]] <- comb[,grep(var,names(comb))]

DISRCA_1YR <- 
DISRCA_3YR
DISRCA_10YR
DISRSA_1YR
DISRSA_3YR
DISRSA_10YR
ADISRCA_1YR
ADISRCA_3YR
ADISRCA_10YR
ADISRSA_1YR
ADISRSA_3YR
ADISRSA_10YR





# -----------------------------------------------------------
# Run the loop to calculate percentages, means, or densities
for (i in 1:nrow(pvar)) {
  #obs.a[,pvar[i,2]] <- (pvar[i,3] / pvar[i,4]) * 100
  obs.a <- propor(obs.a, pvar[i,3], pvar[i,4], pvar[i,2])
}

# -----------------------------------------------------------
# Fix the NA sample dates and related cols
# use datestr values

fix <- comb[c("DATE", "Date", "datestr", "Year_Sampled", "Month_Sampled", "YEAR")]

comb <- within(comb, rm(DATE, Date, Year_Sampled, Month_Sampled, YEAR))
comb$DATE <- as.POSIXlt(comb$datestr,format="%Y-%m-%d")
comb$YEAR <-  comb$DATE$year+1990
comb$MONTH <- comb$DATE$mon+1 # +1 because it is zero-indexed

rm(fix)
# -----------------------------------------------------------
# This generates a csv with each of the final col colnames.
# and number or NAs
# This file must be edited so 
# RF_Keep = 1 for col that will go to RF
# RF_KEEP = 0 = for col excluded from RF
# Resave the file as VarNames_RF.csv

na.list <- colSums(is.na(comb))
na.df <- t(as.data.frame(t(na.list),row.names = c("na.count")))
na.df <- data.frame(rownames(na.df),na.df,row.names = NULL)
colnames(na.df)[1]="var"
na.df$fss.rf_keep <- NA
na.df$pctfn.rf_keep <- NA
na.df$pctsa.rf_keep <- NA
na.df$pctsafn.rf_keep <- NA
na.df$oe.rf_keep <- NA

# -----------------------------------------------------------
# Write the files

write.csv(na.df, 'VarNACount.csv',row.names = TRUE)
write.csv(comb, 'ssn_RF_data.csv', row.names = FALSE)
