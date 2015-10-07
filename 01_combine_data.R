#############################################################
# This script reads in various tables, cleans them, and
# makes numerous calculations to preperae the data for the
# random forest and ssn analysis. This file must be run in
# 32-bit R.

# Ryan Michie & Peter Bryant
#############################################################

library(plyr)
library(RODBC)
library(SSN)
library(foreign)
library(raster)
library(stringr)
options(stringsAsFactors = FALSE)
options("scipen"=100)

# -----------------------------------------------------------
# Read in data 

indb <- "//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/Tables.mdb"
tablename1 <- "ssn_edges_table_final"
tablename2 <- "ssn_sites_table_final"
tablename3 <- "FSS_by_SVN"
tablename4 <- "PPRCA_PPRSA_Disturbance"
#tablename4 <- "tbl_HASLIDAR_Station_watershed"
tablename5 <- "tb_PPT_annual_avg_by_STATION_KEY"
tablename6 <- "tbl_POP_Dasy"
tablename7 <- "PPRCA_PPRSA_ZonalStats"
tablename8 <- 'FPA_Streams'
channel <-odbcConnectAccess2007(indb)
edgedf <- sqlFetch(channel, tablename1)
obs <- sqlFetch(channel, tablename2)
fss <- sqlFetch(channel, tablename3)
dis <- sqlFetch(channel, tablename4)
#haslidar <- sqlFetch(channel, tablename4)
ppt <- sqlFetch(channel, tablename5)
pop <- sqlFetch(channel, tablename6)
pp.zstats <- sqlFetch(channel, tablename7)
fpa <- sqlFetch(channel, tablename8)
close(channel)
rm(indb, tablename1, tablename2, tablename3, tablename4, tablename5, tablename6, 
   tablename7, channel, tablename8)

# Read in physical habitat data
phab.bugs <- read.csv('//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/
                      Watershed_Characteristics/Station_Selection/
                      Watershed_Char_phab_bugs_merge_SSN_FINAL.csv')

# Read in table describing how to calculate proportions
pvar <- read.csv("var_proportion_table.csv")

# Read in SVNs to remove per Shannon Hubler comments (see Dealing with Low Counts_SH_4 8 14_RM.xlsx)
svn.rm <- read.csv('//deqhq1/TMDL/TMDL_WR/MidCoast/Data/BenthicMacros/
                   Raw_From_Shannon/SVNs_to_Remove_2014_08_04.csv')

# Read in precip data
precip <- read.csv('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/
                   Watershed_Characteristics/Precip/
                   R_output_Precip_samples_2014-09-08.csv')

# Read in NHD ComIDs. We use the COMIDs to map to NHD catchment. This layer was created by doing a spatial join of NHDPlus21 Catchments with the updated obs layer that 
# resolved issues with snapping the obs to LSN05 vs LSN04. Also the flow values and the cumulative areas.
nhd <- read.dbf('C:/users/pbryant/desktop/midcoasttmdl-gis/
                lsn05_watersheds/NHD21_obs_up.dbf', as.is = TRUE)
nhd.flow <- read.dbf('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/
                     Watershed_Characteristics/NHDplus_21/
                     EROMExtension/EROM_MA0001.DBF')
nhd.area <- read.dbf('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/
                     Watershed_Characteristics/NHDplus_21/
                     Attributes/CumulativeArea.dbf')

#Read in slope files
slope <- read.csv('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/
                  Watershed_Characteristics/SLOPES/Final/slopesmerge.csv')
ryan.slope <- read.csv('//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/
                       Watershed_Characteristics/SLOPES/Final/slopes_ryan.txt')

# -----------------------------------------------------------
# Accumulate variables to each observation. This involves two steps.
# 1. The first step relies on the watershed characteristics being calculated to each edge/catchment rid
# and accumulated downstream to each edge/catchment rid. 
# This step deletes the most downstream edge/catchment rid.
# 2. The second step relies on deriving pour point watershed areas within each catchment area (i.e. clipped whole pour point watersheds to the catchment polygon).
# This step takes the upstream accumulated value and adds the pour point catchment area to get the adjusted accumulated pour point watershed area
# accurately represented. 

# Note: This method was employed becuase we had already calculated and accumulated all the values in the old method. It may 
# be worth considering just tabulating all the variables to the whole pour point watershed instead of doing this subtraction and addition business.
# You would still need to clip and tabulate variables to each pour point watershed clipped to the catchment.

#Attach the pop
pop <- rename(pop, c('POPARCA' = 'APOPRCA'))
edgedf <- merge(edgedf, pop[,c('rid','POPRCA','APOPRCA')], by="rid", all.x=TRUE)

#Attach the FPA Stream accumulations
edgedf <- merge(edgedf, fpa[,c('rid','RCA_FPA','ARCA_FPA')], 
                by='rid', all.x=TRUE)

#Derive pprca and pprsa areas in square meters
pp.zstats$PPRCA_SQM <- pp.zstats$PPRCA_COUNT * 900
pp.zstats$PPRSA_SQM <- pp.zstats$PPRSA_COUNT * 900
pp.zstats$PPRCA_NASQM <- pp.zstats$PPRCA_NACOUNT * 900
pp.zstats$PPRSA_NASQM <- pp.zstats$PPRSA_NACOUNT * 900

# remove all the NAs in the accumulated fields except fishpres, replace with zero
edgedf[!(names(edgedf) %in% "fishpres")][is.na(edgedf[!(names(edgedf) %in% 
                                                          "fishpres")])] <- 0

#In order to calculate pour point specific accumulations we need
#to remove the lowest rca from each existing accumulation
accumulated <- colnames(edgedf)[grep("^A",colnames(edgedf))]
rca <- gsub("^A","",colnames(edgedf)[grep("^A",colnames(edgedf))])

for (arca in accumulated) {
  rca.each <- gsub("^A","",arca)
  newcol <- paste("US_",arca,sep="")
  edgedf[,newcol] <- edgedf[,arca] - edgedf[,rca.each]
}

#Create a dataframe of just the US areas
edgedf.wo <- edgedf[,names(edgedf[,setdiff(names(edgedf),c(accumulated,rca))])]

#Bring together US areas with each station
obs.a <- merge(edgedf.wo, obs[,c('rid','STATION_KEY', 'LAT_RAW', 'LONG_RAW')], 
               by = 'rid')
#This one just brings in the 'SVN' column to allow for matching to the 
#disturbance areas which are to the sample date and not the station
obs.a <- merge(obs.a, phab.bugs[,c('STATION_KEY','SVN','Date')], 
               by = 'STATION_KEY', all.y = TRUE)
#first we need to trim the svn string column cause there are extra spaces in there
obs.a$SVN <- str_trim(obs.a$SVN)
#Bring in the year specific disturbances
obs.a <- merge(obs.a, dis, by = 'SVN', all.y = TRUE)
#resolve station key
obs.a <- within(obs.a, rm(STATION_KE))
#Bring in the clipped pour point areas for each variables
obs.a <- merge(obs.a, pp.zstats, by.x = 'STATION_KEY', by.y = 'STATION_KE', 
               all.x = TRUE)

# -----------------------------------------------------------
# Calculate Year specfic disturbance

disvar <- c("US_ADISRCA_1YR","US_ADISRCA_3YR","US_ADISRCA_10YR",
            "US_ADISRSA_1YR","US_ADISRSA_3YR","US_ADISRSA_10YR")

for (i in 1:nrow(obs.a)) {
  for (d in 1:length(disvar)) {
    if(!(is.na(obs.a[i,"Year_Sampled"]))) {
      y <- ifelse(obs.a[i,"Year_Sampled"] >2008, 2008, obs.a[i,"Year_Sampled"])
      var <- paste0(disvar[d],"_",as.character(y))
      obs.a[i,disvar[d]] <- ifelse(is.na(obs.a[i,var]),NA,obs.a[i,var])
    } else 
    {
      obs.a[i,disvar[d]] <- NA
    }
  }
}

rm(i,d,y,disvar)

#Now we can remove all the year specific disturbances
obs.a <- obs.a[, setdiff(names(obs.a), names(obs.a)[grep("[0-9]{4}$", 
                                                        names(obs.a))])]
#We can also remove variables we don't have complete data for lidar
obs.a <- obs.a[, setdiff(names(obs.a), names(obs.a)[grep("LI$", names(obs.a))])]
#We decided in the first round that we would only include high susceptibility
obs.a <- obs.a[, setdiff(names(obs.a), c('US_ASUSCEP1_DE',
                                        'US_ASUSCEP2_DE',
                                        'US_ASUSCEP3_DE'))]
#Align the disturbance column names
names(obs.a) <- gsub('DISTURB', 'DIS', names(obs.a))
names(obs.a)[grep('DIS', names(obs.a))] <- toupper(grep('DIS', 
                                                        names(obs.a), 
                                                        value = TRUE))

#Create vector of accumulated areas for looping
edge.accum <- grep("^US_",names(obs.a),value = TRUE)
edge.accum <- edge.accum[!edge.accum %in% c('US_AEDGELEN','US_ATYPEN')]

#Set NA values to 0 for the math to work below
obs.a[!(names(obs.a) %in% "fishpres")][is.na(obs.a[!(names(obs.a) %in% 
                                                       "fishpres")])] <- 0

#Loop through each variable at each accumulated scale and add the clipped pour
#point to the accumulated area that has had the lowest catchment area subtracted
for (ra in c('PPRCA', 'PPRSA', 'OTHER')){
  if (ra == 'OTHER') {
    edge.accum.sub <- edge.accum[!grepl('RCA', edge.accum) & 
                                   !grepl('RSA', edge.accum)]
  } else {
    edge.accum.sub <- edge.accum[grep(substr(ra, 3, 5), edge.accum)]
  }
  for (char in edge.accum.sub){
    if (char %in% grep('SQM|COUNT|FPA', edge.accum, value=TRUE)) {
      search <- gsub("US|ARCA|ARSA|RCA|RSA|LITH|_A|_RCA|^A|_", "", char)
    } else {
      search <- gsub("^(_A)", "", gsub("^US|ARCA|ARSA|RCA|RSA|LITH|
                                       (_A){0}|_RCA|^A|M$|(_DE)$", "", char))
    }
    if (search == 'SILT') {
      ppcol <- paste(ra, "_SILT", sep = "")
    } else if (ra == 'OTHER') {
      ppcol <- paste('PPRCA', search, sep="_")
      }
    else {
      ppcol <- grep(paste(ra, search, sep = '_'), names(obs.a), value=TRUE)  
    }
    if (ra == 'OTHER') {
      newcol <- paste("ARCA", search, sep='_')
    } else {
      newcol <- paste("A", gsub("PP", "", ra), "_", search, sep='')
    }
    obs.a[, newcol] <- obs.a[, char] + obs.a[, ppcol]
  }
}

#Now remove the US accumulations since we don't need them anymore
obs.a <- obs.a[, names(obs.a)[!names(obs.a) %in% grep("^US_", 
                                                      names(obs.a), 
                                                      value = TRUE)]]

# These are edgedf col that we want to delete from obs
edge.rm <- c("OBJECTID", "arcid", "from_node", 
             "to_node","HydroID", "GridID","NextDownID",
             "DrainID", "Shape_Length", "RCA_PI")

# clean up
obs.a <- obs.a[, !colnames(obs.a) %in% edge.rm]
names(obs.a) <- gsub("^PP","",names(obs.a))

#Switch the naming around so the names are easier to read with the variable name up front and the 
#scale at which is was calculated after (e.g. COMP_ARCA, COMP_ARSA, COMP_RCA, COMP_RSA)
names(obs.a)[grep("^ARCA|^ARSA|^RCA|^RSA", names(obs.a))] <- sapply(
  strsplit(names(obs.a)[grep("^ARCA|^ARSA|^RCA|^RSA", names(obs.a))], "A_"),
  function(x) {paste(x[2], "_", x[1], "A", sep = "")})


# -----------------------------------------------------------
#PRECIP

#Identify columns from precip to remove
rm.col <- c("EXCLUDE", "NewSample", "Bug_RefApril13", "TMDL",
            "MIDCOAST", "Samples", "Ref_Samples", "REF",
            "FLAG", "FLAG_REASON", "REF_shubler", "Bug_RefMay05",
            "SVN2KEY", "Date", "Month_Sampled", "Day_Sampled", 
            "Year_Sampled", "HabitatSampled", "SamplingAgency",
            "SamplingProtocol", "Field_QAQC",
            "Lab_QAQC", "Project_name", "TS_May05", "FSS_May05", 
            "PREDATOR_Nov05_model", "PREDATOR_Nov05_score",
            "PREDATOR_Nov05_Condition", "PREDATOR_outlier",
            "PREDATOR_Integrated_Rpt_2010", "Bug_Count_RIV",
            "PREDATOR_reference_model", "EcoT20_FURR",
            "EMAP_REF", "FURR_Score", "HDI_FURR", "SITE_NAME", 
            "Shannon_Original", "Shannon.First.Download")
precip <- precip[, !colnames(precip) %in% rm.col]

# These are samples to remove because of low bug counts. 
# They are only in the precip table.
svns.precip.rm <- c("00054CSR", "00105CSR", "01006CSR", "01014CSR",
                    "01020CSR", "01026CSR", "02125CSR", "03042CSR",
                    "04017CSR", "98086CSR", "99020CSR", "99051CSR",
                    "H108328", "H108330", "H108331", "H109844",
                    "H109845", "H109848", "H109875", "H109885",
                    "H109887", "H111024", "H113071", "H113082")

precip <- precip[!(precip$SVN %in% svns.precip.rm),]

#Need to clean up the SVN for those MIX samples 
precip$SVN <- str_trim(precip$SVN)

#We also need the annual averages
ppt <- within(ppt, rm(OBJECTID))

#Merge precip values with the rest of the obs
obs.a <- merge(obs.a, precip[, names(precip)[!names(precip) %in% 
                                               c('STATION_KEY', 
                                                 'days.f.origin')]], 
               by = 'SVN', all.x = TRUE)
obs.a <- merge(obs.a, ppt, by = 'STATION_KEY', all.x = TRUE)

rm(rm.col, svns.precip.rm)
# -----------------------------------------------------------
# Pull in the slope data, fix the col names, and append it together
ryan.slope <- rename(ryan.slope, c('SLOPE_AVG_MAP' = 'XSLOPE_MAP', 
                                   'Z_Min' = 'MIN_Z', 
                                   'Z_Max' = 'MAX_Z', 
                                   "REACHLEN" = 'RchLenFin'))

ryan.slope <- cbind(ryan.slope, 
                    data.frame("Source" = rep(NA, nrow(ryan.slope)), 
                               "CONFIDENCE" = rep(NA, nrow(ryan.slope)), 
                               "REACHLEN" = rep(NA, nrow(ryan.slope)),
                               "XSLOPE" = rep(NA, nrow(ryan.slope))))

ryan.slope <- within(ryan.slope, rm(OBJECTID, TYPE, Shape_Length))
ryan.slope$MIN_Z<- ryan.slope$MIN_Z * 0.3048
ryan.slope$MAX_Z <- ryan.slope$MAX_Z * 0.3048
slope <- rename(slope, c("SLOPE_AVG" = 'XSLOPE', 
                         "SLOPE_AVG_MAP" = 'XSLOPE_MAP'))
slope.all <- rbind(slope, ryan.slope)
slope.all <- slope.all[!duplicated(slope.all$STATION_KEY),]

slope.all$XSLOPE_MAP <- ((slope.all$MAX_Z - slope.all$MIN_Z) / 
                           slope.all$RchLenFin) * 100

# Station 13121 is missing from the slopes table but it's ok because 
# it's at the same location as 33298.
# Station 13121 and 33298 should be merged but c'est la vie
# it will cause problems with the other tables so 
# we just duplicate the slope values from 33298 and give it to 13121
X13121 <- slope.all[slope.all$STATION_KEY == "33298",]
X13121$STATION_KEY <- "13121"
slope.all <- rbind(slope.all, X13121)

#Merge with the rest of the obs
obs.a <- merge(obs.a, slope.all[,c('STATION_KEY','XSLOPE_MAP','MIN_Z')], 
               by = 'STATION_KEY', all.x = TRUE)

#Fix the negative stream slope values
obs.a[obs.a$XSLOPE_MAP < 0,'XSLOPE_MAP'] <- 0.001

rm(slope, ryan.slope, X13121)

# -----------------------------------------------------------
# Pull in the measured habitat data

# This would be interesting to explore further but since we only have this phab data at a handful of sites
# we are not incorpoating it at this time

#Not sure we need the phab.bugs since we don't have data at all the sites
#comb <- merge(phab.bugs, obs.a, by = 'SVN', all = TRUE)

# -----------------------------------------------------------
# Calculate Stream power 
# fist pull in the nhd comiD
nhd <- rename(nhd, c('FEATUREID' = 'NHDP21_COMID'))
nhd <- nhd[, c('SVN', 'NHDP21_COMID')]
obs.a <- merge(obs.a, nhd, by = 'SVN', all.x = TRUE)

# Calculate stream power by first pulling in the flow data from nhdplus v21
nhd.flow <- nhd.flow[, c('Comid', 'Q0001E')]

#Some corrections to incorrect NHD ComID mapping
nhd.correction <- c(23872143, 23876155, 23881226, 23915183)
names(nhd.correction) <- c(34617, 22504, 34673, 34682)
obs.a[obs.a$STATION_KEY %in% names(nhd.correction), 'NHDP21_COMID'] <- obs.a[
  obs.a$STATION_KEY %in% names(nhd.correction), 'STATION_KEY']
obs.a$NHDP21_COMID <- mapvalues(obs.a$NHDP21_COMID, 
                                from = names(nhd.correction), 
                                to = nhd.correction)

#fill in based on similar area, rainfall and location for stations that fall on tribs within large catchments
obs.a[obs.a$STATION_KEY %in% c('23817', '33323', '33333', '33355', '35786', 
                               'dfw_39723', 'dfw_49477', 'dfw_795'), 
      'NHDP21_COMID'] <- c(23876327, 23890092, 23890512, 23882258, 23876785, 
                           23882086, 23920738, 23876439)
obs.a <- merge(obs.a, nhd.flow, by.x = "NHDP21_COMID", 
               by.y = 'Comid', all.x = TRUE)

#Bring in NHD cumulative area
obs.a <- merge(obs.a, nhd.area, by.x = 'NHDP21_COMID', 
               by.y = 'ComID', all.x = TRUE)

#convert SQM_RCA to square kilometers
obs.a$SQKM_ARCA <- obs.a$SQM_ARCA / 1e6

#This one had the catchment drawn wrong so this is just forcing it to the right value
#comb[comb$STATION_KEY == 34617,'ppSQkm'] <- comb[comb$STATION_KEY == 34617,'TotDASqKM']

#Then we need to derive the scaling ratio based on watershed area (instead of using flow at the base of the NHD and because we don't have/need M values)
obs.a$nhd_ratio <- 1 - ((obs.a$TotDASqKM - obs.a$SQKM_ARCA) / (obs.a$TotDASqKM))

#Use the ratio to scale the flow estimate
obs.a$Q0001E_adj <- obs.a$Q0001E * obs.a$nhd_ratio

#Calculate stream power
obs.a$STRMPWR <- obs.a$XSLOPE_MAP * obs.a$Q0001E_adj

rm(nhd.flow, nhd)
# -----------------------------------------------------------
#Pull in updated FSS values and remove the SVNs that have low counts
#Post revision we are using the output from the list of stations already filtered for low count
obs.a <- merge(obs.a, fss[, c('SVN', 'FSS_26Aug14')], by = 'SVN', all.x = TRUE)

rm(fss)
# -----------------------------------------------------------
# Clean up date columns and formatting

obs.a$DATE <- as.POSIXlt(obs.a$Date, format="%m/%d/%Y")
obs.a$YEAR <-  obs.a$DATE$year + 1900
obs.a$MONTH <- obs.a$DATE$mon + 1 # +1 because it is zero-indexed
obs.a <- within(obs.a, rm(datestr, Q0001E, Year_Sampled))

# -----------------------------------------------------------
# Run the loop to calculate percentages, means, or densities

pvar2 <- pvar[pvar$numerator %in% names(obs.a),]

for (i in 1:nrow(pvar2)) {
  obs.a[, pvar2[i, 1]] <- (obs.a[, pvar2[i, 2]] / 
                             obs.a[, pvar2[i, 3]]) * pvar2[i, 4]
  #comb <- propor(comb, pvar[i,3], pvar[i,4], pvar[i,2])
}

# -----------------------------------------------------------

#Clean up the columns in obs.a a bit
obs.a <- within(obs.a, rm(fishpres, nhd_ratio, TotDASqKM, DivDASqKM, Date, 
                          COUNT_ARCA, COUNT_ARSA, COUNT_RCA, COUNT_RSA,
                          NACOUNT_ARCA, NACOUNT_RCA, NACOUNT_RSA, NASQM_ARCA, 
                          NASQM_RCA, NASQM_RSA, NHDP21_COMID))

#The previous workflow identified columns with NA data and based on the 
# number of NAs in the variable the determination for inclusion/exclusion was 
#made. The variables in these current obs are the result of those determinations.

#There is still a need to clean up the dataframe to just the variables we are relating
#Here we will output a csv and mark which columns to keep and which to exclude for subsequent steps
# Leave commented unless you changed something above because this file has been edited outside of the script.
# old <- read.csv('VarNames_RF_v2.csv')
# new <- data.frame(var = names(obs.a))
# ready <- merge(new, old, by = 'var', all.x = TRUE)
# x <- 1:length(names(obs.a))
# names(obs.names) <- x
# obs.names <- names(obs.a)
# ready$col.num <- mapvalues(ready$var, from = obs.names, to = names(obs.names))
# ready[ready$var %in% c('LAT_RAW', 'LONG_RAW'),'keep'] <- 1
# write.csv(names(obs.a), 'VarNames_RF_v2.csv')


# -----------------------------------------------------------
# Write the files

#write.csv(obs.a, 'ssn_RF_data.csv', row.names = FALSE)



