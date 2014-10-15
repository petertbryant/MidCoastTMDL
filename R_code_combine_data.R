library(plyr)
library(RODBC)
library(SSN)
options(stringsAsFactors = FALSE)

## READ DATA FROM ACCESS 2007
indb <- "//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/Tables.mdb"
tablename <- "edge_table_final"
channel <-odbcConnectAccess2007(indb)
edgedf <- sqlFetch(channel, tablename)
close(channel)
rm(indb, tablename, channel)

#read in the SSN object
bugs <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/lsn.ssn", o.write=FALSE)
obs<- getSSNdata.frame(bugs, Name = "Obs")

#put edge accumulations with the station
obs.a <- merge(obs, edgedf, by = 'rid', all.x = TRUE)
obs.a <- rename(obs.a, c('STATION_KE' = 'STATION_KEY'))

precip <- read.csv('//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Precip/R_output_Precip_samples_2014-09-08.csv')
phab.bugs <- read.csv('//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Station_Selection/Watershed_Char_phab_bugs_merge_SSN_FINAL.csv')
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

precip2 <- precip[,names(precip)[!names(precip) %in% names(phab.bugs)]]
precip <- cbind(data.frame('SVN' = precip[,'SVN']),precip2)
rm(precip2)

comb <- merge(phab.bugs, precip, by = 'SVN', all = TRUE)

comb2 <- comb[,names(comb)[!names(comb) %in% names(slope.all)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']),comb2)
rm(comb2)

comb <- merge(comb, slope.all, by = 'STATION_KEY', all = TRUE)

comb2 <- comb[,names(comb)[!names(comb) %in% names(obs.a)]]
comb <- cbind(data.frame('STATION_KEY' = comb[,'STATION_KEY']),comb2)
rm(comb2)

comb <- merge(comb, obs.a, by = 'STATION_KEY', all.y = TRUE)

View(comb[,1:99])
View(comb[,100:199])
View(comb[,200:299])
View(comb[,300:399])
View(comb[,c(1:2,400:496)])
View(comb[,500:505])

na.list <- apply((apply(comb, 2, is.na)),2,table)
na.df <- data.frame('col' = names(na.list))
na.df$na.count <- NA
for (i in 1:length(na.list)) {
  if(!is.na((na.list[i][[1]][2]))) {
    na.df$na.count[i] <- na.list[i][[1]][2]
  }
}