library(RODBC)
library(plyr)
library(foreign)
library(stringr)

options(stringsAsFactors = FALSE)

#Establish data connections
con <- odbcConnectAccess("C:/users/pbryant/desktop/midcoasttmdl-gis/PPRCA_PPRSA_Disturbance.mdb")
con2 <- odbcConnectAccess("//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/ssn/lsn04/tables.mdb")

#brind in data
svn <- read.csv("//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Station_Selection/Watershed_Char_phab_bugs_merge_SSN_FINAL.csv")
svn$SVN <- str_trim(svn$SVN)
obs <- read.dbf("//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/ssn/lsn05/lsn.ssn/sites.dbf")
rid.cross <- sqlFetch(con2, "LSN0405_Crosswalk")
tables <- sqlTables(con)
tblnames <- tables$TABLE_NAME[grep("DISTURB",tables$TABLE_NAME)]

#Combine obs with svn to get year sampled associated with each rid
ry <- merge(obs, svn[,c('SVN','Year_Sampled')], by = 'SVN', all.x = TRUE)
ry04 <- merge(ry, rid.cross, by.x = 'rid', by.y = 'rid_LSN05', all.x = TRUE)

#Put all the individual tabulations back into a single table first and then put the complete
#tables together into a single dataframe
for (yr in 1996:2008) {
  tbls.by.yr <- tblnames[grep(yr,tblnames)]
  for (tp in c('10yr', '1yr', '3yr')) {
    tbls.by.yr.tp <- tbls.by.yr[grep(tp,tbls.by.yr)]
    for (ra in c('PPRCA', 'PPRSA')) {
      tbls <- tbls.by.yr.tp[grep(ra, tbls.by.yr.tp)]
      for (i in 1:length(tbls)) {
        tmp <- sqlFetch(con, tbls[i])
        if (nrow(tmp) > 0) {
          newcol.name <- paste(strsplit(tbls[i],"_")[[1]][c(1,4,5,6)],collapse="_")
          tmp <- rename(tmp, c('VALUE_1' = newcol.name))
          tmp <- within(tmp, rm(OBJECTID))
          ifelse(i == 1, tblFull <- tmp, tblFull <- rbind(tblFull, tmp))
        }
      }
      ifelse(ra == 'PPRCA', tblFull2 <- tblFull, tblFull2 <- merge(tblFull2, tblFull, by = 'STATION_KE', all = TRUE))
    }
    ifelse(tp == '10yr', int <- tblFull2, int <- merge(int, tblFull2, by = 'STATION_KE', all = TRUE))
  }
  ifelse(yr == 1996, ddf <- int, ddf <- merge(ddf, int, by = 'STATION_KE', all = TRUE))
}

#Put all the disturbance values together into a single dataframe
# for (i in 1:length(tblnames)) {
#       tmp <- sqlFetch(con, tblnames[i])
#       newcol.name <- paste(strsplit(tblnames[i],"_")[[1]][c(1,4,5,6)],collapse="_")
#       tmp <- rename(tmp, c('VALUE_1' = newcol.name))
#       tmp <- within(tmp, rm(OBJECTID))
#       
#       ifelse(i == 1, ddf <- tmp, ddf <- merge(ddf, tmp, by = 'STATION_KE', all = TRUE))
# }

#Combine the obs with year sampled and the disturbances
#ry <- data.frame(rid_LSN04 = sample(ddf$rid_LSN04, 780, replace = TRUE), Year_Sampled = rep(1996:2008,783/13)) #for testing
ryd <- merge(ry04, ddf, by = "STATION_KE", all.x = TRUE)

newcols <- c()
for (scale in c('PPRCA','PPRSA')) {
  for (yr in c('1yr', '3yr', '10yr')) {
    newcol <- paste(scale, "_DISTURB_", yr, sep = "")
    newcols <- c(newcols,newcol)
    for (i in 1:nrow(ryd)) {
      lu_year <- ifelse(ryd$Year_Sampled[i] > 2008,2008,ryd$Year_Sampled[i])
      ryd[i,newcol] <- ryd[i,grep(paste(newcol, "_", lu_year, sep = ""),names(ryd),value = TRUE)]
    }
  }
}

out <- ryd[c('SVN','STATION_KE','rid_LSN04','Year_Sampled',newcols)]

sqlSave(con2, out, tablename = 'PPRCA_PPRSA_Disturbance', rownames = FALSE, varTypes = c('SVN' = 'VARCHAR(255)',
                                                                                         'STATION_KE' = 'VARCHAR(255)',
                                                                                         'rid_LSN04' = 'INTEGER',
                                                                                         'Year_Sampled' = 'INTEGER'))
