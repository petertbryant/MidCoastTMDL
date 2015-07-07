#Compile pour point watershed PPARCA and PPARSA values
library(RODBC)
library(plyr)

options(stringsAsFactors = FALSE)

con <- odbcConnectAccess("C:/users/pbryant/desktop/midcoasttmdl-gis/pprca_zonalstats.mdb")
con2 <- odbcConnectAccess("//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/ssn/lsn04/tables.mdb")

tablenames.all <- sqlTables(con, tableType = "TABLE")
tablenames <- tablenames.all[grep("^PP",tablenames.all$TABLE_NAME),'TABLE_NAME']

rid.change <- c('VALUE_' = 'rid_LSN04')

for (char in c('CLAY','COMP','EROD','KFACTWS','pop','SAND','SILT','SILT_CLAY','SUSCEP','OWNCLASS', 'BLM')) {
  tbls.by.char <- tablenames[grep(char,tablenames)]
  if (char == 'CLAY') {
    tbls.by.char <- tbls.by.char[!grepl('SILT',tbls.by.char)]
  }
  if (char == 'SILT') {
    tbls.by.char <- tbls.by.char[!grepl('CLAY',tbls.by.char)]
  }
  for (ra in c('PPRCA','PPRSA')) {
    tbls <- tbls.by.char[grep(ra, tbls.by.char)]
    for (i in 1:length(tbls)) {
      tmp <- sqlFetch(con, tbls[i])
      if (nrow(tmp) > 0) {
        if (char == 'SILT_CLAY') {
          newcol.name1 <- paste(strsplit(tbls[i],"_")[[1]][c(1,4,5)],collapse="_")
          newcol.name2 <- paste(ra,"NACOUNT",sep="_")
          tmp <- rename(tmp, c('MEAN' = newcol.name1, "COUNT_" = newcol.name2), warn_missing = FALSE)
          tmp <- tmp[,c('STATION_KE',newcol.name1,newcol.name2)]
        } else if (char == 'SUSCEP') {
          newcol.name1 <- paste(ra, "_SUSCEP4",sep="")
          newcol.name2 <- paste(ra, "_SUSCEP5",sep="")
          tmp <- rename(tmp, c('VALUE_4' = newcol.name1, 
                               'VALUE_5' = newcol.name2), warn_missing = FALSE)
          if (newcol.name1 %in% names(tmp) & newcol.name2 %in% names(tmp)) {
            tmp <- tmp[,c('STATION_KE',newcol.name1,newcol.name2)]
          } else if (newcol.name1 %in% names(tmp) & !newcol.name2 %in% names(tmp)) {
            tmp <- tmp[,c('STATION_KE',newcol.name1)]
            tmp[,newcol.name2] <- NA
          } else {
            next()
          }
        } else if (char == 'OWNCLASS') {
          odf <- paste(ra, "_OWN_ODF", sep = "")
          fed <- paste(ra, "_OWN_FED", sep = "")
          urb <- paste(ra, "_OWN_URB", sep = "")
          agr <- paste(ra, "_OWN_AGR", sep = "")
          pri <- paste(ra, "_OWN_PRI", sep = "")
          tmp <- rename(tmp, c('VALUE_2' = odf, 
                               'VALUE_8' = fed, 
                               "VALUE_9" = urb,
                               "VALUE_10" = agr), warn_missing = FALSE)
          tmp[,pri] <- tmp$VALUE_1 + tmp$VALUE_3
          if (!odf %in% names(tmp)) {
            tmp[,odf] <- 0
          } 
          if (!fed %in% names(tmp)) {
            tmp[,fed] <- 0
          } 
          if (!urb %in% names(tmp)) {
            tmp[,urb] <- 0 
          } 
          if (!agr %in% names(tmp)) {
            tmp[,agr] <- 0
          }
          tmp <- tmp[,c("STATION_KE",pri,odf,fed,urb,agr)]
        } else if (char == 'pop') {
          count.col <- paste(ra, "COUNT", sep = "_")
          newcol.name <- paste(strsplit(tbls[i],"_")[[1]][c(1,4)],collapse="_")
          tmp <- rename(tmp, c("SUM_" = newcol.name,
                               "COUNT_" = count.col), warn_missing = FALSE)
          tmp <- tmp[,c('STATION_KE',newcol.name,count.col)]
        } else {
          newcol.name <- paste(strsplit(tbls[i],"_")[[1]][c(1,4)],collapse="_")
          tmp <- rename(tmp, c('VALUE_1' = newcol.name, 
                               'MEAN' = newcol.name, 
                               "SUM_" = newcol.name,
                               "LENGTH" = newcol.name), warn_missing = FALSE)
          tmp <- tmp[,c('STATION_KE',newcol.name)]
        }
        
        ifelse(i == 1, tblFull <- tmp, tblFull <- rbind(tblFull, tmp))
      }
    }
    ifelse(ra == 'PPRCA', tblFull2 <- tblFull, tblFull2 <- merge(tblFull2, tblFull, by = 'STATION_KE', all = TRUE))
  } 
  ifelse(char == 'CLAY', wcdf <- tblFull2, wcdf <- merge(wcdf, tblFull2, by = 'STATION_KE', all = TRUE))
}

wcdf <- rename(wcdf, c('PPRCA_s1' = 'PPRCA_ROADLEN',
                     'PPRSA_s1' = 'PPRSA_ROADLEN'))

roadx <- sqlFetch(con, 'PPRCA_ROADX')
roadx <- roadx[,c('STATION_KE','PNT_COUNT')]
roadx <- rename(roadx, c('PNT_COUNT' = 'PPRCA_ROADX'))
wcdf <- merge(wcdf, roadx, by = 'STATION_KE', all = TRUE)

splash <- sqlFetch(con, 'PPRCA_SPLASH')
splash <- splash[,c('STATION_KE','PNT_COUNT')]
splash <- rename(splash, c('PNT_COUNT' = 'PPRCA_SPLASH'))
wcdf <- merge(wcdf, splash, by = 'STATION_KE', all = TRUE)

fish <- sqlFetch(con, 'PPRCA_TYPEF')
typef <- fish[fish$Fishpres == 'Fish',c('STATION_KE','LENGTH')]
typef <- rename(typef, c('LENGTH' = 'PPRCA_TYPEF'))
wcdf <- merge(wcdf, typef, by = 'STATION_KE', all = TRUE)

fish <- ddply(fish, .(STATION_KE), summarize, PPRCA_FPA = sum(LENGTH))
wcdf <- merge(wcdf, fish[,c("STATION_KE","PPRCA_FPA")], by = 'STATION_KE', all = TRUE)

con2.tbls <- sqlTables(con2, tableType = "TABLE")

if ("PPRCA_PPRSA_ZonalStats" %in% con2.tbls$TABLE_NAME) {
  sqlDrop(con2, "PPRCA_PPRSA_ZonalStats")
}

sqlSave(con2, wcdf, tablename = 'PPRCA_PPRSA_ZonalStats', varTypes = c('STATION_KE' = 'VARCHAR(255)'), rownames = FALSE)
