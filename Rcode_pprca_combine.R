#Compile pour point watershed PPARCA and PPARSA values
library(RODBC)
library(plyr)

options(stringsAsFactors = FALSE)

con <- odbcConnectAccess("C:/users/pbryant/desktop/midcoasttmdl-gis/pprca_zonalstats.mdb")

tablenames <- sqlTables(con, tableType = "TABLE")
tablenames <- tablenames[grep("^PP",tablenames$TABLE_NAME),'TABLE_NAME']

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
          newcol.name <- paste(strsplit(tbls[i],"_")[[1]][c(1,4,5)],collapse="_")
          tmp <- rename(tmp, c('VALUE_1' = newcol.name, 'MEAN' = newcol.name, "SUM_" = newcol.name), warn_missing = FALSE)
          tmp <- tmp[,c('STATION_KE',newcol.name)]
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
        } 
        else {
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
  ifelse(char == 'CLAY', int <- tblFull2, int <- merge(int, tblFull2, by = 'STATION_KE', all = TRUE))
}

int <- rename(int, c('PPRCA_s1' = 'PPRCA_ROADLEN',
                     'PPRSA_s1' = 'PPRSA_ROADLEN'))

roadx <- sqlFetch(con, 'PPRCA_ROADX')
roadx <- roadx[,c('STATION_KE','PNT_COUNT')]
roadx <- rename(roadx, c('PNT_COUNT' = 'PPRCA_ROADX'))
int <- merge(int, roadx, by = 'STATION_KE', all = TRUE)

splash <- sqlFetch(con, 'PPRCA_SPLASH')
splash <- splash[,c('STATION_KE','PNT_COUNT')]
splash <- rename(splash, c('PNT_COUNT' = 'PPRCA_SPLASH'))
int <- merge(int, splash, by = 'STATION_KE', all = TRUE)

fish <- sqlFetch(con, 'PPRCA_TYPEF')
fish <- fish[fish$Fishpres == 'Fish',c('STATION_KE','LENGTH')]
fish <- rename(fish, c('LENGTH' = 'PPRCA_TYPEF'))
int <- merge(int, fish, by = 'STATION_KE', all = TRUE)

rm(list = setdiff(ls(),'int'))
