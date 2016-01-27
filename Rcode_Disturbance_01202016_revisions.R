library(RODBC)

con2 <- odbcConnectAccess('C:/users/pbryant/desktop/midcoasttmdl-gis/lsn05_watersheds/lsn04/Tables.mdb')
to_append <- sqlFetch(con2, 'ssn_edges_table_final')
to_append <- to_append[,-grep("DIS", names(to_append))]

edges <- read.csv('C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/LSN05_Watersheds/LSN04/ADIS_revised_01202016.txt', stringsAsFactors = FALSE)

dis <- edges[,c('rid',sort(grep('DIS',names(edges), value = TRUE)))]
dis <- as.data.frame(lapply(dis, function (x) as.numeric(gsub(",",'',x))))

to_save <- merge(to_append, dis, by = 'rid')

sqlSave(con2, to_save, tablename = "ssn_edges_table_final_dis_update")

conT <- odbcConnectAccess('T:/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN04/Tables.mdb')

sqlDrop(conT, sqtable = 'ssn_edges_table_final')
sqlSave(channel = conT, 
        dat = to_save, 
        tablename = 'ssn_edges_table_final', 
        rownames = FALSE,)

odbcCloseAll()
