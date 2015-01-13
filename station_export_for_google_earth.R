export <- obs.vars[,c('SVN','STATION_KEY','SITE_NAME','DATE',"FSS_26Aug14",
                 "sum_1095_days", 
                 "PDISRSA_1YR",
                 "PALITHERODRCA", "PASILTRCA", "PACLAYRCA",
                 "DAPOPRCA2010","POWNRCA_PRI","PAOWNRCA_AGR","DAROADX",
                 "PASUSCEP5_DE",
                 "STRMPWR", "XSLOPE_MAP","MIN_Z","LAT_RAW")]

impaired <- data.frame(STATION_KEY = c(21842,34660,21792,33361,26818,33418,33417,34695,26822,33320,33333,30403,34665,26816,25297,26964,29906,33327),
                       TMDL_Target = c(14,14,14,3,7,rep(14,6),8,14,8,14,8,14,14))

export <- merge(export, impaired, by = 'STATION_KEY', all.x = TRUE)

export$Impaired <- ifelse(is.na(export$TMDL_Target),'No','Yes')

write.csv(export, "E:/Github/KML-R-tool/stns_for_ge.csv")
