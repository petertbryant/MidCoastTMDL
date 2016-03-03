sam_df <- impaired
precip1 <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/watershed_characteristics/precip/r_output_ws_dfl_dmean_precip_daily_mean_2014-09-08.csv')

# This is the interval of days prior to the sample we make calculations over
sum_it <- c(1095)
sum_it_col <- paste0("sum_",sum_it,"_days")

for(i in 1:nrow(sam_df)) {
  station <- sam_df$STATION_KEY[i]
  dfs <- precip1[precip1$STATION_KEY == station,]
  for(s in 1:nrow(dfs)) {
    day_start <- s
    day_end <- s + 1095
    dfr <- dfs[dfs$days.f.origin %in% c(day_start:day_end),]
    newrow <- data.frame(STATION_KEY = station, datestr = dfr[1096,'datestr'], sum_1095_days = sum(dfr$value))
    ifelse(s == 1, dfd <- newrow, dfd <- rbind(dfd, newrow))
  }
  ifelse(i == 1, dfdall <- dfd, dfdall <- rbind(dfdall, dfd))
}

dfdall <- dfdall[!is.na(dfdall$datestr),]
dfdall$datestr <- factor(as.character(dfdall$datestr))

save(dfdall, file = 'C:/users/pbryant/desktop/midcoasttmdl-gis/precip_daily_sum_1095_days.Rdata')

#plot of all sites
boxplot(dfdall$sum_1095_days~dfdall$STATION_KEY,las=2)
points(obs.complete[obs.complete$STATION_KEY %in% unique(dfdall$STATION_KEY),'sum_1095_days'] ~ 
         factor(as.character(obs.complete[obs.complete$STATION_KEY %in% unique(dfdall$STATION_KEY),'STATION_KEY'])),
       pch=19,cex=1.2,col='red')

#plot for each site alone
for (i in 1:length(unique(dfdall$STATION_KEY))) {
  main_text = paste('sum_1095_days at',unique(dfdall$STATION_KEY)[i])
  boxplot(dfdall[which(dfdall$STATION_KEY == unique(dfdall$STATION_KEY)[i]),'sum_1095_days'],main=main_text,ylim=c(3000,14000))
  ind.obs <- obs.complete[obs.complete$STATION_KEY == unique(dfdall$STATION_KEY)[i],'sum_1095_days']
  for (j in 1:length(ind.obs)) {
    points(ind.obs[j],pch=19,cex=1.2,col='red')
  }
}




