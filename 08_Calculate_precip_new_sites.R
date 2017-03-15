# Read Daily PRISIM precip rasters
# read the 24 hour freq rasters
# read the watershed shapefile
# read the sample data

#iterate through each watershed
## output matrix/vector of daily precip averaged for each watershed
## output matrix/vector of the average 24 hour freq interval ppt for each watershed

#iterate through each sample
## grab the date/year of sample

## sum of average daily upstream precip for the following time periods preceeding the sample date
### 60 days  (2 months)
### 180 days (6 months)
### 365 days (1 year)
### 1095 days(3 years)

## Calculate the number of days that exceed
### 24 hour 6 month event
### 24 hour 2 year event
### 24 hour 10 year event
### 24 hour 25 year event
### 24 hour 50 year event
### 24 hour 100 year event
### for the following time periods preceeding the sample date
### 60 days  (2 months)
### 180 days (6 months)
### 365 days (1 year)
### 1095 days(3 years)

## calculate how many days have passed since the following precipitaiton event threshold using the average daily upstream precip:
### 24 hour 6 month event
### 24 hour 2 year event
### 24 hour 10 year event
### 24 hour 25 year event
### 24 hour 50 year event
### 24 hour 100 year event

# output the results 

library(rgdal)
library(raster)
library(reshape)
thedate <-Sys.Date()

output_path <- "//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Precip/"
path_24ppt <- "D:/ArcGIS/Layers/climate/24_hour_precip/rasters/img/"
path_dppt <- "D:/ArcGIS/Layers/PRISM/climate/precipitation/Daily/"
path_ws <- "C:/WorkSpace/Biocriteria/WatershedCharaterization/SSN/Stations/watersheds/watersheds_ssn_GCS_NAD83.shp"
path_sta <- "C:/WorkSpace/Biocriteria/WatershedCharaterization/SSN/Stations/Stations_SSN_20140822_GCS_NAD83.shp"
path_samples <-"//Deqhq1/tmdl/TMDL_WR/MidCoast/Models/Sediment/Watershed_Characteristics/Station_Selection/Watershed_Char_Bug_Samples_ADD_SSN_20140822_FINAL.csv"

ffiles <- c("6month.img","2year.img","10year.img","25year.img","50year.img","100year.img")
rfiles <- list.files(path_dppt, pattern="\\.bil$")
rnames <- substr(rfiles,24,31)

# This is the projection of the various raster data I am working with
proj_harnlambert <-CRS("+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999999999 +y_0=0 +datum=NAD83 +units=ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
proj_GCSNAD83 <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

########################################################################################################
# import watershed shapefile
ws <- shapefile(path_ws,stringsAsFactors=FALSE)

#convert to a data frame and fix the station col name
ws_df <- data.frame(ws)
colnames(ws_df)[2] <- "STATION_KEY"

#import station shapefile
sta <- shapefile(path_sta,stringsAsFactors=FALSE)

#convert to a data frame and fix the station col name
sta_df <- data.frame(sta)
colnames(sta_df)[2] <- "STATION_KEY"
colnames(sta_df)[20] <- "LON_FACv21"

# get the x/y for each of the stations into the correct format for extract().
sta_xy <- cbind(sta_df$LON_FACv21,sta_df$LAT_FACv21)
sta_xy <- SpatialPoints(sta_xy, proj4string=CRS(proj4string(sta)), bbox = NULL)

# import samples table
sam_df <- read.table(path_samples, header=TRUE, sep=",", na.strings= c("NA") ,check.names=FALSE, stringsAsFactors=FALSE, fill=TRUE)


########################################################################################################
########################################################################################################
# import PRISIM daily precip rasters (data is in mm) and
# calculate the means for each freq interval for each watershed

# Get the 1st daily precip raster and reproject watershed shapefile into that projection if they are different
r1 <- raster(paste0(path_dppt,rfiles[1]))
if(proj4string(ws) != proj4string(r1)) {
  ws <- spTransform(ws,CRS(proj4string(r1))) }

# number of iterations (THIS TALKES AWHILE)
rit <- 1:length(rfiles)

for(i in rit) { 
  r1 <- raster(paste0(path_dppt,rfiles[i]))
  r1 <- crop(r1, extent(ws))
  # create stack for the daily precip
  if(i == 1) {
    s1 <-stack(r1)
  } else {
    s1 <- addLayer(s1,r1)
  }
}


names(s1) <- rnames[rit] #add the date as layer name

#########################################
save.image(file="precip_data.RData")
#load("precip_data.RData")
#########################################

# Since some of the watersheds are smaller than a the precip raster cell size the extraction will fail. To fix this I'm 
# disaggregating the raster into a smaller cell size by a factor of 5 (4KM -> 0.8km).
# The cell values are the same as the larger raster since I am not using bilinear interpolation
# THIS TAKES A LONG TIME
s_ws <- disaggregate(s1, fact=5)

# Calculate the mean daily precip for each watershed and output as a dataframe 
# THIS TAKES A LONG TIME
ws_dmean <- extract(s_ws, ws, weights=FALSE, small= TRUE, fun=mean, df=TRUE)

#########################################
save.image(file="precip_data.RData")
#load("precip_data.RData")
#########################################

########################################################################################################
## Just some general formatting

# Add the station ID to the data frame
ws_dfw_24mean <- cbind(ws_df[2],ws_24mean)
ws_dfw_dmean <- cbind(ws_df[2],ws_dmean)

# Drop the ID column name for the melt
ws_dfw_24mean$ID <- NULL
ws_dfw_dmean$ID <- NULL

# put everthing in long format
ws_dfl_24mean <- melt(ws_dfw_24mean, id.vars=c("STATION_KEY"))
ws_dfl_dmean <- melt(ws_dfw_dmean, id.vars=c("STATION_KEY"))

# remove the X in front of all the dates
ws_dfl_dmean$datestr <- gsub("[X]","",ws_dfl_dmean$variable)
ws_dfl_dmean$datestr <- as.POSIXlt(as.character(ws_dfl_dmean$datestr),format="%Y%m%d")

# Get the sample date into POSIX format
sam_df$datestr <- strptime(sam_df$Date,format="%m/%d/%Y")

# This is start date of the daily prcip data
origin_sec <- as.POSIXlt(as.character("19950101"),format="%Y%m%d")

#calculate the number of days since 1/1/1995 (start of daily precip data)
sam_df$days.f.origin <- as.integer(((as.numeric(sam_df$datestr) - as.numeric(origin_sec)) / 86400) +.5)
ws_dfl_dmean$days.f.origin <- as.integer(((as.numeric(ws_dfl_dmean$datestr) - as.numeric(origin_sec)) / 86400) +.5)

########################################################################################################
## sum of average daily  precip for the following time periods preceeding the sample date
### 60 days  (2 months)
### 180 days (6 months)
### 365 days (1 year)
### 1095 days(3 years)
## Calculate the number of days that exceed the
### 24 hour 6 month event
### 24 hour 2 year event
### 24 hour 10 year event
### 24 hour 25 year event
### 24 hour 50 year event
### 24 hour 100 year event
### for the following time periods preceeding the sample date
### 60 days  (2 months)
### 180 days (6 months)
### 365 days (1 year)
### 1095 days(3 years)


# this sets up the column names and the various iteration counts for the for routines

# This is the interval of days prior to the sample we make calculations over
sum_it <- c(60,180,365,1095)
sum_it_col <- paste0("sum_",sum_it,"_days")
# Name of the columns from ws_24mean that hold the 24 hour metrics
count_it <- c("X6month","X2year", "X10year", "X25year", "X50year", "X100year")
# Build the count column names
count_it2 <- paste0("count_",sum_it,"_days")
count_it_col <-as.vector(outer(count_it, count_it2, paste, sep="_"))

# Add the columns to the sample dataframe
sam_df[,sum_it_col] <- NA
sam_df[,count_it_col] <- NA

#i <- 1 #TEST
#day_end <- 355 # TEST
#s <- 2 # TEST
#c <-1 # TEST

for(i in 1:nrow(sam_df)) {
  station <- sam_df$STATION_KEY[i]
  day_end <- sam_df$days.f.origin[i]
  x <- 1
  for(s in 1:length(sum_it)) {
    day_start <- day_end - sum_it[s]
    df <- subset(ws_dfl_dmean,days.f.origin %in% c(day_start:day_end) & STATION_KEY == station)
    sam_df[,sum_it_col[s]][i] <- sum(df$value)
    
    for(c in 1:length(count_it)) {
      # Get the average daily upstream precip metric for the sample we are on
      metric <- ws_dfw_24mean[ws_dfw_24mean$STATION_KEY %in% station, count_it[c]]
      
      # sum the number of days that exceed that metric and add it back into the sample data frame under that metric column.
      sam_df[,count_it_col[x]][i] <- sum(df$value > metric)
      x <- x +1
    }
  }
}

########################################################################################################
## calculate how many days have passed since the following precipitaiton event threshold using the average daily upstream precip:
### 24 hour 6 month event
### 24 hour 2 year event
### 24 hour 10 year event
### 24 hour 25 year event
### 24 hour 50 year event
### 24 hour 100 year event


# Build the metric count column names and add it to the sample dataframe
days_since_col <- paste0("days_since_",count_it)
sam_df[,days_since_col] <- NA

for(i in 1:nrow(sam_df)) {
  station <- sam_df$STATION_KEY[i]
  day_end <- sam_df$days.f.origin[i]
  df2 <- subset(ws_dfl_dmean,days.f.origin %in% c(0:day_end) & STATION_KEY == station)
  df2$days.f.sam <- df2$days.f.origin - max(df2$days.f.origin)
  
  for(c in 1:length(count_it)) {
    metric <- ws_dfw_24mean[ws_dfw_24mean$STATION_KEY %in% station, count_it[c]]
    df3 <- subset(df2, value >= metric)
    
    # sum the number of days that equal or exceed that metric and add it back into the sample data frame under that metric column. 
    # INF means it never happened before 1995
    sam_df[,days_since_col[c]][i] <- max(df3$days.f.sam) * -1
    
  }
}


rm(df,df2,df3,c,i,x)


#########################################
save.image(file="precip_data.RData")
#load("precip_data.RData")
#########################################

########################################################################################################
## Add percent of 24 hour metric as a daily timeseries. This will go into the long format dataframe

df1 <- ws_dfl_dmean
df1$PX6month <- NA
df1$PX2year <- NA
df1$PX10year <- NA
df1$PX25year <- NA
df1$PX50year <- NA
df1$PX100year <- NA

for(i in 1:nrow(sta_df)) {
  station <- sta_df$STATION_KEY[i]
  df2 <- subset(df1, STATION_KEY == station)
  
  for(c in 1:length(count_it)) {
    metric <- ws_dfw_24mean[ws_dfw_24mean$STATION_KEY %in% station, count_it[c]]
    # calcualte what percent the daily mean precip is of the 24 hour metric
    df2[5 + c] <- df2$value / metric * 100
  }
  if(i==1) {
    df3 <- df2
  } else {
    df3 <- rbind(df3,df2)
  }
}

ws_dfl_dmean <- df3
rm(df1,df2,df3)

########################################################################################################
## write data to ouput
write.csv(sam_df, file=paste0(output_path,"R_output_Precip_samples_",thedate,".csv"),row.names = FALSE)
write.csv(ws_dfw_24mean, file=paste0(output_path,"R_output_dfw_24mean_24hour_Precip_metrics_",thedate,".csv"),row.names = FALSE)
write.csv(ws_dfw_dmean, file=paste0(output_path,"R_output_ws_dfl_dmean_precip_daily_mean_",thedate,".csv"),row.names=FALSE)
write.csv(ws_dfl_dmean, file=paste0(output_path,"R_output_ws_dfl_dmean_precip_daily_mean",thedate,".csv"),row.names=FALSE)


#########################################
save.image(file="precip_data.RData")
load("precip_data.RData")
#########################################

## Make a plot of daily rainfall for a station
library(lattice)
dmean.F <- ws_dfl_dmean[ws_dfl_dmean$STATION_KEY %in% c(11850),]
dmean.F$STATION_KEY <- as.factor(dmean.F$STATION_KEY)
dmean.F$date2 <- as.POSIXct(dmean.F$datestr) 

xax<-c(as.POSIXct("1994-01-01"),as.POSIXct("2014-01-01"))
yaxlab <- list(at=seq(0,125,25))
xaxlab <- list(at=c(as.POSIXct("1994-01-01"),
                    as.POSIXct("1995-01-01"),
                    as.POSIXct("1996-01-01"),
                    as.POSIXct("1997-01-01"),
                    as.POSIXct("1998-01-01"),
                    as.POSIXct("1999-01-01"),
                    as.POSIXct("2000-01-01"),
                    as.POSIXct("2001-01-01"),
                    as.POSIXct("2002-01-01"),
                    as.POSIXct("2003-01-01"),
                    as.POSIXct("2004-01-01"),
                    as.POSIXct("2005-01-01"),
                    as.POSIXct("2006-01-01"),
                    as.POSIXct("2007-01-01"),
                    as.POSIXct("2008-01-01"),
                    as.POSIXct("2009-01-01"),
                    as.POSIXct("2010-01-01"),
                    as.POSIXct("2011-01-01"),
                    as.POSIXct("2012-01-01"),
                    as.POSIXct("2013-01-01"),
                    as.POSIXct("2014-01-01")),
               labels=c("",
                        "1995","","","","",
                        "2000","","","","",
                        "2005","","","","",
                        "2010","","","",""))

xyplot(value ~ date2 | STATION_KEY, 
       data = dmean.F, type="l",
       xlab="Date", 
       ylab="Precipitaiton (mm)",
       scales=list(x=xaxlab,y=yaxlab))
