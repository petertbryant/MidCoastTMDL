library(sp)
library(rgdal)
library(plyr)
library(foreign)
library(raster)
library(rgeos)
library(maptools)

options(scipen = 100)

#### Load Data ####
#All stations in SSN
stns <- readShapeSpatial('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/ssn/stations/stations_ssn_20140822_gcs_nad83.shp', 
                         proj4string=CRS('+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs'))
stns <- spTransform(stns, CRS("+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999984 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048 +no_defs"))

#the huc 10 field is truncated for some reason
stns$HU_10 <- substr(stns$HU_12,1,10)

#Midcoast specific elevation - 
mc_elev <- raster("//deqhq1/tmdl/tmdl_wr/midcoast/gis/raster/dem/midcoast_10m/dem10m_midcst")
crs(mc_elev) <- "+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999984 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048 +no_defs"

#Statewide precipitation 30 year annual average over time period 1971 to 2000
precip <- raster("//deqhq1/tmdl/tmdl_wr/midcoast/gis/raster/climate/precipitation_annual_1971-2000_mmx100/ppt_71-00_lam")
crs(precip) <- "+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999984 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048 +no_defs"

#Statewide fine sediment potential geology classification 
fsp <- raster("//deqhq1/tmdl/tmdl_wr/midcoast/gis/raster/geology_erod_res/fsp/fsp_state")
crs(fsp) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
pct_fsp <- readShapeSpatial('//deqhq1/tmdl/tmdl_wr/midcoast/gis/shapefiles/sediment/target/or_huc_5th_fsp_dissolve.shp',
                            proj4string=CRS('+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999984 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048 +no_defs'))

#NHD attributes
flow <- read.dbf("//deqhq1/tmdl/tmdl_wr/midcoast/gis/NHD_attributes/NHDPlus17/flowlineattributesflow.dbf")
nhd_precip <- read.dbf("//deqhq1/tmdl/tmdl_wr/midcoast/gis/NHD_attributes/NHDPlus17/flowlineattributestempprecip.dbf")

#FSS values 
fss <- read.xlsx("//deqhq1/tmdl/tmdl_wr/midcoast/data/benthicmacros/data_from_ctsi/Siletz tribes_2012 PREDATOR bio indices.xlsx", sheetName = 'Fines')
fss[fss$CTSI.ID == 'CTSI_57','CTSI.ID'] <- 29898
fss$FSS <- round(fss$FSS)

#### Derive attributes ####
#These stations are on the SSN network and not the NHD but for the CTSI stations the max distance to the NHD is 40 feet. 
#I don't think this distance will affect the elevation too dramatically so we should be good.
stns_e <- extract(mc_elev, stns, sp = TRUE)

#Precip data needs to be divided by 100
stns_ep <- extract(precip, stns_e, sp = TRUE)
stns_ep$ppt_71.00_lam <- stns_ep$ppt_71.00_lam/100

stns_epf <- merge(stns_ep, pct_fsp@data, by.x = 'HU_10', by.y = "WATERSHED")
stns_epff <- extract(fsp, stns_epf, sp = TRUE)
stns_epff$ERODRES_CAT <- ifelse(stns_epff$fsp_state == 3,'E','R')

stns_attr <- merge(stns_epff, flow, by.x = 'NHDP21_COM', by.y = 'COMID', all.x = TRUE)
stns_attr <- merge(stns_attr, nhd_precip, by.x = 'NHDP21_COM', by.y = 'COMID', all.x = TRUE)
stns_attr$STRMPOWER <- stns_attr$SLOPE * stns_attr$MAFLOWU

stns_attr <- merge(stns_attr, fss[,c('CTSI.ID','FSS')], by.x = 'STATION_KE', by.y = 'CTSI.ID', all.x = TRUE)

CTSI_chars <- rename(stns_attr, c('STATION_KE' = 'STATION_KEY','ppt_71.00_lam' = 'PRECIPSITE_MM', 'LONG_FACv2' = 'LONG_NHD', 
                                 'LAT_FACv21' = 'LAT_NHD', 'dem10m_midcst' = 'ELEV_FT', 'FSP3_PCT' = 'PCT_FSP3',
                                 'FSS' = 'FSS_May05'))

rm(list = setdiff(ls(), c('CTSI_chars')))

######################################################
#### Checking the additional data by comparing to ####
 ### previously used CART input data              ###
  ##                                              ##
   #                                              #

# #These are the input data for the CART calculated for ALL samples (except CTSI sites)
# chars <- read.csv('//deqhq1/tmdl/tmdl_wr/midcoast/data/benthicmacros/stationwork/r_inputs/r_input_cart_2013_06_15.csv')
# chars$SVN <- str_trim(chars$SVN)
# 
# #checking additional data additions
# check <- merge(stns_attr@data, chars, by.x = 'STATION_KE', by.y = 'STATION_KEY', all.x = TRUE)
# 
# #NHD Attrs
# check[,c('LAT_NHD','LAT_FACv21')] #pretty close
# check[,c('AREAWTMAP.x','AREAWTMAP.y')] #NHD attrs look good
# 
# #Preciptitation
# check[,c('ppt_71.00_lam','PRECIPSITE_MM')] #off by quite a bit in some locations - this is due to the points snapped to SSN and not to NHD
# check$precip_resid <- check[,'ppt_71.00_lam'] - check[,'PRECIPSITE_MM']
# plot(check[!is.na(check$precip_resid) & grepl('Siletz',check$HU_8_NAME.x),'precip_resid'])
# 
# #Elevation
# check[,c('dem10m_midcst','ELEV_FT')]
# check$elev_resid <- check[,'dem10m_midcst'] - check[,'ELEV_FT']
# plot(check$elev_resid) #slight diff due to moving stations
# 
# #Geology - all good for CTSI sites. some oddities for some of the other sites but outside the scope of this update
# all(check[,'FSP3_PCT'] == check[,'PCT_FSP3'], na.rm = TRUE)
