library(sp)
library(rgdal)
library(plyr)
library(foreign)
library(raster)
library(rgeos)
library(maptools)

options(scipen = 100)

stns <- readShapeSpatial('//deqhq1/tmdl/tmdl_wr/midcoast/models/sediment/ssn/stations/stations_ssn_20140822_gcs_nad83.shp')
mc_elev <- raster("//deqhq1/tmdl/tmdl_wr/midcoast/gis/raster/dem/midcoast_10m/dem10m_midcst")

extract(mc_elev,stns, sp = TRUE)
