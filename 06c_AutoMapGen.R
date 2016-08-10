library(ggmap)
library(ggplot2)
library(rgdal)
library(RgoogleMaps)
library(maptools)

streams <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'edges')
streams <- spTransform(streams, CRS("+proj=longlat +datum=WGS84"))
streams <- fortify(streams)

stns <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'sites')
stn <- stns[stns$STATION_KE == '21792',][1,]
stn <- spTransform(stn, CRS("+proj=longlat +datum=WGS84"))
stn <- data.frame(stn)

arca <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'LSN05_Watersheds')
stn_arca <- arca[arca$rid_LSN04 == unique(obs.complete[obs.complete$STATION_KEY == 21792, "rid_LSN04"]),]
stn_arca <- spTransform(stn_arca,  CRS("+proj=longlat +datum=WGS84"))
stn_arca <- fortify(stn_arca)

roadx <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'BLM_ROADX_SP')
roadx <- spTransform(roadx,  CRS("+proj=longlat +datum=WGS84"))
roadx <- data.frame(roadx)

wqlim <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", "Biocrit_Sed_303d")
wqlim <- spTransform(wqlim, CRS("+proj=longlat +datum=WGS84"))
wqlim <- fortify(wqlim)

dmas.shp <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", "DMAs")
dmas.shp <- spTransform(dmas.shp, CRS("+proj=longlat +datum=WGS84"))
dmas.shp@data$id <- rownames(dmas.shp@data)
dmas <- fortify(dmas.shp)
dmas <- merge(dmas, dmas.shp@data, by = 'id')

center <- c(mean(stn_arca$long), mean(stn_arca$lat))

zoom <- min(MaxZoom(range(stn_arca$lat),
                    range(stn_arca$long)))

gm <- get_googlemap(center = center, zoom = zoom, maptype = "roadmap")

gmap <- ggmap(gm, extent = "device") 
#gmap <- gmap + geom_path(data = streams, aes(x = long, y = lat, group = group), color = 'blue') 
gmap <- gmap + geom_path(data = stn_arca, aes(x = long, y = lat, group = group), colour = 'black')
gmap <- gmap + geom_point(data = stn, aes(x = coords.x1, y = coords.x2)) + 
  geom_text(data = stn, aes(x = coords.x1, y = coords.x2, label = STATION_KE),
            nudge_x = -.0018, nudge_y = .00019)
gmap_arca <- gmap + ggtitle("ARCA")
gmap_arca

gmap_roadx <- gmap_arca + geom_point(data = roadx, aes(x= coords.x1, y = coords.x2), colour = 'green')
gmap_roadx

gmap_wqlim <- gmap_arca + geom_path(data = wqlim, aes(x = long, y = lat, group = group), colour = 'orange')
gmap_wqlim

gmap <- ggmap(gm, extent = "normal", maprange = FALSE) 
#gmap <- gmap + geom_path(data = streams, aes(x = long, y = lat, group = group), color = 'blue') 
gmap <- gmap + geom_path(data = stn_arca, aes(x = long, y = lat, group = group), colour = 'black')
gmap <- gmap + geom_point(data = stn, aes(x = coords.x1, y = coords.x2)) + 
  geom_text(data = stn, aes(x = coords.x1, y = coords.x2, label = STATION_KE),
            nudge_x = -.0018, nudge_y = .00019)
gmap_arca <- gmap + ggtitle("ARCA")
gmap_arca
gmap_dmas <- gmap_arca + geom_polygon(data = dmas, aes(x = long, y = lat, group = group, fill = DMA_RP)) +
  coord_map(projection = 'mercator', xlim = range(stn_arca$long),
            ylim = range(stn_arca$lat)) + theme_nothing()
gmap_dmas
