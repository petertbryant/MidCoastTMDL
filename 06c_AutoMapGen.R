library(ggmap)
library(ggplot2)
library(rgdal)
library(RgoogleMaps)
library(maptools)
library(rgeos)
library(raster)
library(ggsn)

#### Bring in the data files for mapping and plotting ####
# streams <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'edges')
# streams <- spTransform(streams, CRS("+proj=longlat +datum=WGS84"))
# streams <- fortify(streams)

stns_shp <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'sites')

arca_shp <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'LSN05_Watersheds')

roadx_shp <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", 'BLM_ROADX_SP')
roadx_shp <- spTransform(roadx_shp,  CRS("+proj=longlat +datum=WGS84"))
roadx_shp <- data.frame(roadx_shp)
roadx_shp$color <- 'Road crossing'

wqlim_shp <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", "Biocrit_Sed_303d")
wqlim_shp <- spTransform(wqlim_shp, CRS("+proj=longlat +datum=WGS84"))
wqlim_shp@data$id <- rownames(wqlim_shp@data)

dmas_shp <- readOGR("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping", "DMAs")
dmas_shp <- spTransform(dmas_shp, CRS("+proj=longlat +datum=WGS84"))

# pop_ras <- raster("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping/coast_dasy")
# pop_ras <- projectRaster(pop_ras, crs = CRS("+proj=longlat +datum=WGS84"))

# suscep_ras <- raster("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping/smc")
# suscep.poly <- rasterToPolygons(suscep_ras, dissolve = TRUE)

#### Define station area #### 
stn <- stns_shp[stns_shp$STATION_KE == 26822,][1,]
stn <- spTransform(stn, CRS("+proj=longlat +datum=WGS84"))
stn <- data.frame(stn)

stn_arca.shp <- arca_shp[arca_shp$rid_LSN04 == unique(
  ssn_RF_data[ssn_RF_data$STATION_KEY == 26822, "rid"]),]
stn_arca.shp <- spTransform(stn_arca.shp,  
                            CRS("+proj=longlat +datum=WGS84"))
stn_arca <- fortify(stn_arca.shp)

wqlim_shp <- wqlim_shp[stn_arca.shp,]

#### Get the base map ####
center <- c(mean(stn_arca$long), mean(stn_arca$lat))

zoom <- min(MaxZoom(range(stn_arca$lat),
                    range(stn_arca$long)))

gm <- get_googlemap(center = center, zoom = zoom, maptype = "roadmap")

#### Set up the station and arca plot ####
gmap <- ggmap(gm, extent = "device") 
gmap <- gmap + geom_path(data = stn_arca, 
                         aes(x = long, 
                             y = lat, 
                             group = group), 
                         colour = 'black')
gmap <- gmap + geom_point(data = stn, 
                          aes(x = coords.x1, 
                              y = coords.x2,
                              color = factor(STATION_KE)),
                          size = 3) +
  scale_color_manual(values = "black") + 
  theme(legend.title = element_blank())
# geom_text(data = stn, aes(x = coords.x1, 
#                           y = coords.x2, 
#                           label = STATION_KE),
#           nudge_x = -.0018, nudge_y = .00019) 
gmap_arca <- gmap + ggtitle("ARCA") + 
  scalebar(data = stn_arca, 
           height = 0.02, 
           dist = 1, 
           dd2km = TRUE, 
           model = 'WGS84', 
           st.size = 3, 
           location = 'topright', 
           anchor = c(x = (attr(gm, "bb")$ur.lon - .001), 
                      y = (attr(gm, "bb")$ur.lat - .001)))
gmap_arca

#### Road crossings ####
gmap_roadx <- gmap_arca + 
  geom_point(data = roadx_shp, 
             aes(x= coords.x1, 
                 y = coords.x2,
                 colour = color),
             size = 3) + 
  ggtitle("Road Crossings") + 
  scale_color_manual(values = c('black','green')) + 
  theme(legend.title = element_blank())
gmap_roadx

#### WQ listings ####
if (nrow(wqlim_shp@data) > 0) {
  wqlim <- fortify(wqlim_shp)
  wqlim <- merge(wqlim, wqlim_shp@data, by = 'id')
  wqlim <- plyr::rename(wqlim, c('Listing_St' = "WQ_limited"))
  #wqlim <- (wqlim[!duplicated(wqlim$Pollutant),])
  
  gmap_wqlim <- gmap_arca + geom_path(data = wqlim, 
                                      aes(x = long, 
                                          y = lat, 
                                          group = group,
                                          color = WQ_limited)) +
    ggtitle("303(d) Listings") + 
    scale_color_manual(values = c('orange','black','blue')) +
    theme(legend.title = element_blank())
} else {
  gmap_wqlim <- gmap_arca + ggtitle("No 303(d) Listings in the ARCA") +
    theme(legend.title = element_blank())
}
gmap_wqlim

#### Set up base map for polygon and raster mapping ####
gm <- get_googlemap(center = center, zoom = zoom, maptype = "roadmap")
gmap_poly <- ggmap(gm, extent = "normal", maprange = FALSE) 
gmap_poly <- gmap_poly + geom_path(data = stn_arca, aes(x = long, 
                                                        y = lat, 
                                                        group = group), 
                                   colour = 'black')
gmap_poly <- gmap_poly + geom_point(data = stn, aes(x = coords.x1, 
                                                    y = coords.x2,
                                                    color = factor(STATION_KE)),
                                    size = 3) + 
  # geom_text(data = stn, aes(x = coords.x1, 
  #                           y = coords.x2, 
  #                           label = STATION_KE),
  #           nudge_x = -.0018, 
  #           nudge_y = .00019) + 
  scale_color_manual(values = 'black') + 
  theme(legend.title = element_blank()) +
  scalebar(data = stn_arca, 
           height = 0.02, 
           dist = 1, 
           dd2km = TRUE, 
           model = 'WGS84', 
           st.size = 3, 
           location = 'topright', 
           anchor = c(x = (attr(gm, "bb")$ur.lon - .01), 
                      y = (attr(gm, "bb")$ur.lat - .01)))

bb_poly <- as(extent(sort(as.vector(as.matrix(attr(gm, "bb"))))), 
              "SpatialPolygons")
proj4string(bb_poly) <- CRS("+proj=longlat +datum=WGS84")

#### DMAs ####
dmas_sub <- raster::intersect(dmas_shp, bb_poly)
dmas_sub@data$id <- rownames(dmas_sub@data)
dmas <- fortify(dmas_sub)
dmas <- merge(dmas, dmas_sub@data, by = 'id')
dmas <- plyr::rename(dmas, c('DMA_RP' = "Designated_Management_Authority"))

gmap_dmas <- gmap_poly + geom_polygon(data = dmas, 
                                      aes(x = long, 
                                          y = lat, 
                                          group = group, 
                                          fill = Designated_Management_Authority), 
                                      alpha = 0.6) + 
  theme(axis.line = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank()) + 
  ggtitle("DMAs") 
gmap_dmas

#### Population ####
pop.ras.sub <- crop(pop_ras, extent(bb_poly), snap = "out")
pop.poly <- rasterToPolygons(pop.ras.sub)
pop.poly@data$id <- 1:nrow(pop.poly@data)
pop.poly.fort <- fortify(pop.poly, data = pop.poly@data)
pop.poly.fort <- merge(pop.poly.fort, pop.poly@data, by = 'id')
if (!all(pop.poly.fort$coast_dasy == 0)) {
  pop.poly.fort[pop.poly.fort$coast_dasy == 0, 'coast_dasy'] <- NA   
  
  gmap_pop <- gmap_poly + geom_polygon(data = pop.poly.fort, 
                                       aes(x = long, 
                                           y = lat, 
                                           group = group, 
                                           fill = coast_dasy), 
                                       alpha = 0.5, 
                                       size = 0) +  ## size = 0 to remove the polygon outlines
    scale_fill_gradientn(colours = topo.colors(255)) + 
    theme(axis.line = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_blank(), 
          axis.ticks = element_blank()) + 
    ggtitle("Population density") 
} else {
  gmap_pop <- gmap_poly + 
    theme(axis.line = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_blank(), 
          axis.ticks = element_blank()) + 
    ggtitle("Population density is 0 in this ARCA") 
}
#gmap_pop

# #### Susceptibility ####
# gm <- get_googlemap(center = center, zoom = zoom, maptype = "terrain")
# gmap_poly <- ggmap(gm, extent = "normal", maprange = FALSE) 
# gmap_poly <- gmap_poly + geom_path(data = stn_arca, aes(x = long, 
#                                                         y = lat, 
#                                                         group = group), 
#                                    colour = 'black')
# gmap_poly <- gmap_poly + geom_point(data = stn, aes(x = coords.x1, 
#                                                     y = coords.x2,
#                                                     color = factor(STATION_KE)),
#                                     size = 3) + 
#   # geom_text(data = stn, aes(x = coords.x1, 
#   #                           y = coords.x2, 
#   #                           label = STATION_KE),
#   #           nudge_x = -.0018, 
#   #           nudge_y = .00019) + 
#   scale_color_manual(values = 'black') + 
#   theme(legend.title = element_blank()) +
#   scalebar(data = stn_arca, 
#            height = 0.02, 
#            dist = 1, 
#            dd2km = TRUE, 
#            model = 'WGS84', 
#            st.size = 3, 
#            location = 'topright', 
#            anchor = c(x = (attr(gm, "bb")$ur.lon - .01), 
#                       y = (attr(gm, "bb")$ur.lat - .01)))
# 
# bb_poly_suscep <- spTransform(bb_poly, crs(suscep.ras))
# suscep.ras.sub <- crop(suscep.ras, extent(bb_poly_suscep))
# suscep.ras.sub <- projectRaster(suscep.ras.sub, 
#                                 crs = CRS("+proj=longlat +datum=WGS84"),
#                                 method = "ngb")
# suscep.poly <- rasterToPolygons(suscep.ras.sub, dissolve = TRUE)
# suscep.poly <- raster::intersect(suscep.poly, bb_poly)
# suscep.poly@data$id <- 1:nrow(suscep.poly@data)
# suscep.poly.fort <- fortify(suscep.poly)
# suscep.poly.fort <- merge(suscep.poly.fort, suscep.poly@data, by = 'id')
# suscep.poly.fort$suscep <- as.factor(suscep.poly.fort$suscep)
# 
# gmap_suscep <- gmap_poly + geom_polygon(data = suscep.poly.fort, 
#                                         aes(x = long, 
#                                             y = lat, 
#                                             group = group, 
#                                             fill = suscep), 
#                                         alpha = 0.3, 
#                                         size = 0) +  ## size = 0 to remove the polygon outlines
# scale_fill_manual(values = c('1' = '#00F500', '2' = '#94F700', '3' = '#F5F500',
#                                '4' = '#FC8B00', '5' = "#F50000")) + 
#   theme(axis.line = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         panel.background = element_blank(), 
#         axis.text = element_blank(),
#         axis.title = element_blank(), 
#         axis.ticks = element_blank()) + 
#   ggtitle("Landslide Susceptibility") 
# gmap_suscep

#### Disturbance ####
yrs <- unique(ssn_RF_data[ssn_RF_data$STATION_KEY == 26822, "YEAR"])
dis_map_list <- list()
for (i in 1:length(yrs)) {
  if (any(yrs > 2008)) {
    yrs[yrs > 2008] <- 2008
    yrs <- unique(yrs)
  }
  dis.ras <- raster(paste0("C:/Users/pbryant/Desktop/MidCoastTMDL-GIS/Shapefiles_for_mapping/dis_1YR_", 
                           yrs[i]))
  bb_poly_dis <- spTransform(bb_poly, crs(dis.ras))
  dis.ras.sub <- crop(dis.ras, extent(bb_poly_dis), snap = "out")
  dis.ras.sub <- projectRaster(dis.ras.sub, 
                               crs = CRS("+proj=longlat +datum=WGS84"))
  dis.poly <- rasterToPolygons(dis.ras.sub)
  dis.poly <- raster::intersect(dis.poly, bb_poly)
  dis.poly@data$id <- 1:nrow(dis.poly@data)
  dis.poly.fort <- fortify(dis.poly, data = dis.poly@data)
  dis.poly.fort <- merge(dis.poly.fort, dis.poly@data, by = 'id')
  dis.poly.fort[dis.poly.fort$dis_1YR_2003 == 0, paste0('dis_1YR_', 
                                                        yrs[i])] <- NA 
  
  gm <- get_googlemap(center = center, zoom = zoom, maptype = "terrain")
  gmap_poly <- ggmap(gm, extent = "normal", maprange = FALSE) 
  gmap_poly <- gmap_poly + geom_path(data = stn_arca, aes(x = long, 
                                                          y = lat, 
                                                          group = group), 
                                     colour = 'black')
  gmap_poly <- gmap_poly + geom_point(data = stn, aes(x = coords.x1, 
                                                      y = coords.x2)) + 
    geom_text(data = stn, aes(x = coords.x1, 
                              y = coords.x2, 
                              label = STATION_KE),
              nudge_x = -.0018, 
              nudge_y = .00019) + scale_color_continuous(guide = FALSE)
  
  
  gmap_dis <- gmap_poly + geom_polygon(data = dis.poly.fort, 
                                       aes_string(x = "long", 
                                                  y = "lat", 
                                                  group = "group", 
                                                  fill = paste0("dis_1YR_", 
                                                                yrs[i])), 
                                       alpha = 0.9, 
                                       size = 0) +  
    theme(axis.line = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_blank(), 
          axis.ticks = element_blank()) + 
    ggtitle("Disturbance - 1 Year Prior") 
  dis_map_list[[i]] <- gmap_dis
}
