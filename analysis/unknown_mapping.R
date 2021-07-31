source(here("analysis/", "Functions.R"))
source(here::here("analysis/", "density_analysis.R"))


library(raster)
library(sf)
library(rgdal)
library(here)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(ggsn)



#-----------------------------------------------------------------
## Make the map
#-----------------------------------------------------------------


us <- raster::getData("GADM", country = c("United States"), level = 1)


clipper_small <- st_polygon(list(rbind(c(-120.7, 34.65),
                                       c(-119.25, 34.65), 
                                       c(-119.25, 33.8), 
                                       c(-120.7, 33.8), 
                                       c(-120.7, 34.65)))) %>% st_sfc() %>% st_set_crs(4326)


shore_small <- us %>% sf::st_as_sf() %>% sf::st_intersection(clipper_small) %>% sf::st_transform(4326) %>% sf::st_union()


patches <- readOGR(here("data/spatial", "patches.shp")) %>% st_as_sf() %>% st_set_crs("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs") %>% st_transform(4326) %>% st_buffer(dist = 0) %>% sf::st_intersection(clipper_small)


depth <- marmap::read.bathy("~/Downloads/sbc.xyz")
coords <- st_coordinates(clipper_small)
depth <- marmap::subsetBathy(depth, x = coords[,1], y = coords[,2], locator = F)
depth.df <- marmap::fortify.bathy(depth)


blues <- colorRampPalette(colors = c("#94AAC3", "#F9FAFB")) #Low #94AAC3, high #F9FAFB
browns <- colorRampPalette(colors = c("#ACD0A5", "#C3A76B"))



zoom_map <- ggplot()+
  geom_tile(data = filter(depth.df, z < 0), aes(x = x, y = y, fill = z))+
  scale_fill_gradientn(colours = blues(10))+
  geom_contour(data = filter(depth.df, z < 10), aes(x = x, y = y, z = z), color = "black", binwidth = 100, alpha = 0.25)+
  geom_sf(data = shore_small, fill = "#596778", lwd = 0.01)+
  geom_sf(data = patches, fill = alpha("#00802b", 0.9), col = alpha("#00802b", 0))+
  scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))+
  labs(x = "", y = "", size = expression(paste("Predicted\nconsumption (g m"^-2,"d"^-1,")")))+
  coord_sf(xlim = c(-120.7, -119.25), ylim = c(33.8, 34.65), expand = F)+
  annotation_scale(location = "bl", style = "ticks",  pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"))+
  annotation_north_arrow(location = "tr", style = north_arrow_nautical, height = unit(0.75, "cm"), width = unit(0.75, "cm"), pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"))+
  guides(color = FALSE, fill = FALSE, size = guide_legend(override.aes = list(shape = 21, color = "black", alpha = 1)))+
  theme(legend.title=element_text(size=14),legend.key=element_blank(),legend.background = element_blank())
#guides(color = FALSE, fill = FALSE, size = guide_legend(direction = "horizontal", label.position = "top", label.vjust = 0, override.aes = list(shape = 21)))+
#theme(legend.position = c(0.4, 0.9), panel.border = element_rect(colour = "black", fill=NA, size=1))

ggsave(filename = here::here("figures/ggmap-bd.png"),plot = zoom_map, device = "png", width = 12*0.85, height = 8*0.85)





clipper_large <- st_polygon(list(rbind(c(-135, 21),
                                       c(-135, 62), 
                                       c(-105, 62), 
                                       c(-105, 21), 
                                       c(-135, 21)))) %>% st_sfc() %>% st_set_crs(4326)

world1 <- sf::st_as_sf(maps::map('world', region = c("USA", "Canada", "Mexico"), plot = FALSE, fill = TRUE)) %>% st_transform(4326)  %>% st_union() %>% st_intersection(clipper_large)


map <- ggplot() + 
  geom_sf(data = world1, fill = "#C2BDBD", col = "#C2BDBD") + 
  geom_sf(data = clipper_small, color = "red", fill = NA, lwd = 1)+
  coord_sf(xlim = c(-135, -105), ylim = c(22,61), expand = F)+
  scale_x_continuous(breaks = c(-110, -120, -130))+
  scale_y_continuous(breaks = c(seq(25, 55, by = 15)))+
  theme(axis.text = element_text(size = 8), panel.background = element_blank(), axis.line = element_line())

ggsave(filename = here::here("figures/ggmap_big-bd.png"),plot = map, device = "png", width = 12*0.25, height = 8*0.85 )

#------------------------------------
## Adrian zoom
#------------------------------------


clipper_AS <- st_polygon(list(rbind(c(-121.506048, 35),
                                    c(-106.374163, 35), 
                                    c(-106.374163, 21.767804), 
                                    c(-121.506048, 21.767804), 
                                    c(-121.506048, 35)))) %>% st_sfc() %>% st_set_crs(4326)


# world1 <- sf::st_as_sf(maps::map('world', region = c("USA", "Canada", "Mexico"), plot = FALSE, fill = TRUE)) %>% st_transform(4326)  %>% st_union() %>% st_intersection(clipper_AS)

mex <- raster::getData("GADM", country = c("Mexico"), level = 1) %>% sf::st_as_sf() %>% sf::st_intersection(clipper_AS) %>% sf::st_transform(4326) 

shore <- us %>% sf::st_as_sf() %>% sf::st_intersection(clipper_AS) %>% sf::st_transform(4326) %>% sf::st_union(mex)

bathy <- marmap::getNOAA.bathy(lon1 = -121.506048, lon2 = -106.374163,
                               lat1 = 21.767804 , lat2 = 35, resolution = 1)

depth.df <- marmap::fortify.bathy(bathy)
breaks <- c(-4000, -3000, -2000, -1000, -500, -100)

p.wide <- ggplot() + 
  geom_tile(data = filter(depth.df, z < 0), aes(x = x, y = y, fill = z), show.legend = F)+
  scale_fill_gradientn(colours = blues(10))+
  geom_contour(data = filter(depth.df, z < 0), aes(x = x, y = y, z = z), color = "black", alpha = 0.25, breaks = breaks)+
  geom_sf(data = shore, fill = "#596778", color = "#596778") + 
  coord_sf(expand = F)+
  geom_sf(data = clipper_small, color = "red", fill = NA, lwd = 1)+
  labs(x = "", y = "")+
  theme(axis.text = element_text(size = 12), panel.background = element_blank(), axis.line = element_line())

ggsave(filename = here::here("figures/ggmap_wide-BD.png"),plot = p.wide, device = "png", width = 12*0.85, height = 8*0.85)

