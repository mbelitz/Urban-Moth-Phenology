library(terra)
library(dplyr)
library(sf)
library(tigris)
library(tidyterra)
library(ggplot2)
library(tmap)

terraOptions(memfrac=0.9, tempdir = "temp")

# read in impervious lc
impervious <- rast("fl_2022_ccap_v2_hires_impervious_20231226.tif")

fl <- tigris::counties(state = "Florida") %>% 
  filter(NAME == "Alachua")
fl <- st_transform(fl, crs = st_crs(impervious))

# site location
# gather sites
sites <- vect("data/Site Locations.shp") %>% 
  project(impervious)

sites_buffer_300 <- st_buffer(st_as_sf(sites), dist = 300)
st_bbox(sites_buffer_300)
bbox <- st_bbox(c(xmin = 1301240.7, ymin = 814903.5, 
                  xmax = 1333952.3  , ymax = 836695.6), 
                crs = st_crs(sites_buffer_300))
# crop landcover 
lc_crop_300 <- crop(impervious, bbox)
plot(lc_crop_300)

impervious <- c(1)
notImpervious <- c(0)
class_matrix <- cbind(impervious,1) %>% 
  rbind(cbind(notImpervious,0))

lc_crop_300_rc <- classify(lc_crop_300,
                           rcl = class_matrix,
                           others = NA)

plot(lc_crop_300_rc)
lc_crop_300_rc

library(ggspatial)
library(tidyterra)

imperviousSurfacePlot <- ggplot() +
  geom_spatraster(lc_crop_300_rc, mapping = aes()) +
  geom_sf(sites_buffer_300, mapping = aes(),fill = NA, color = "#3A0CA3", linewidth =1, alpha = 0.6) + 
  scale_fill_gradient(breaks = c(0,1), low = "white", high = "#7A7D7D") +
  coord_sf() +
  guides(fill = guide_legend()) +
  ggtitle("Impervious surface") +
  annotation_scale() + # add scale
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=13.5))

imperviousSurfacePlot

################################################################################
canopy <- rast("fl_2022_ccap_v2_hires_canopy_20231226.tif")

# crop land cover
lc_canopy_300 <- crop(canopy, bbox)
forest <- c(1)
notForest <- c(0,2)
canopy_matrix <- cbind(forest,1) %>% 
  rbind(cbind(notForest,0))

lc_canopy_300 <- classify(lc_canopy_300,
                          rcl = canopy_matrix,
                          others = NA)
plot(lc_canopy_300)

canopyCoverPlot <- ggplot() +
  geom_spatraster(lc_canopy_300, mapping = aes()) +
  geom_sf(sites_buffer_300, mapping = aes(),fill = NA, color = "#3A0CA3", linewidth =1, alpha = 0.6) + 
  scale_fill_gradient(breaks = c(0,1), low = "white", high = "#5FB05A") +
  coord_sf() +
  guides(fill = guide_legend()) +
  ggtitle("Canopy cover") +
  annotation_scale() + # add scale
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=13.5))
canopyCoverPlot