library(terra)
library(dplyr)
library(sf)
library(tigris)

# read in impervious lc
impervious <- rast("g:/BigData/fl_2022_ccap_v2_hires_impervious_20231226.tif")

fl <- tigris::counties(state = "Florida") %>% 
  filter(NAME == "Alachua")
fl <- st_transform(fl, crs = st_crs(impervious))


impervious_crop <-  crop(impervious,fl)
plot(impervious_crop)
impervious_alachua <- mask(impervious_crop, fl)

impervious_focal <- focal(impervious_alachua, w = 29, fun = "mean")
