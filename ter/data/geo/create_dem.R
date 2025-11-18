library(terra)

#Ter in before Ripoll
#dem <- rast("~/Documents/NEREIDA/main/ter/data/geo/dem_ter_up.tif")
dem_mask <- rast("~/Documents/NEREIDA/ter/data/geo/dem_ter_up.tif")
#basin <- vect("catchment_boundary.shp")
#dem_crop <- crop(dem, basin)
#dem_mask <- mask(dem_crop, basin)

# Compute hypsometric data
alt <- values(dem_mask, na.rm = TRUE)
breaks <- quantile(alt, seq(0, 1, by = 0.01)) # every 5%
HypsoData <- data.frame(
  HypsoPerc = seq(0, 100, by = 1),
  HypsoAlt  = as.numeric(breaks)
)

write.csv(HypsoData, "~/Documents/NEREIDA/main/ter/data/geo/hypso_ter_up.csv",
          quote = F, row.names = F)

#Ter in before Sau
dem_mask <- rast("~/Documents/NEREIDA/ter/data/geo/dem_ter_pasteral.tif")
#basin <- vect("catchment_boundary.shp")
#dem_crop <- crop(dem, basin)
#dem_mask <- mask(dem_crop, basin)

# Compute hypsometric data
alt <- values(dem_mask, na.rm = TRUE)
breaks <- quantile(alt, seq(0, 1, by = 0.01)) # every 5%
HypsoData <- data.frame(
  HypsoPerc = seq(0, 100, by = 1),
  HypsoAlt  = as.numeric(breaks)
)

write.csv(HypsoData, "~/Documents/NEREIDA/main/ter/data/geo/hypso_ter_pasteral.csv",
          quote = F, row.names = F)





#Ter in before Sau and after Ripoll
dem_mask <- rast("~/Documents/NEREIDA/ter/data/geo/dem_ter_afterRipoll.tif")
#basin <- vect("catchment_boundary.shp")
#dem_crop <- crop(dem, basin)
#dem_mask <- mask(dem_crop, basin)

# Compute hypsometric data
alt <- values(dem_mask, na.rm = TRUE)
breaks <- quantile(alt, seq(0, 1, by = 0.01)) # every 5%
HypsoData <- data.frame(
  HypsoPerc = seq(0, 100, by = 1),
  HypsoAlt  = as.numeric(breaks)
)

write.csv(HypsoData, "~/Documents/NEREIDA/main/ter/data/geo/hypso_ter_afteripoll.csv",
          quote = F, row.names = F)





#Ter in Gola
dem_mask <- rast("~/Documents/NEREIDA/ter/data/geo/dem_ter_down.tif")
#basin <- vect("catchment_boundary.shp")
#dem_crop <- crop(dem, basin)
#dem_mask <- mask(dem_crop, basin)

# Compute hypsometric data
alt <- values(dem_mask, na.rm = TRUE)
breaks <- quantile(alt, seq(0, 1, by = 0.01)) # every 5%
HypsoData <- data.frame(
  HypsoPerc = seq(0, 100, by = 1),
  HypsoAlt  = as.numeric(breaks)
)

write.csv(HypsoData, "~/Documents/NEREIDA/main/ter/data/geo/hypso_ter_down.csv",
          quote = F, row.names = F)




#Ter in Sau
dem_mask <- rast("~/Documents/NEREIDA/ter/data/geo/dem_ter_sau.tif")
#basin <- vect("catchment_boundary.shp")
#dem_crop <- crop(dem, basin)
#dem_mask <- mask(dem_crop, basin)

# Compute hypsometric data
alt <- values(dem_mask, na.rm = TRUE)
breaks <- quantile(alt, seq(0, 1, by = 0.01)) # every 5%
HypsoData <- data.frame(
  HypsoPerc = seq(0, 100, by = 1),
  HypsoAlt  = as.numeric(breaks)
)

write.csv(HypsoData, "~/Documents/NEREIDA/main/ter/data/geo/hypso_ter_sau.csv",
          quote = F, row.names = F)
