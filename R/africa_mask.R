africa_mask <- function(){
  Pf_5k <- terra::rast("data/rasters/MAP_Regions_Pf_5k.tif")
  
  Pf_5k_africa <- (Pf_5k == 1)
  Pf_5k_africa[Pf_5k_africa == 0] <- NA
  Pf_5k_africa[] <- as.integer(Pf_5k_africa[])
  africa_ext <- ext(-33,63,-41,38)
  Pf_5k_africa <- terra::crop(Pf_5k_africa,africa_ext)
 
  Pf_5k_africa 
}