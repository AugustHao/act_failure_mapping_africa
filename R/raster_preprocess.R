
#' for a folder of MAP rasters, calculate the mean across years, and crop to
#' africa
#'
#' @param filepath where the raster files live
#' @param remove_every_second_layer if every second layer is to be removed, to
#'   deal with some quirks in MAP layers
#'
#' @return a single layer raster file
#' @export
#'
#' @examples
make_synoptic_mean_rast <- function(filepath = "data/rasters/", 
                                    name,
                                    remove_every_second_layer = TRUE) {
  
  r_stack <- terra::rast(list.files(paste0(filepath,name,"/"),full.names = TRUE))
  africa <- terra::rast("data/rasters/africa_mask.tif")
  # renive the second layer, some LatAm specific thing?
  if (remove_every_second_layer) {
    r_stack <- r_stack[[seq(1,nlyr(r_stack),by = 2)]]
  }
  
  stack_africa <- terra::crop(r_stack,africa)
  stack_africa <- terra::mask(stack_africa,africa)
  africa_mean <- terra:::mean(stack_africa)
  names(africa_mean) <- name
  return(africa_mean)
}


