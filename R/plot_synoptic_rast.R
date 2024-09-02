default_colours <- function(){rev(idpalette("iddu", 100))}

plot_synoptic_rast <- function(raster,
                               plot_name = NULL,
                               plot_desc = NULL,
                               colours = default_colours(),
                               scale_name = waiver(),
                               limits = NULL,
                               overlay_sites = FALSE,
                               site_data = NULL,
                               write_to_disk = TRUE,
                               filename = NULL
) {
  
  p <- ggplot() +
    geom_spatraster(
      data = raster
    ) + 
    scale_fill_gradientn(
      limits = limits,
      name = scale_name,
      na.value = "transparent",
      colours = colours)  +
    ggtitle(
      label = plot_name,
      subtitle = plot_desc
    ) +
    theme_snp_maps() 
  
  if (overlay_sites) {
    p <- p + 
      geom_point(aes(x = longitude, 
                     y = latitude), 
                 data = site_data,
                 shape = "+", 
                 inherit.aes = TRUE,
                 size = 2) 
  }

  if (write_to_disk) {
    ggsave(paste0("figures/",filename,".png"),
           p,
           width = 10, 
           height = 6,
           units = "in")
  } else {
    p
  }
  
}