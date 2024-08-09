plot_synoptic_rast <- function(raster,
                               name = NULL,
                               desc = NULL,
                               overlay_sites = FALSE,
                               site_data = NULL,
                               write_to_disk = TRUE,
                               filename = NULL
) {
  
  p <- ggplot() +
    geom_spatraster(
      data = raster
    ) +
    #facet_wrap(~lyr, nrow = 1, ncol = 3) +
    scale_fill_gradientn(
      #labels = scales::percent,
      # name = "value",
      limits = c(0, 1),
      na.value = "transparent",
      colours = rev(idpalette("iddu", 100)))  +
    ggtitle(
      label = name,
      subtitle = desc
    ) +
    theme_snp_maps() 
  
  if (overlay_sites) {
    p <- p + 
      geom_point(aes(x = longitude, 
                     y = latitude), 
                 data = snp_data,
                 #col = "white",
                 shape = "+", 
                 inherit.aes = TRUE,
                 size = 1) 
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