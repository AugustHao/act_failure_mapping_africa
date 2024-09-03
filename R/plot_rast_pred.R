plot_rast_pred <- function(focus_ext = TRUE,
                           sample_raster,
                           snp_table,
                           pixel_sims,
                           snp_data = snp_data) {
  
  
  if (focus_ext) {
    # zoom in to data rich region
    new_ext <- ext(25,37,-4,6)
    sample_raster <- sample_raster %>% 
      crop(new_ext)
    }
  
  snp_freq_mean_rast <- sample_raster * 0
  
  # plot focus
  outline_vector <- gadm(
    country = c("Uganda",
                "Rwanda",
                "Burundi",
                "Kenya",
                "Tanzania",
                "South Sudan",
                "Ethiopia"),
    level = 0,
    path = "~/not_synced/africa/"
  )
  
  outline_vector <- crop(outline_vector,sample_raster)
  outline_vector <- st_cast(st_as_sf(outline_vector),"LINESTRING")
  
  for (i in snp_table$id[snp_table$valid == TRUE]) {
    
    snp_freq_mean_rast[terra::cells(snp_freq_mean_rast)] <- pixel_sims[,i]
    names(snp_freq_mean_rast) <- "mean"
    
    plot_dat <- snp_data %>% 
      filter(snp_id == i) %>% 
      mutate(prevalence = snp_count/sample_size) %>% 
      group_by(longitude,latitude) %>% 
      summarise(prevalence = mean(prevalence))
    
    plot_dat <- plot_dat[!is.na(
      terra::extract(sample_raster,
                     cbind(plot_dat$longitude,plot_dat$latitude))),
    ]
    
    # plot_dat <- plot_dat %>% 
    #   mutate(across(longitude:latitude, \(x) jitter(x, factor = 45))) 
    
    p <- ggplot() +
      geom_spatraster(
        data = snp_freq_mean_rast
      ) +
      scale_fill_gradientn(
        name = "mean prevalence",
        limits = c(0, 0.5),
        na.value = "transparent",
        colours = rev(idpalette("idem", 100)))  +
      ggtitle(
        label = snp_table$snp[i],
        subtitle = "prevalence in artemisinin resistance marker data"
      ) +
      theme_snp_maps()
    
    if (focus_ext) {
      p <- p + 
        geom_spatvector(data = outline_vector,
                        linewidth = 0.2,
                        col = "grey")
      
      ggsave(paste0("figures/",
                    snp_table$snp[i],
                    "_pred_focus.png"),
             p,
             width = 10, height = 6, units = "in")
      
      # plot data on its own
      ggplot() + 
        geom_spatvector(data = outline_vector,
                        linewidth = 0.2,
                        col = "grey") +
        theme_snp_maps() + 
        geom_point(aes(x = longitude,
                       y = latitude,
                       col = prevalence),
                   data = plot_dat,
                   #col = "white",
                   shape = "+",
                   inherit.aes = TRUE,
                   size = 10) +
        scale_color_gradientn(name = "observed prevalence",
                              limits = c(0, 0.5),
                              na.value = "transparent",
                              colours = alpha(rev(idpalette("iddu", 100)[80:100]),1)) +
        ggtitle(
          label = snp_table$snp[i],
          subtitle = "observed prevalence in sites (averaged across samples and years)"
        )
      
      ggsave(paste0("figures/",
                    snp_table$snp[i],
                    "_data_focus.png"),
             width = 10, height = 6, units = "in")
      
    } else {
      ggsave(paste0("figures/",
                    snp_table$snp[i],
                    "_pred.png"),
             p,
             width = 10, height = 6, units = "in")
    }

  }
}