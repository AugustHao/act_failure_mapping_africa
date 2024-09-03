gp_check <- function(spatial_latents,
                     target_raster,
                     posterior_sims = NULL,
                     plot_var = FALSE) {
  
  coords_pixel <- terra::xyFromCell(target_raster,
                                    cell = terra::cells(target_raster))
  
  latents_pixel <- greta.gp::project(spatial_latents, coords_pixel)
  
  post_pixel_sims <- calculate(latents_pixel,
                               nsim = 1e3,
                               values = posterior_sims)
  
  post_mean_pixel <- apply(post_pixel_sims[[1]],2:3,mean)
  post_var_pixel <- apply(post_pixel_sims[[1]],2:3,var)
  
  post_mean_rast <- target_raster * 0
  names(post_mean_rast) <- "prior mean"
  
  post_var_rast <- target_raster * 0
  names(post_var_rast) <- "prior variance"
  
  n_latents <- ncol(spatial_latents)
  
  png(paste0("figures/gp_",
             ifelse(is.null(posterior_sims),
                    "prior",
                    "posterior"),
             "_check.png"),
      width = 250 * n_latents, height = 300)
  
  if (plot_var) {
    par(mfcol = c(2,n_latents))
    for (i in seq_len(n_latents)) {
      post_mean_rast[terra::cells(post_mean_rast)] <- post_mean_pixel[,i]
      plot(post_mean_rast, main = paste0("latent group mean ",i))
      post_var_rast[terra::cells(post_var_rast)] <- post_var_pixel[,i]
      plot(post_var_rast, main = paste0("latent group variance ",i))
    }
  } else {
    par(mfcol = c(1,n_latents))
    for (i in seq_len(n_latents)) {
      post_mean_rast[terra::cells(post_mean_rast)] <- post_mean_pixel[,i]
      plot(post_mean_rast, main = paste0("latent group mean ",i))
    }
  }
  dev.off()
}