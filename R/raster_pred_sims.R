raster_pred_sims <- function(focus_ext = TRUE,
                             sample_raster,
                             latents_obs = model$model_arrays$gp,
                             beta = model$model_arrays$parameters$beta,
                             loadings = model$model_arrays$parameters$loadings,
                             covariates = covariate_rast,
                             draws = draws
                             ) {
  
  if (focus_ext) {
    # zoom in to data rich region
    new_ext <- ext(25,37,-4,6)
    focus_rast <- sample_raster %>% 
      crop(new_ext)
    coords_pixel <- terra::xyFromCell(focus_rast,
                                      cell = terra::cells(focus_rast))
  } else {
    coords_pixel <- terra::xyFromCell(sample_raster,
                                      cell = terra::cells(sample_raster))
  }
  
  # get latent factors and design matrix
  latents_pixel <- greta.gp::project(latents_obs, coords_pixel)
  X_pixel <- build_design_matrix(covariates, coords_pixel)
  
  if (anyNA(X_pixel)) {
    stop("NAs in design matrix, check gap in rasters")
  }
  
  # predict SNP frequencies
  snp_freq_logit_pixel <- X_pixel %*% t(beta) +
    t(loadings %*% t(latents_pixel))
  snp_freq_pixel <- ilogit(snp_freq_logit_pixel)
  
  # # compute posterior samples of SNPs
  post_pixel_sims <- calculate(snp_freq_pixel,
                               values = draws,
                               nsim = 1e2)
  
  # return mean at pixel level 
  apply(post_pixel_sims$snp_freq_pixel,2:3,mean)
  
}