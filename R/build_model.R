build_snp_model <- function(snp_data, 
                            n_latent = 4,
                            n_snp = n_snp,
                            X_obs = X_obs,
                            coords = snp_coords) {
  
  # define model parameters
  parameters <- define_greta_parameters(n_snp = n_snp,
                                        n_latent = n_latent)
  
  # define Gaussian process for latent factors over SNP and TF locations
  # matern 5/2 isotropic kernel
  kernel_lengthscale <- normal(5, 1, truncation = c(0, Inf))
  kernel_sd <- normal(0, 1, truncation = c(0, Inf))
  kernel <- mat52(lengthscales = c(kernel_lengthscale, kernel_lengthscale),
                  variance = kernel_sd ^ 2)
  
  # define knots for reduced-rank GP approximation
  kmn <- kmeans(coords, centers = 25)
  
  # define GPs over spatial latent factors, evaluated at all data locations
  latents_obs <- gp(x = coords,
                    kernel = kernel,
                    inducing = kmn$centers,
                    n = n_latent)
  
  # combine these with parameters to get matrices SNP frequencies at the SNP data locations
  snp_freq_logit_obs <- X_obs %*% t(parameters$beta) +
    t(parameters$loadings %*% t(latents_obs))
  snp_freq_obs <- ilogit(snp_freq_logit_obs)
  
  snp_data_index <- cbind(seq_len(nrow(snp_data)), snp_data$snp_id)
  
  # model overdispersion in the data via an overdispersion parameter rho. This
  # prior makes rho approximately uniform, but fairly nicely behaved
  # normal(0, 1.6) is similar to logit distribution; with hierarcichal prior, set
  # overall sd to: sqrt((1.6 ^ 2) - 1)
  logit_rho_mean <- normal(0, 1.3)
  logit_rho_sd <- normal(0, 1, truncation = c(0, Inf))
  logit_rho_raw <- normal(0, 1, dim = n_snp)
  logit_rho_snps <- logit_rho_mean + logit_rho_raw * logit_rho_sd
  rho_snps <- ilogit(logit_rho_snps)
  
  distribution(snp_data$snp_count) <- betabinomial_p_rho(N = snp_data$sample_size,
                                                         p = snp_freq_obs[snp_data_index],
                                                         rho = rho_snps[snp_data$snp_id])
  
  # trace only snps with data
  snp_with_data <- snp_data$snp_id %>% unique()
  
  # traced obj in model
  m <- model(kernel_lengthscale, 
             kernel_sd,
             logit_rho_mean,
             logit_rho_raw[snp_with_data],
             parameters$beta[snp_with_data,])
  
  # save all other other model obj too
  model_arrays <- list(
    parameters = parameters,
    gp = latents_obs,
    kernel_lengthscale = kernel_lengthscale,
    kernel_sd = kernel_sd,
    logit_rho_mean = logit_rho_mean,
    logit_rho_raw = logit_rho_raw,
    rho_snps = rho_snps,
    mean_snps = snp_freq_obs
  )
  
  list(
    snp_data_index = snp_data_index,
    X_obs = X_obs,
    snp_data = snp_data, 
    n_latent = n_latent,
    n_snp = n_snp,
    X_obs = X_obs,
    coords = coords,
    snp_with_data = snp_with_data,
    model_arrays = model_arrays,
    model = m
  )
  
  
}