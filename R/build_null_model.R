build_null_model <- function(snp_data,
                             hierarchical = TRUE,
                             n_snp = n_snp) {
  
  if (hierarchical) {
    # define model parameters
    intercept <- define_greta_parameters(n_snp = n_snp,
                                         n_latent = 1)$beta[,1]
  } else {
    intercept <- greta::normal(0,3,dim = n_snp)
  }

  snp_freq_obs <- ilogit(intercept)
  
  snp_data_index <- snp_data$snp_id
  
  if (hierarchical) {
    logit_rho_mean <- normal(0, 1.3)
    logit_rho_sd <- normal(0, 1, truncation = c(0, Inf))
    logit_rho_raw <- normal(0, 1, dim = n_snp)
    logit_rho_snps <- logit_rho_mean + logit_rho_raw * logit_rho_sd
    rho_snps <- ilogit(logit_rho_snps)
  } else {
    rho_snps <- ilogit(normal(0, 1.6,dim = n_snp))
  }
  
  distribution(snp_data$snp_count) <- betabinomial_p_rho(N = snp_data$sample_size,
                                                         p = snp_freq_obs[snp_data_index],
                                                         rho = rho_snps[snp_data$snp_id])
  
  # traced obj in model
  m <- model(intercept,rho_snps)
  
  # save all other other model obj too
  model_arrays <- list(
    rho_snps = rho_snps,
    mean_snps = snp_freq_obs
  )
  
  list(
    model_arrays = model_arrays,
    model = m
  )
  
  
}