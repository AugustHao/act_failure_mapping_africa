build_null_model <- function(snp_data,
                             n_snp = n_snp) {
  

  intercept <- greta::normal(0,3,dim = n_snp)
  snp_freq_obs <- ilogit(intercept)
  
  snp_data_index <- snp_data$snp_id
  
  rho_snps <- ilogit(normal(0, 1.6,dim = n_snp))
  
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