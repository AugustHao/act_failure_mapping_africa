cross_validate <- function(snp_data_idx,
                           n_latent = n_latent,
                           n_snp = n_snp,
                           X_obs = X_obs,
                           coords = snp_coords,
                           null_model = FALSE) {
  
  posterior_deviance_result <- c()
  deviance_posterior_result <- c()
  
  for (fold in unique(snp_data_idx$block_id)) {
    train <- snp_data_idx$block_id != fold
    test <- snp_data_idx$block_id == fold
    
    if (null_model) {
      train_model <- build_null_model(snp_data_idx[train,],
                                      n_snp = n_snp)
    } else {
      train_model <- build_snp_model(snp_data = snp_data_idx[train,], 
                                     n_latent = n_latent,
                                     n_snp = n_snp,
                                     X_obs = X_obs[train,],
                                     coords = coords[train,])
    }
    
    draws <- mcmc(train_model$model, 
                  one_by_one = !null_model, # only need this for full model
                  warmup = 1e3,
                  n_samples = 1e3)
    
    # predict
    
    if (null_model) {
      snp_freq_test <- X_obs[test,1] %*% t(train_model$model_arrays$mean_snps)
    } else {
      test_gp <- greta.gp::project(train_model$model_arrays$gp,
                                   x_new = coords[test,])
      
      snp_freq_logit_test <- X_obs[test,] %*% t(train_model$model_arrays$parameters$beta) +
        t(train_model$model_arrays$parameters$loadings %*% t(test_gp))
      snp_freq_test <- ilogit(snp_freq_logit_test)
    }
    
    pred_rho <- train_model$model_arrays$rho_snps
    
    p_test <- snp_freq_test[cbind(1:nrow(snp_freq_test),
                                  pull(snp_data_idx[test,"snp_id"]))]
    rho_test <- pred_rho[pull(snp_data_idx[test,"snp_id"])]
    
    posterior_deviance <- betabinomial_deviance_p_rho(pull(snp_data_idx[test,"snp_count"]),
                                                      pull(snp_data_idx[test,"sample_size"]),
                                                      p = p_test,
                                                      rho = rho_test
                                                      )
    
    posterior_deviance_mean <- calculate(posterior_deviance,values = draws,nsim = 100)[[1]] %>% mean
    
    p_test_post_mean <- calculate(p_test,values = draws,nsim = 100)[[1]] %>% mean
    
    p_test_post_mean <- calculate(p_test,values = draws,nsim = 100)[[1]] %>% mean
    
    deviance_posterior_mean <- betabinomial_deviance_p_rho(pull(snp_data_idx[test,"snp_count"]),
                                                           pull(snp_data_idx[test,"sample_size"]),
                                                           p = p_test_post_mean,
                                                           rho = p_test_post_mean
    ) %>% mean
    
    posterior_deviance_result <- c(posterior_deviance_result,posterior_deviance_mean)
    
    deviance_posterior_result <- c(deviance_posterior_result,deviance_posterior_mean)
  }

  list(
    posterior_deviance = posterior_deviance_result,
    deviance_posterior = deviance_posterior_result
  )
}