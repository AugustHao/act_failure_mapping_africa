plot_marginal_response <- function(X_obs,
                                   snp_data,
                                   snp_table,
                                   snp_data_index,
                                   draws,
                                   beta = model$model_arrays$parameters$beta,
                                   rho_snps = model$model_arrays$rho_snps) {
  
  
  X_tib <- as_tibble(X_obs)
  colnames(X_tib)[2] <- "treatment_EFT"
  covars_names <- colnames(X_tib)[-1]
  # covars_names <- c(covars_names,"landuse")
  response_data <- cbind(snp_data,X_tib)
  
  # make a prev out of 100 for plotting
  response_data$prev_100 <- response_data$snp_count/response_data$sample_size * 100
  # check residual for each snp
  for (i in snp_table$id[snp_table$valid == TRUE]) {
    
    this_snp_data_index <- snp_data_index[snp_data_index[,2] == i,]
    response_data_this <- response_data[this_snp_data_index[,1],]
    
    for (j in 1:length(covars_names)) {
      p <- marginal_response_plot(covars_names[j],
                                  covariates = X_tib,
                                  beta = beta,
                                  rho_snps = rho_snps,
                                  snp_idx = i,
                                  transformations = NULL,
                                  draws = draws,
                                  response_data = response_data_this,
                                  response_is_rug = TRUE,
                                  response = "prev_100")
      ggsave(paste0("figures/",
                    snp_table$snp[i],
                    "_to_",
                    covars_names[j],
                    "_marginal_response.png"),
             width = 10,
             height = 5,
             units = "in",
             plot = p)
    }
    
  }
}