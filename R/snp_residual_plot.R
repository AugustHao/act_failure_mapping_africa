snp_residual_plot <- function(snp_table,
                              snp_data_index,
                              snp_data,
                              snp_freq_obs,
                              rho_snps,
                              draws = draws
                              ) {
  
  # check residual for each snp
  for (i in snp_table$id[snp_table$valid == TRUE]) {
    
    if (snp_table$count[i] < 3) {
      cat("skip, too few obs \n")
    } else {
      this_snp_data_index <- snp_data_index[snp_data_index[,2] == i,]
      pred_count_snp <- betabinomial_p_rho(N = snp_data$sample_size[this_snp_data_index[,1]],
                                           p = snp_freq_obs[this_snp_data_index],
                                           rho = rho_snps[this_snp_data_index[,2]])

      pred_count_snp <- calculate(pred_count_snp,
                                  nsim = 100,
                                  values = draws)
      
      png(paste0("figures/",snp_table$snp[i],"_residual_diag.png"),width = 800, height = 400)
      qq_residual_diag(pred_count_snp,
                       snp_data$snp_count[this_snp_data_index[,1]],
                       title = snp_table$snp[i])
      dev.off()
    }
  }
}
