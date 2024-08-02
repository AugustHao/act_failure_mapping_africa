loadings_plot <- function(loadings_posterior,
                          snp_names) {


  #calculate CI values for coloured grids
  ci_90_lo <- apply(loadings_posterior, 2:3, quantile, c(0.05),na.rm = TRUE)
  ci_90_hi <- apply(loadings_posterior, 2:3, quantile, c(0.95),na.rm = TRUE)
  ci_mean <- apply(loadings_posterior, 2:3, mean,na.rm = TRUE)
  # ci_50_hi <- apply(loadings_posterior, 2:3, quantile, c(0.75),na.rm = TRUE)
  # ci_50_lo <- apply(loadings_posterior, 2:3, quantile, c(0.25),na.rm = TRUE)
  
  # colnames(ci_90_hi) <- paste0("latent_", seq_len(ncol(ci_90_lo)))
  # rownames(ci_90_hi) <- snp_names
  # ci_90_hi_plot <- corrplot(ci_90_hi,
  #                           title = "95% quantile of SNP to spatial latent factor loading posterior",
  #                           is.corr = FALSE,
  #                           method="color",
  #                           #col = idpalette::idem(100),
  #                           tl.col="black", 
  #                           tl.srt=45)
  # 
  # colnames(ci_90_lo) <- paste0("latent_", seq_len(ncol(ci_90_lo)))
  # rownames(ci_90_lo) <- snp_names
  # ci_90_lo_plot <- corrplot(ci_90_lo,
  #                           title = "5% quantile of SNP to spatial latent factor loading posterior",
  #                           is.corr = FALSE,
  #                           method="color",
  #                           #col = idpalette::idem(100),
  #                           tl.col="black", 
  #                           tl.srt=45)
  
  colnames(ci_mean) <- paste0("latent_", seq_len(ncol(ci_mean)))
  rownames(ci_mean) <- snp_names
  corrplot(ci_mean,
           mar = c(0, 0, 2, 0),
           title = "mean SNP to spatial latent factor loading posterior",
           is.corr = FALSE,
           method="color",
           #col = idpalette::idem(100),
           tl.col="black", 
           tl.srt=45,
           cl.length = 10,
           cl.align.text = 'l')

  }


