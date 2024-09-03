plot_loading <- function(loadings_posterior,
                         snp_names,
                         save_local = TRUE) {
  
  #calculate CI values for coloured grids
  # ci_90_lo <- apply(loadings_posterior, 2:3, quantile, c(0.05),na.rm = TRUE)
  # ci_90_hi <- apply(loadings_posterior, 2:3, quantile, c(0.95),na.rm = TRUE)
  ci_mean <- apply(loadings_posterior, 2:3, mean,na.rm = TRUE)
  
  colnames(ci_mean) <- paste0("latent_", seq_len(ncol(ci_mean)))
  rownames(ci_mean) <- snp_names
  
  if (save_local) {
    png("figures/loadings_posterior.png",width = 600, height = 800)
    corrplot::corrplot(ci_mean,
                       mar = c(0, 0, 2, 0),
                       title = "mean SNP to spatial latent factor loading",
                       is.corr = FALSE,
                       method="color",
                       #col = idpalette::idem(100),
                       tl.col="black", 
                       tl.srt=45,
                       cl.length = 10,
                       cl.align.text = 'l')
    dev.off()
  } else {
    corrplot::corrplot(ci_mean,
                       mar = c(0, 0, 2, 0),
                       title = "mean SNP to spatial latent factor loading",
                       is.corr = FALSE,
                       method="color",
                       #col = idpalette::idem(100),
                       tl.col="black", 
                       tl.srt=45,
                       cl.length = 10,
                       cl.align.text = 'l')
  }
}


