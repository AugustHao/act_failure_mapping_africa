#' make marginal response plots wrt a chosen covariate, parameter and draws
#' variable names, and model specification are hard coded, headache to fix I
#' think but can deal with later
#'
#' @param target_covariate
#' @param covariates
#'
#' @return
#' @export
#'
#' @examples
marginal_response_plot <- function(target_covariate,
                                   snp_idx,
                                   covariates,
                                   transformations = NULL,
                                   draws = draws,
                                   response_data = NULL,
                                   response = "snp_count",
                                   discrete_covariate = FALSE) { 
  # marginal response plots
  x <- marginal_design_matrix(covariates,
                              target_covariate = target_covariate,
                              transformations = transformations)
  
  marginal_pred <- x %*% t(parameters$beta[snp_idx,]) 
  marginal_pred <- ilogit(marginal_pred)
  
  marginal_pred <-  betabinomial_p_rho(N = 100,
                                    p = marginal_pred,
                                    rho = rho_snps[snp_idx])
  
  marginal_pred <- calculate(marginal_pred,
                             values = draws,
                             nsim = 1e2)[[1]]
  
  #calculate CI values for ribbon plot
  ci_90_lo <- apply(marginal_pred, 2:3, quantile, c(0.05),na.rm = TRUE)
  ci_90_hi <- apply(marginal_pred, 2:3, quantile, c(0.95),na.rm = TRUE)
  ci_50_hi <- apply(marginal_pred, 2:3, quantile, c(0.75),na.rm = TRUE)
  ci_50_lo <- apply(marginal_pred, 2:3, quantile, c(0.25),na.rm = TRUE)
  
  vals <- data.frame(ci_90_lo,ci_90_hi,ci_50_hi,ci_50_lo, x = x[target_covariate])
  # hard code name, probably not necessary
  names(vals)[5] <- "x"
  ribbon_colour <- idpalette::idem(1)
  #base_colour <- grey(0.3)
  
  p <- ggplot(vals, aes(x =x)) + 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_90_lo,
                                      ymax = ci_90_hi),
                         fill = ribbon_colour,
                         alpha = 0.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_50_lo,
                                      ymax = ci_50_hi),
                         fill = ribbon_colour,
                         alpha = 0.5) +
    # ggplot2::geom_line(ggplot2::aes(y = ci_90_lo),
    #                    colour = base_colour,
    #                    alpha = 0.8) +
    # ggplot2::geom_line(ggplot2::aes(y = ci_90_hi),
    #                    colour = base_colour,
    #                    alpha = 0.8) +
    scale_y_continuous(name = "marginal response") + 
    cowplot::theme_cowplot() +
    cowplot::panel_border(remove = TRUE) +
    ggplot2::theme(legend.position = "none",
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, face = "bold"),
                   axis.title.y.right = ggplot2::element_text(vjust = 0.5, angle = 90),
                   panel.spacing = ggplot2::unit(1.2, "lines")
    ) +
    ggtitle(
      label = "marginal response with respect to one covariate",
      subtitle = "light and dark ribbons represent 90% and 50% credible intervals based on 1000 posterior simulations")
  
  if (discrete_covariate) { 
    p <- p + 
      scale_x_discrete(name = target_covariate) 
  } else {
    p <- p + 
      scale_x_continuous(name = target_covariate,n.breaks = 10) 
  }
  
  if (!is.null(response_data)) {
    
    response_data <- response_data %>% 
      select(c(response,target_covariate)) %>% 
      rename(x = target_covariate,
             y= response)
    p <- p + 
      geom_point(aes(x = x, 
                     y = y), 
                 data = response_data,
                 col = alpha(idpalette::iddu(5)[4],0.4),
                 shape = "+", 
                 inherit.aes = FALSE,
                 size = 3) 
    
  }
  # base plot code for debug purpose
  # 
  # plot(test$temp_mean,marginal_pred[1,,],type = "n",xlab = "temp", ylab = "response")
  # 
  # for (i in 1:10) {
  #     points(test$temp_mean,marginal_pred[i,,],type = "l")
  # }
  p
  
}