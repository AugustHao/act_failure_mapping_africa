#' makes a DHARMa style scaled residual QQ plot to check for
#' over/underdispersion and other model diagnostic needs
#'
#' @param count_posterior
#' @param count_true
#'
#' @return
#' @export
#'
#' @examples
qq_residual_diag <- function(count_posterior,
                             count_true,
                             title = NULL,
                             plot_only = TRUE){
  
  count_posterior <- t(count_posterior[[1]][,,1])
  
  count_posterior_median <- apply(count_posterior,1,median)
  
  DHARMa_obj <- DHARMa::createDHARMa(simulatedResponse = count_posterior,
                                     observedResponse = count_true,
                                     fittedPredictedResponse = count_posterior_median,
                                     integerResponse = TRUE)
  
  if (plot_only) {
    DHARMa:::plot.DHARMa(DHARMa_obj, title = paste(title,"DHARMa residual",sep = " "))
  } else {
    return(DHARMa_obj)
  }

  
}