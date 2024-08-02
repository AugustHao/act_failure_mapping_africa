#' takes a covariate data frame and a name arg for the target covariate, make a
#' sensible marginal response design matrix to predict to for marginal response
#' plots
#'
#' @param covariates
#' @param target_covariate
#'
#' @return
#' @export
#'
#' @examples
marginal_design_matrix <- function(covariates,
                                   target_covariate,
                                   transformations = NULL) {
  
  # extract range to plot
  target_range <- range(covariates %>% dplyr::pull(target_covariate))
  
  # target column
  target_x <- seq(target_range[1],
                  target_range[2],
                  length.out = 1e3)
  
  # if transformations are applied, match to the respective column names
  if (!is.null(transformations)) {
    
    # apply transformations, dunno how to do this elegantly other than coding
    # case by case, so just quadratic implementation now
    if (transformations == "^2") {
      target_x_transform <- target_x^2
    }
    
    target_covariate <- c(target_covariate,paste0(target_covariate,transformations))
  }
  
  # get the non-target column names to replicate
  non_target_names <- names(covariates)[!(names(covariates) %in% target_covariate)]
  n_non_targets <- length(non_target_names)
  
  # make them an empty matrix
  non_target_mat <- matrix(0,nrow = 1e3, ncol = n_non_targets)
  
  # make design matrix
  if (!is.null(transformations)) {
    x_df <- data.frame(cbind(target_x,
                             target_x_transform,
                             non_target_mat))
    
  } else {
    x_df <- data.frame(cbind(target_x,non_target_mat))
  }
  
  # fix names
  names(x_df) <- c(target_covariate,non_target_names)
  
  # finally, get in the same col order as input matrix
  x_df<-x_df[names(covariates)]
  
  x_df
}