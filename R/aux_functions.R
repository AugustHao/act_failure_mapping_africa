# plotting themes for snp prediction maps
theme_snp_maps <- function() {
  theme_minimal() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_rect(fil = "white", linetype = "blank"))
}


# given a vector of sample sizes N, predicted proportions p, and overdispersion
# parameter rho (on unit interval, where 0 means no overdispersioan and 1 means
# maximum), returna betatbinomial-distributed variable greta array
betabinomial_p_rho <- function(N, p, rho) {
  
  # model the observation (betabinomial) sd as a multiplier on the binomial sd,
  # accounting for additional error due to nonindependent sampling of individuals
  # from the population. This is based on the INLA parameterisation
  
  # solve for a and b:
  #   p = a / (a + b)
  #   rho = 1 / (a + b + 1)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  
  # define betabinomial according to the greta interface
  beta_binomial(size = N, alpha = a, beta = b)
  
}

# markers validated by WHO:
# https://www.who.int/news-room/questions-and-answers/item/artemisinin-resistance
# note some of these have associated markers, so future job is to code up marker
# by association matrix
who_markers <- function(){
  c("F446I",
    "N458Y",
    "C469Y",
    "M476I",
    "Y493H",
    "R539T",
    "I543T",
    "P553L",
    "R561H",
    "P574L",
    "C580Y",
    "R622I",
    "A675V")
}

# filter snp data to exclude country level coordinates, those too old and those
# outside the valid snp table
filter_snp_data <- function(snp_data, 
                            year_range = 6, 
                            valid_snps,
                            filter_rast = FALSE,
                            rast = NULL) {
  year_cutoff <- (snp_data %>% pull(year_start) %>% max) - year_range
  snp_data <- snp_data %>% 
    filter(
      site_type != "Country", 
      year_start > year_cutoff,
      snp_name %in% valid_snps) 
  
  if (filter_rast) {
    coords <- snp_data %>%
      select(longitude,latitude) %>%
      as.matrix()
    cell_ids <- terra::cellFromXY(rast, coords)
    vals <- terra::extract(rast, cell_ids)
    valid_ids <- complete.cases(vals)
    snp_data <- snp_data[valid_ids,]
  }
  snp_data
}

# build covariate stack, with interpolation over NA and scaling
build_covariate_stack <- function(treatment_EFT_rast, PfPR_rast) {
  # interpolate NAs in EFT layer just for extracting data
  treatment_EFT_no_na <- terra::focal(treatment_EFT_rast,
                                      w = 21,
                                      fun = "modal",
                                      na.policy = "only")
  names(treatment_EFT_no_na) <- "treatment_EFT"
  
  # interpolate NAs in PfPR layer just for extracting data
  PfPR_no_na <- terra::focal(PfPR_rast,
                             w = 3,
                             fun = "mean",
                             na.policy = "only")
  names(PfPR_no_na) <- "PfPR"
  
  # build covariate stack
  covariates <- c(treatment_EFT_no_na,PfPR_no_na)
  # note need to scale the covariate layer but keen the original scale version for visualisation too
  covariates <- scale(covariates)
  covariates
}

get_cv_blocks <- function(snp_data,covariate_rast) {
  # turn data into sf for blockCV
  snp_data_sf <- st_as_sf(snp_data,
                          coords = c("longitude","latitude"),
                          crs = crs(covariate_rast))
  # get cv blocks
  cv_blocks <- blockCV::cv_spatial(snp_data_sf,
                                   r = covariate_rast,
                                   k = 5)
  cv_blocks
}

# # hierarchical rho prior for beta binom distribution
# rho_prior <- function(n_snp) { 
#   # model overdispersion in the data via an overdispersion parameter rho. This
#   # prior makes rho approximately uniform, but fairly nicely behaved
#   # normal(0, 1.6) is similar to logit distribution; with hierarcichal prior, set
#   # overall sd to: sqrt((1.6 ^ 2) - 1)
#   logit_rho_mean <- normal(0, 1.3)
#   logit_rho_sd <- normal(0, 1, truncation = c(0, Inf))
#   logit_rho_raw <- normal(0, 1, dim = n_snp)
#   logit_rho_snps <- logit_rho_mean + logit_rho_raw * logit_rho_sd
#   rho_snps <- ilogit(logit_rho_snps)
#   
#   rho_snps
# }

# # spatial gp prior
# gp_objs <- function(coords, n_latent) {
#   # define Gaussian process for latent factors over SNP locations
#   # matern 5/2 isotropic kernel
#   kernel_lengthscale <- normal(5, 1, truncation = c(0, Inf))
#   kernel_sd <- normal(0, 1, truncation = c(0, Inf))
#   kernel <- mat52(lengthscales = c(kernel_lengthscale, kernel_lengthscale),
#                   variance = kernel_sd ^ 2)
#   
#   # define knots for reduced-rank GP approximation
#   kmn <- kmeans(coords, centers = 50)
#   
#   # define GPs over spatial latent factors, evaluated at all data locations
#   latents_obs <- gp(x = coords,
#                     kernel = kernel,
#                     inducing = kmn$centers,
#                     n = n_latent)
#   
#   list(k_l = kernel_lengthscale,
#        k_sd = kernel_sd,
#        latents_obs = latents_obs)
# }

# compute observation-level log like for beta binomial with p rho parameterisation
betabinomial_ll_p_rho <- function(y, N, p, rho) {
  
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  # get ll
  # https://cran.r-project.org/web/packages/extras/vignettes/beta-binomial-deviance-residuals.html
  lgamma(N + 1) - lgamma(y + 1) - lgamma(N - y + 1) + lgamma(y + a) + lgamma(N - y + b) - lgamma(N + a + b) + lgamma(a + b) - lgamma(a) - lgamma(b) 
}

# fudged deviance for beta binom, saturated model does not have a ll of 0 but
# for purpose of model comparison it's enough
betabinomial_deviance_p_rho <- function(y, N, p, rho) {
  -2 * betabinomial_ll_p_rho(y, N, p, rho)
}