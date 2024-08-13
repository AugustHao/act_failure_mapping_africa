. <- lapply(list.files("R", full.names = TRUE), source)
set.seed(2024-04-26)

# load in MAP data layers
PfPR <- make_synoptic_mean_rast(name = "PfPR")
Pf_incidence <- make_synoptic_mean_rast(name = "Pf_IncidenceRate")
treatment_EFT <- make_synoptic_mean_rast(name = "treatment_EFT",
                                         remove_every_second_layer = FALSE)
pop <- make_synoptic_mean_rast(name = "pop",
                               remove_every_second_layer = FALSE)

# lapply(c(PfPR,Pf_incidence,treatment_EFT), FUN =  plot)

# # coarsened version of treatment for prediction
treatment_EFT_coarse <- terra::aggregate(treatment_EFT,
                                         5,
                                         fun="modal"
                                         # think this should be treated as categorical to avoid messy country border effect?
                                         )
# plot(treatment_EFT_coarse)

# load marker prevalence data
snp_data <- readr::read_csv("data/nvs_africa.csv")

# filter to specific coordinates and select just the last 5 years
year_range <- 6
year_cutoff <- (snp_data %>% pull(year_start) %>% max) - year_range
snp_data <- snp_data %>% 
  filter(
    site_type != "Country", 
    year_start > year_cutoff)

# make table of unique snps
snp_table <- table(snp_data$snp_name) %>% 
  sort(decreasing = TRUE) %>% 
  as.data.frame() %>% 
  rename(snp = Var1,
         count = Freq) %>% 
  mutate(id = row_number())
# note we count all SNPs in n_snp, even those not modelled, to keep the id consistent
n_snp <- n_distinct(snp_table$snp)
# markers validated by WHO:
# https://www.who.int/news-room/questions-and-answers/item/artemisinin-resistance
# note some of these have associated markers, so future job is to code up marker
# by association matrix
who_markers <- c("F446I",
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

snp_table$valid <- snp_table$snp %in% who_markers

# filter data to just valid markers with at least 5 obs
valid_snps <- snp_table %>% 
  filter(valid == TRUE & count >= 5) %>% 
  pull(snp) %>% 
  as.character()
snp_data <- snp_data %>% 
  filter(snp_name %in% valid_snps)

# coordinates for sample sites
coords <- snp_data %>%
  select(longitude,latitude) %>%
  as.matrix()

# # plot raw rasters and data for vis
# ggplot() +
#   geom_spatraster(
#     data = PfPR
#   ) +
#   #facet_wrap(~lyr, nrow = 1, ncol = 3) +
#   scale_fill_gradientn(
#     #labels = scales::percent,
#     # name = "value",
#     limits = c(0, 1),
#     na.value = "transparent",
#     colours = rev(idpalette("iddu", 100)))  +
#   ggtitle(
#     label = "Plasmodium falciparum Parasite Rate (PfPR)",
#     subtitle = "Proportion of Children 2 to 10 years of age showing, on a given year, detectable Plasmodium falciparum parasite 2017-2022"
#   ) +
#   theme_snp_maps() +
#   geom_point(aes(x = longitude,
#                  y = latitude),
#              data = snp_data,
#              #col = "white",
#              shape = "+",
#              inherit.aes = TRUE,
#              size = 2)
# 
# ggsave(paste0("figures/PfPR_vis.png"),width = 10, height = 6, units = "in")
# 
# ggplot() +
#   geom_spatraster(
#     data = Pf_incidence
#   ) +
#   #facet_wrap(~lyr, nrow = 1, ncol = 3) +
#   scale_fill_gradientn(
#     #labels = scales::percent,
#     # name = "value",
#     limits = c(0, 0.7),
#     na.value = "transparent",
#     colours = colorRampPalette(c("#B9DDF1FF","#2A5783FF"))(100)
#     # note to self check if can get MAP official colours
#       )  +
#   ggtitle(
#     label = "Plasmodium falciparum Incidence Rate",
#     subtitle = "Number of newly diagnosed Plasmodium falciparum cases per 1,000 population, on a given year 2017-2022"
#   ) +
#   theme_snp_maps() +
#   geom_point(aes(x = longitude,
#                  y = latitude),
#              data = snp_data,
#              #col = "white",
#              shape = "+",
#              inherit.aes = TRUE,
#              size = 2)
# 
# ggsave(paste0("figures/Pf_incidence_vis.png"),width = 10, height = 6, units = "in")
# 
# ggplot() +
#   geom_spatraster(
#     data = treatment_EFT
#   ) +
#   #facet_wrap(~lyr, nrow = 1, ncol = 3) +
#   scale_fill_gradientn(
#     #labels = scales::percent,
#     # name = "value",
#     #limits = c(0, 1),
#     na.value = "transparent",
#     colours = colorRampPalette(c("#FFC685FF","#9E3D22FF"))(100)
#     )  +
#   ggtitle(
#     label = "Effective Treatment (EFT)",
#     subtitle = "Proportion of Malaria Cases receiving Effective Treatment with an Antimalarial Medicine 2017-2022"
#   ) +
#   theme_snp_maps() +
#   geom_point(aes(x = longitude,
#                  y = latitude),
#              data = snp_data,
#              #col = "white",
#              shape = "+",
#              inherit.aes = TRUE,
#              size = 2)
# 
# ggsave(paste0("figures/treatment_vis.png"),width = 10, height = 6, units = "in")

ggplot() +
  geom_spatraster(
    data = log(pop)#terra::aggregate(pop,20,fun = "sum", na.rm = TRUE)/1e6
  ) +
  #facet_wrap(~lyr, nrow = 1, ncol = 3) +
  scale_fill_gradientn(
    #labels = scales::percent,
    name = "log-population density",
    #limits = c(0, 1),
    na.value = "transparent",
    colours = colorRampPalette(c("#B3E0A6FF","#24693DFF"))(100)
  )  +
  ggtitle(
    label = "Human population density",
    subtitle = "natural-log transformed population count on a 5km grid averaged over 2017-2020"
  ) +
  theme_snp_maps() +
  geom_point(aes(x = longitude,
                 y = latitude),
             data = snp_data,
             #col = "white",
             shape = "+",
             inherit.aes = TRUE,
             size = 2)

ggsave(paste0("figures/pop_vis_log.png"),width = 10, height = 6, units = "in")

# extract out the design matrices (pre-scaled)

# interpolate NAs in EFT layer just for extracting data
treatment_EFT_no_na <- terra::focal(treatment_EFT,
                                    w = 21,
                                    fun = "modal",
                                    na.policy = "only")
names(treatment_EFT_no_na) <- "treatment_EFT"
plot(treatment_EFT_no_na)
points(coords)

# interpolate NAs in PfPR layer just for extracting data
PfPR_no_na <- terra::focal(PfPR,
                           w = 3,
                           fun = "mean",
                           na.policy = "only")
names(PfPR_no_na) <- "PfPR"
plot(PfPR_no_na)
points(coords)

# build covariate stack
covariates <- c(treatment_EFT_no_na,PfPR_no_na)
# note need to scale the covariate layer but keen the original scale version for visualisation too
covariates <- scale(covariates)

X_obs <- build_design_matrix(covariates = covariates, 
                             coords,
                             scale = FALSE)
# check for NAs in pred values
anyNA(X_obs)
# define number of latents
n_latent <- 4
# define model parameters
parameters <- define_greta_parameters(n_snp = n_snp,
                                      n_latent = n_latent)

# note we just use two betas (intercept + treatment EFT), because we do not want
# to fit to the two covariates specification yet

# define Gaussian process for latent factors over SNP and TF locations
# matern 5/2 isotropic kernel
kernel_lengthscale <- normal(5, 1, truncation = c(0, Inf))
kernel_sd <- normal(0, 1, truncation = c(0, Inf))
kernel <- mat52(lengthscales = c(kernel_lengthscale, kernel_lengthscale),
                variance = kernel_sd ^ 2)

# define knots for reduced-rank GP approximation
kmn <- kmeans(coords, centers = 25)

# define GPs over spatial latent factors, evaluated at all data locations
latents_obs <- gp(x = coords,
                  kernel = kernel,
                  inducing = kmn$centers,
                  n = n_latent)

# # prior predictive check
# gp_check(spatial_latents = latents_obs,
#          target_raster = treatment_EFT_coarse)

# combine these with parameters to get matrices SNP frequencies at the SNP data locations
snp_freq_logit_obs <- X_obs %*% t(parameters$beta) +
  t(parameters$loadings %*% t(latents_obs))
snp_freq_obs <- ilogit(snp_freq_logit_obs)

# define the likelihood over the SNPs
snp_data <- snp_data %>% 
  mutate(coord_id = row_number(),
         snp_id = match(snp_name, snp_table$snp))

snp_data_index <- cbind(snp_data$coord_id, snp_data$snp_id)


# model overdispersion in the data via an overdispersion parameter rho. This
# prior makes rho approximately uniform, but fairly nicely behaved
# normal(0, 1.6) is similar to logit distribution; with hierarcichal prior, set
# overall sd to: sqrt((1.6 ^ 2) - 1)
logit_rho_mean <- normal(0, 1.3)
logit_rho_sd <- normal(0, 1, truncation = c(0, Inf))
logit_rho_raw <- normal(0, 1, dim = n_snp)
logit_rho_snps <- logit_rho_mean + logit_rho_raw * logit_rho_sd
rho_snps <- ilogit(logit_rho_snps)

distribution(snp_data$snp_count) <- betabinomial_p_rho(N = snp_data$sample_size,
                                                       p = snp_freq_obs[snp_data_index],
                                                       rho = rho_snps[snp_data$snp_id])
# fit the model
m <- model(kernel_lengthscale, 
           kernel_sd,
           logit_rho_mean,
           logit_rho_raw,
           parameters$beta)
draws <- mcmc(m, 
              one_by_one = TRUE,
              warmup = 2e3,
              n_samples = 2e3)

# check convergence is not too horrible
r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)


# check overdispersion param in beta-binom
rho_snps_posterior <- calculate(rho_snps[snp_table$valid == TRUE],
                                nsim = 100, 
                                values = draws)[[1]]

rho_snps_posterior <- apply(rho_snps_posterior,2:3,mean)
rho_snps_posterior

# check loading
loadings_posterior <- calculate(parameters$loadings[snp_table$valid == TRUE,],
                                nsim = 100, 
                                values = draws)[[1]]

png("figures/loadings_posterior.png",width = 600, height = 800)
loadings_plot(loadings_posterior,
              snp_names = valid_snps)
dev.off()
# 
# loadings_posterior <- apply(loadings_posterior,2:3,mean)
# loadings_posterior

# # check betas
# beta_posterior <- calculate(parameters$beta[snp_table$valid == TRUE,],
#                             nsim = 100, 
#                             values = draws)[[1]]
# 
# beta_posterior <- apply(beta_posterior,2:3,mean)
# beta_posterior
# some negative effects, a bit weird
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
# 
# 
# 
# # check overall residual patterns
# pred_count <-  betabinomial_p_rho(N = snp_data$sample_size,
#                                   p = snp_freq_obs[snp_data_index],
#                                   rho = rho_snps[snp_data$snp_id])
# pred_count <- calculate(pred_count,
#                         nsim = 100, 
#                         values = draws)
# overall_residual <- qq_residual_diag(pred_count,
#                                      snp_data$snp_count,
#                                      title = "across snp",plot_only = FALSE)
# 
# plot(overall_residual)
# residual_pts <- data.frame(value = overall_residual$scaledResiduals,
#                            x = snp_data$longitude,
#                            y = snp_data$latitude)
# 
# lm(value ~ x:y, data = residual_pts) %>% summary
# residual_sf <- st_as_sf(residual_pts %>% 
#                           mutate(across(x:y, \(x) jitter(x, factor = 20))), 
#                         coords = c("x","y"), 
#                         crs = crs(treatment_EFT_coarse))
# plot(residual_sf)
# # no spatial relationship in overall residual pattern
# 
# # check residual for each snp
# for (i in snp_table$id[snp_table$valid == TRUE]) {
#   
#   if (snp_table$count[i] < 3) {
#     cat("skip, too few obs \n")
#   } else {
#     this_snp_data_index <- snp_data_index[snp_data_index[,2] == i,]
#     pred_count_snp <- betabinomial_p_rho(N = snp_data$sample_size[this_snp_data_index[,1]],
#                                          p = snp_freq_obs[this_snp_data_index],
#                                          rho = rho_snps[this_snp_data_index[,2]])
#     # binomial(snp_data$sample_size[this_snp_data_index[,1]],
#     #                          snp_freq_obs[this_snp_data_index])
#     pred_count_snp <- calculate(pred_count_snp,
#                                 nsim = 100, 
#                                 values = draws)
#     
#     png(paste0("figures/",snp_table$snp[i],"_residual_diag.png"),width = 800, height = 400)
#     qq_residual_diag(pred_count_snp,
#                      snp_data$snp_count[this_snp_data_index[,1]],
#                      title = snp_table$snp[i])
#     dev.off()
#   }
# }


# posterior predictive check
png(paste0("figures/latent_group_realisation.png"),width = 1000, height = 300)
gp_check(spatial_latents = latents_obs,
         target_raster = treatment_EFT_coarse,
         posterior_sims = draws)
dev.off()

# mean and var patterns look a bit odd, but maybe I'm more visually used to
# smooth rbf results

# make posterior prediction maps

# focus extent
focus_ext <- ext(25,37,-4,6)
focus_rast <- treatment_EFT %>% 
  crop(focus_ext)
coords_pixel <- terra::xyFromCell(focus_rast,
                                  cell = terra::cells(focus_rast))

# remember to sample coords from treatment layer cause it's got more holes in it
# coords_pixel <- terra::xyFromCell(treatment_EFT_coarse,
#                                   cell = terra::cells(treatment_EFT_coarse))

# get latent factors and design matrix
latents_pixel <- greta.gp::project(latents_obs, coords_pixel)
X_pixel <- build_design_matrix(covariates, coords_pixel)
anyNA(X_pixel)
# predict SNP frequencies
snp_freq_logit_pixel <- X_pixel %*% t(parameters$beta) +
  t(parameters$loadings %*% t(latents_pixel))
snp_freq_pixel <- ilogit(snp_freq_logit_pixel)

# # compute posterior samples of SNPs
post_pixel_sims <- calculate(snp_freq_pixel,
                             values = draws,
                             nsim = 1e2)

snp_freq_post_mean_pixel <- apply(post_pixel_sims$snp_freq_pixel,2:3,mean)
snp_freq_post_sd_pixel <- apply(post_pixel_sims$snp_freq_pixel,2:3,sd)
snp_freq_post_mean <- focus_rast * 0
snp_freq_post_sd <- focus_rast * 0

# plot focus
outline_vector <- gadm(
  country = c("UGA",
              "Rwanda",
              "Burundi",
              "Kenya",
              "Tanzania",
              "South Sudan",
              "Ethiopia"),
  level = 0,
  path = "~/not_synced/uga/"
  )
  
outline_vector <- crop(outline_vector,focus_rast)
outline_vector <- st_cast(st_as_sf(outline_vector),"LINESTRING")

for (i in snp_table$id[snp_table$valid == TRUE]) {
  
 snp_freq_post_mean[terra::cells(snp_freq_post_mean)] <- snp_freq_post_mean_pixel[,i]
  snp_freq_post_sd[terra::cells(snp_freq_post_sd)] <- snp_freq_post_sd_pixel[,i]
  snp_stack <- c(snp_freq_post_mean,snp_freq_post_sd)
  names(snp_stack) <- c("mean","std. dev.")
 
  plot_dat <- snp_data %>% 
    filter(snp_id == i) %>% 
    mutate(prevalence = snp_count/sample_size) %>% 
    mutate(prediction = terra::extract(snp_freq_post_mean,data.frame(longitude,latitude), ID = FALSE, raw = TRUE)) %>% 
    group_by(longitude,latitude) %>% 
    summarise(prevalence = mean(prevalence))
 
  plot_dat <- plot_dat[!is.na(
    terra::extract(focus_rast,
                   cbind(plot_dat$longitude,plot_dat$latitude))),
  ]
  
  # plot_dat <- plot_dat %>% 
  #   mutate(across(longitude:latitude, \(x) jitter(x, factor = 45))) 
  
  ggplot() +
    geom_spatraster(
      data = snp_stack$mean
    ) +
    #facet_wrap(~lyr, nrow = 1, ncol = 2) +
    scale_fill_gradientn(
      #labels = scales::percent,
      name = "mean prevalence",
      limits = c(0, 0.5),
      na.value = "transparent",
      colours = rev(idpalette("idem", 100)))  +
    ggtitle(
      label = snp_table$snp[i],
      subtitle = "prevalence in artemisinin resistance marker data"
    ) +
    theme_snp_maps() + 
    geom_spatvector(data = outline_vector,
                    linewidth = 0.2,
                    col = "grey")
  
  ggsave(paste0("figures/",snp_table$snp[i],"_pred_focus.png"),width = 10, height = 6, units = "in")
  
  ggplot() + 
    geom_spatvector(data = outline_vector,
                    linewidth = 0.2,
                    col = "grey") +
    theme_snp_maps() + 
    geom_point(aes(x = longitude,
                   y = latitude,
                   col = prevalence),
               data = plot_dat,
               #col = "white",
               shape = "+",
               inherit.aes = TRUE,
               size = 10) +
    scale_color_gradientn(name = "observed prevalence",
                          limits = c(0, 0.5),
                          na.value = "transparent",
                          colours = alpha(rev(idpalette("iddu", 100)[80:100]),1)) +
    ggtitle(
      label = snp_table$snp[i],
      subtitle = "observed prevalence in sites (averaged across samples and years)"
    )
  
  ggsave(paste0("figures/",snp_table$snp[i],"_data_focus.png"),width = 10, height = 6, units = "in")
}


select_snp <- c("C469Y","A675V","R622I","R561H")
select_id <- snp_table$id[snp_table$snp %in% select_snp]

  
  plot_dat <- snp_data %>% 
    filter(snp_id %in% select_id) %>% 
    mutate(prevalence = snp_count/sample_size) %>% 
    mutate(prediction = terra::extract(snp_freq_post_mean,data.frame(longitude,latitude), ID = FALSE, raw = TRUE)) %>% 
    group_by(longitude,latitude,snp_name) %>% 
    summarise(prevalence = mean(prevalence))


  ggplot(data = plot_dat) +
    geom_spatraster(
      data = pop/pop,
    ) +
    facet_wrap(~snp_name, nrow = 2, ncol = 2) +
    scale_fill_gradientn(
      #labels = scales::percent,
      name = NULL,
      na.value = "transparent",
      colours = grey(0.95),
      guide="none")  + 
    theme_snp_maps() +     
    geom_point(aes(x = longitude,
                   y = latitude,
                   col = prevalence),
               data = plot_dat,
               #col = "white",
               shape = 1,
               inherit.aes = TRUE,
               size = 3) +
    scale_color_gradientn(name = "observed prevalence",
                          limits = c(0, 0.5),
                          na.value = "transparent",
                          colours = alpha(rev(idpalette("iddu", 100)[80:100]),1)) +
    ggtitle(
      label = "data visualisation for selected SNPs",
      subtitle = "observed prevalence in sites (averaged across samples and years)"
    )
  ggsave(paste0("figures/selected_data_vis.png"),width = 10, height = 6, units = "in")

