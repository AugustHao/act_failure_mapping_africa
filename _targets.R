# Load packages required to define the pipeline:
library(targets)
library(geotargets)

# Set target options:
tar_option_set(
  packages = c("tidyverse",
               "terra",
               "tidyterra",
               "greta",
               "greta.gp",
               "geodata",
               "DHARMa",
               "mapview",
               "sf",
               "patchwork",
               "idpalette",
               "corrplot",
               "blockCV"),
  seed = 2024-04-26
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Replace the target list below with your own:
list(
  # load raster of all covariates and vis layers
  tar_terra_rast(
    name = PfPR_rast,
    command = make_synoptic_mean_rast(name = "PfPR")
  ),
  tar_terra_rast(
    name = Pf_incidence_rast,
    command = make_synoptic_mean_rast(name = "Pf_IncidenceRate")
  ),
  tar_terra_rast(
    name = treatment_EFT_rast,
    command = make_synoptic_mean_rast(name = "treatment_EFT",
                                      remove_every_second_layer = FALSE)
  ),
  tar_terra_rast(
    name = pop_rast,
    command = make_synoptic_mean_rast(name = "pop",
                                      remove_every_second_layer = FALSE)
  ),
  
  # coarsen treatment EFT layer to deal with holes in sampled values
  tar_terra_rast(
    treatment_EFT_coarse,
    terra::aggregate(treatment_EFT_rast,
                     5,
                     fun="modal"
                     # treat as discrete so use modal
    )
  ),
  
  # load response data
  tar_target(
    snp_data_raw,
    readr::read_csv("data/nvs_africa.csv")
  ),
  
  # make a table of all snps
  tar_target(
    snp_table,
    table(snp_data_raw$snp_name) %>% 
      sort(decreasing = TRUE) %>% 
      as.data.frame() %>% 
      rename(snp = Var1,
             count = Freq) %>% 
      mutate(id = row_number(),
             valid = snp %in% who_markers())
  ),
  
  # count snps
  tar_target(
    n_snp,
    n_distinct(snp_table$snp)
  ),
  
  # filter snps to WHO valid markers and at least 5 obs
  tar_target(
    valid_snps,
    snp_table %>% 
      filter(valid == TRUE & count >= 5) %>% 
      pull(snp) %>% 
      as.character()
  ),
  
  # filter rows in the dataset
  tar_target(
    snp_data,
    filter_snp_data(snp_data_raw, valid_snps = valid_snps)
  ),
  
  # get the snp coords
  tar_target(
    snp_coords,
    snp_data %>%
      select(longitude,latitude) %>%
      as.matrix()
  ),
  
  # plot raster vis
  tar_target(
    save_PfPR_rast_vis,
    plot_synoptic_rast(raster = PfPR_rast,
                       plot_name = "Plasmodium falciparum Parasite Rate (PfPR)", 
                       plot_desc = "Proportion of Children 2 to 10 years of age showing, on a given year, detectable Plasmodium falciparum parasite 2017-2022",
                       colours = default_colours(),
                       overlay_sites = TRUE, 
                       write_to_disk = TRUE,
                       limits = c(0,1),
                       filename = "PfPR_vis",
                       site_data = snp_data
                       )
  ),
  
  tar_target(
    save_incidence_rast_vis,
    plot_synoptic_rast(raster = Pf_incidence_rast,
                       plot_name = "Plasmodium falciparum Incidence Rate", 
                       plot_desc = "Number of newly diagnosed Plasmodium falciparum cases per 1,000 population, on a given year 2017-2022",
                       colours = colorRampPalette(c("#B9DDF1FF","#2A5783FF"))(100),
                       overlay_sites = TRUE, 
                       write_to_disk = TRUE,
                       limits = c(0,0.7),
                       filename = "Pf_incidence_vis",
                       site_data = snp_data
    )
  ),
  
  tar_target(
    save_treatment_EFT_rast_vis,
    plot_synoptic_rast(raster = treatment_EFT_rast,
                       plot_name = "Effective Treatment (EFT)", 
                       plot_desc = "Proportion of Malaria Cases receiving Effective Treatment with an Antimalarial Medicine 2017-2022",
                       colours = colorRampPalette(c("#FFC685FF","#9E3D22FF"))(100),
                       overlay_sites = TRUE, 
                       write_to_disk = TRUE,
                       filename = "treatment",
                       site_data = snp_data
    )
  ),
  
  tar_target(
    save_pop_rast_vis,
    plot_synoptic_rast(raster = log(pop_rast),
                       plot_name = "Human population density", 
                       plot_desc = "natural-log transformed population count on a 5km grid averaged over 2017-2020",
                       colours = colorRampPalette(c("#B3E0A6FF","#24693DFF"))(100),
                       overlay_sites = TRUE, 
                       write_to_disk = TRUE,
                       filename = "log_pop",
                       site_data = snp_data
    )
  ),
  
  # prepare covariate rast stack
  tar_terra_rast(
    covariate_rast,
    build_covariate_stack(treatment_EFT_rast,PfPR_rast)
  ),
  
  # get design matrix
  tar_target(
    X_obs, 
    build_design_matrix(covariates = covariate_rast, 
                        coords = snp_coords,
                        scale = FALSE)
  ),
  
  # get spatial CV blocks
  tar_target(
    cv_blocks,
    get_cv_blocks(snp_data,covariate_rast)
  ),
  
  # incorporate indices in the snp data
  tar_target(snp_data_idx,
             snp_data %>% 
               mutate(coord_id = row_number(),
                      snp_id = match(snp_name, snp_table$snp),
                      block_id = cv_blocks$folds_ids)
  ),
  
  # cross validate with null model
  tar_target(
    cv_null_model_result,
    cross_validate(snp_data_idx, 
                   n_latent = 4,
                   n_snp = n_snp,
                   X_obs = X_obs,
                   coords = snp_coords,
                   null_model = TRUE)
  ),
  
  # cross validate
  tar_target(
    cv_result,
    cross_validate(snp_data_idx, 
                   n_latent = 4,
                   n_snp = n_snp,
                   X_obs = X_obs,
                   coords = snp_coords)
  ),
  
  # define model
  tar_target(
    model,
    build_snp_model(snp_data = snp_data_idx, 
                    n_latent = 4,
                    n_snp = n_snp,
                    X_obs = X_obs,
                    coords = snp_coords)
  ),
  
  # diagnostic plot for gp
  tar_target(
    gp_prior_check_plot,
    gp_check(model$model_arrays$gp,
             treatment_EFT_coarse)
  ),
  
  # fit model
  tar_target(
    draws,
    mcmc(model$model, 
         one_by_one = TRUE,
         warmup = 1e3,
         n_samples = 1e3)
  ),
  
  # check convergence is not too horrible
  tar_target(
    r_hats,
    summary(coda::gelman.diag(draws,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf)
  ),
  
  tar_target(
    trace_plot,
    bayesplot::mcmc_trace(draws)
  ),
  
  tar_target(
    save_trace_plot,
    ggsave("figures/mcmc_trace.png",
           plot = trace_plot,
           width = 10,
           height = 7,
           bg = "white")
  )
)
