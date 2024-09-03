read_wwarn <- function(file) {
  
  file <- readxl::read_excel(file)
  
  # subset to africa
  wwarn_africa <- file %>% 
    filter(continent == "Africa")
  
  # process
  wwarn_africa <- wwarn_africa %>% 
    ## filter to exact locs only
    #filter(estLoc == 0) %>%  
    # select to match NVS structure 
    select(year_start = year,
           longitude = lon,
           latitude = lat,
           sample_size = tested,
           snp_count = present,
           snp_name = mutation,
           site = site,
           site_type = estLoc) %>% 
    # force type correction
    mutate(across(year_start:snp_count, as.numeric))
  
  wwarn_africa
}