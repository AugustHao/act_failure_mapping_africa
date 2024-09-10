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
  
  # expand to fill in 0 counts
  wwarn_africa <- wwarn_africa %>% 
    group_by(longitude,latitude,year_start,sample_size) %>% 
    pivot_wider(values_from = snp_count,names_from = snp_name, 
                values_fn = function(x)unique(x)[1],
                values_fill = 0) %>% 
    pivot_longer(cols = !c("year_start","longitude","latitude","sample_size","site_type","site"),
                 names_to = "snp_name",
                 values_to = "snp_count") %>% 
    ungroup()
  wwarn_africa
}