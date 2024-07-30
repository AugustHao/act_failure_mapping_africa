# library(readxl)
# library(tidyverse)
# library(geodata)
# 
# nvs_dat <- read_xlsx("data/Final-NVS-database-February-2024_plus 2 TZA papers-02Jul2024.xlsx",
#                      sheet = 1)
# 
# # some notes for later stuff with WWARN
# table(nvs_dat$`Information level`,nvs_dat$`present in WWARN`)
# # seems to indicate a fair number of WWARN data are not specific - check with
# # WWARN assumption with other work
# 
# nvs_k13_africa <- nvs_dat %>%
#   # filter to K13
#   filter(gene %in% c("pfkelch13","Pfkelch13")) %>%
#   # filter to Africa
#   filter(continent == "Africa")
# 
# table(nvs_k13_africa$country) %>% sort(decreasing = TRUE)
# # # pull out Uganda first for a test run
# # nvs_k13_uganda <- nvs_k13_africa %>%
# #   filter(country == "Uganda")
# 
# # # get extent map for our basic analyses - using outline of Uganda
# # uganda_vector <- gadm(
# #   country = "UGA",
# #   level = 0,
# #   path = "~/not_synced/uga/"
# # )
# # 
# # plot(uganda_vector)
# # points(nvs_k13_uganda$longitude,nvs_k13_uganda$latitude,pch = 19)
# # checked range of long lat, looks good
# 
# # check if need to match to admin units
# table(nvs_k13_africa$`Information level`)
# # yes, meaning need to retain site name to match to polygons later
# 
# # check SNPs
# table(nvs_k13_africa$`SNP for mapping`) %>% sort(decreasing = TRUE)
# table(nvs_k13_africa$`SNP reported`) %>% sort(decreasing = TRUE)
# 
# # `SNP for mapping` aggregates across those ambiguous markers, use this one
# 
# # extract set of columns needed for analysis
# nvs_k13_africa <- nvs_k13_africa %>%
#   select(year_start = `start of data collection`,
#          year_end = `end of data collection`,
#          longitude,
#          latitude,
#          site = Site,
#          site_type = `Information level`,
#          sample_size = tested,
#          snp_name = `SNP for mapping`,
#          snp_count = `present incl mix`,
#          title, # for de-dup
#          author = `first author`# for de-dup
#          )
# 
# # set correct column type
# nvs_k13_africa <- nvs_k13_africa %>%
#   mutate(across(year_start:latitude, as.numeric)) %>%
#   mutate(across(c(sample_size,snp_count), as.numeric))
# 
# table(nvs_k13_africa$year_start,nvs_k13_africa$year_end)
# 
# # check for duplication
# nvs_k13_africa %>%
#   group_by(snp_name,year_start,year_end,longitude,latitude,sample_size,snp_count,title,author) %>%
#   summarise(count = n()) %>% filter(count > 1) %>% nrow()
# 
# # duplications result from non-unique mapping between SNP recorded and SNP for
# # mapping, nothing to worry, and only one duplicated entry has non-zero
# # prevalence (checked in detail, indeed resulting from the above inconsistency).
# # Safe to remove duplcates
# nvs_k13_africa <- nvs_k13_africa[!duplicated(nvs_k13_africa),]
# # check again
# nvs_k13_africa %>%
#   group_by(snp_name,year_start,year_end,longitude,latitude,sample_size,snp_count,title,author) %>%
#   summarise(count = n()) %>% filter(count > 1) %>% nrow()
# 
# # drop biblio columns fo dedup
# nvs_k13_africa <- nvs_k13_africa %>%
#   select(-c(title,author))
# 
# # save uganda data
# write.csv(nvs_k13_africa,"data/nvs_africa.csv",row.names = FALSE)
