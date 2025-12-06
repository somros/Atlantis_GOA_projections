# 12/05/2025
# Alberto Rovellini
# Create ex-vessel price estimates
# multiply catch in projection by these values
# exploratory analysis

# ex-vessel data based on https://www.federalregister.gov/documents/2020/12/18/2020-27824/fisheries-of-the-exclusive-economic-zone-off-alaska-north-pacific-observer-program-standard
# 2020 estimates based on landings and volumes from 2017, 2018, 2019
# Using values for area "GOA"

# mapping species to Atlantis groups
# weighting by proportions for multispecies functional groups based on Rovellini et al. (2024), which was based on BTS data
# The most difficult step is identifying what gear we should use
# The model has no fleets for now, only one F. 
# For most of the rockfish there are substantial differences between gear types
# for now just use the max, but this will need to be revisited
# The right approach would be weighting the prices based on proportion of total mortality from each gear
# Can probably get that from catch data / the metier analysis

library(readxl)
library(tidyverse)

price_dat <- read_xlsx("data/ex_vessel_prices_2020.xlsx", sheet = 1, range = "A1:E43")
key <- read_xlsx("data/ex_vessel_prices_2020.xlsx", sheet = 3, range = "A1:C40")
key <- key %>% filter(Prop > 0)

# select highest price for each species
# TODO: weigh this by mortality per gear

price_dat_2 <- price_dat %>%
  pivot_longer(-Species, names_to = "gear", values_to = "price") %>%
  filter(!is.na(price)) %>%
  group_by(Species) %>%
  slice_max(price) %>%
  ungroup() %>%
  select(-gear)

# aggregate by grp and weight
# some groups have no price data (SCU, DFD, etc.). Price those at 0
price_dat_3 <- price_dat_2 %>%
  left_join(key, by = "Species") %>%
  filter(!is.na(Code)) %>% # this drops codes that were 0's in the biomass allocation from BTS
  group_by(Code) %>%
  summarise(mean_price = weighted.mean(price, Prop)) %>%
  ungroup()

# write out
write.csv(price_dat_3, "data/price.csv", row.names = F)
  

