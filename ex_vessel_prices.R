# ALberto Rovellini
# 07/07/2025
# a look at the ex-vessel prices
# These are prices published by NOAA in the Federal Register annually to determine observer fees
# They come from fish ticket data and market data bases
# Downloaded from https://www.fisheries.noaa.gov/action/north-pacific-observer-program-standard-ex-vessel-prices-groundfish-and-halibut-federal
# I picked 2020, which is based on volumes and values 2017-2019
# I digitized Table 1 and kept only GOA-wide values
# this table does not have sablefish and halibut IFQ and CDQ prices

library(tidyverse)
dat <- read.delim("data/ex_vessel_2020.txt", sep = ",")

# note that these figures reflect volume and value, hence the difference in price between gears
# for the purposes of this, let's average across gears
dat1 <- dat %>%
  select(-Port.area) %>%
  pivot_longer(-Species..species.code., names_to = "gear", values_to = "price") %>%
  group_by(Species..species.code.) %>%
  summarise(mean_price = mean(price, na.rm = T)) %>%
  arrange(-mean_price)

dat1

# this paints a different picture compared to the attainment lens
# ex-vessel prices in 2017-2019 were highest for sablefish, then lots of rockfish follows
# much of this rockfish was in RFD, which scores 10th in attainment
# thornyhead is 3rd highest and scored 11th in attainment
# where it really breaks down is where you have cod and POP lower than rex sole and skates
# pollock is the third lowest price, just one notch up from ATF
#

# All in all I am finding actual ex-vessel prices useful to understand something about the attainament but not too helpful to determine weights
# pollock is a very large volume and it is processed industrially, so the price is low 
# the atf that is attained is still sold it seems, also for industrial processing and possible overseas
# market saturation is possibly important too
# flatfish is not targeted, thus has low attainment, but what is caught can command relatively high prices
# if more was caught however, its price may go down
# so using ex-vessel prices without accounting for volume and attainment seems misleading

# Ono et al. (2018)
# In this paper they use perceived and realize net prices for pollock, cod, yellowfin sole, halibut
# but it is not from data, instead they estimate it so that quotas and catches predicted by their model are closest to 2010-2014 observations
# in other words they think of price as what motivates catch allocations and then catch itself

