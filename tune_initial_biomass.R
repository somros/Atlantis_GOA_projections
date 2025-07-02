# Alberto Rovellini
# 07/01/2025
# GOAL: take groundfish biomass closer to 2020 values than 1990 (calibration) values
# WHY: because if we are running projections starting in 2020 people will look at starting biomass before we begin the projections
# KEY ISSUE: we are not doing this with any other group, so the model would not be representative of 2020 conditions but some sort of hybrid

# Workflow:
# read in stock biomass from assessments
# compute mean 2015-2019 biomass
# compute correction factors for initial biomass
# next step will be bhalpha calibration

# do everything on the alaska only biomass
# this way it does not matter whether it refers to females only, how many age classes, etc.
# it is only a relative scaling
# Assumptions: constant sex ratios, constant age class composition, constant proportions between AK and BC
# These are all false assumptions, but it helps keep things tractable and we are dealing with averages

# for now using the biomass from authors we got way back when, but by now those are old values

library(tidyverse)
library(readxl)

# read groups and get codes for ordering
grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% select(Code)

biom <- read_xlsx("data/biomass_from_stock_assessment.xlsx", range = "A1:AM51")

biom_long <- biom %>%
  pivot_longer(-Year, names_to = "Stock", values_to = "biomass")

# map to Atlantis groups
key <- data.frame("Stock" = unique(biom_long$Stock),
                  "Code" = c("POP","FHS","POL","RFP","RFS","RFS","SBF","FFS","FFS","FFS","FFS","FFS","FFS","FFS",
                             "COD","DFD","ATF","HAL","REX","FFD","THO","DOG","SKB","SKL","SKO","RFS","SCU",
                             "SCU","SCU","SCU","RFD","RFD","RFD","RFD","RFD","RFD","RFD","RFD"))

biom_long <- biom_long %>%
  left_join(key, by = "Stock")

# add up across atlantis stocks
biom_long <- biom_long %>%
  group_by(Year, Code) %>%
  summarise(biomass = sum(biomass, na.rm = T)) %>%
  ungroup()

# 1990 values
biom_1990 <- biom_long %>% filter(Year == 1990) %>% select(-Year)
# 2015-2019 means
biom_end <- biom_long %>% filter(Year >= 2015, Year <= 2019) %>%
  group_by(Code) %>%
  summarise(biomass_end = mean(biomass, na.rm = T))

corr <- biom_1990 %>%
  left_join(biom_end, by = "Code") %>%
  mutate(corr = biomass_end / biomass)

# write vector
corr_vec <- codes %>%
  left_join(corr) %>%
  pull(corr) %>%
  signif(3)

corr_vec[is.na(corr_vec)] <- 1
