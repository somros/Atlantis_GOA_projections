# debug COD productivity
library(tidyverse)

# focus on run 2254:
# highest F, lowest steepnedd (4x BHbeta)
# see if we can get a different value of f

# then repeat for 2252 (original beta)

yr_end <- 30
avg_period <- 5

simdir <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_2252"

biom <- read.delim(here::here(simdir, "outputGOA02252_testAgeBiomIndx.txt"), sep = " ")
catch <- read.delim(here::here(simdir, "outputGOA02252_testCatch.txt"), sep = " ")

# subset time early on if necessary
biom <- biom %>% filter(Time/365 <= yr_end)
catch <- catch %>% filter(Time/365 <= yr_end)

# Process biomass data once, outside the loop
# Convert to long format once and use faster string splitting
biom_long <- biom %>%
  pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")

code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
biom_long$Code <- sapply(code_age_split, `[`, 1)
biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))

# Pre-filter to only relevant species and times
biom_filtered <- biom_long %>%
  filter(Code == "COD") %>%#,
  #Time > 0, 
  #Time < max(Time)) %>%
  select(-Code.Age)  # Remove original column

# Process all species at once using vectorized operations
# Create biomass summaries for all species
biom_selex_all <- biom_filtered %>%
  left_join(selex_df, by = "Code") %>%
  filter(Age >= startage) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_selex = sum(mt), .groups = 'drop')

# inspect the biomass output
# biom_selex_all is what is used to compute f
biom_selex_all %>%
  mutate(Year = Time / 365) %>%
  filter(Year > 25) %>% # using the last 5 years for the f calcs
  ggplot(aes(x = Year, y = biom_mt_selex))+
  geom_line()+
  theme_bw()

# there are very substantial swings in biomass throughout the year!!
# this does highlight the problem of excessive stock productivity, but it does mean that computed F may be very different

# let's use the maximum biomass record for the year instead of the mid-winter data point
# you need to sum over age classes first, or else you may double count biomass over the year 
biom_selex_old <- biom_filtered %>%
  left_join(selex_df, by = "Code") %>%
  filter(Age >= startage) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_selex = sum(mt), .groups = 'drop') %>%
  mutate(Year = Time / 365) %>%
  filter(Year %% 1 == 0) %>%
  rename(winter_mt = biom_mt_selex) %>%
  select(-Time)

biom_selex_max <- biom_filtered %>%
  left_join(selex_df, by = "Code") %>%
  filter(Age >= startage) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_selex = sum(mt), .groups = 'drop') %>%
  mutate(Year = Time / 365) %>%
  mutate(Year = floor(Year)) %>%
  group_by(Year, Code) %>%
  slice_max(biom_mt_selex) %>%
  ungroup() %>%
  rename(max_mt = biom_mt_selex) %>%
  select(-Time)

# join these two
biom_selex_all <- biom_selex_old %>%
  left_join(biom_selex_max, by = c("Year","Code")) %>%
  select(Year, Code, winter_mt, max_mt)

# jopin catch
# Process catch data for all species
catch_long <- catch %>%
  select(Time, all_of("COD")) %>%
  pivot_longer(-Time, names_to = "Code", values_to = "catch_mt") %>%
  mutate(Year = Time / 365) %>%
  select(-Time)


# join operation
# to compute F from exploitation rates we need to stagger catch and biomass
# catch at t356 / biomass at t0
# this is not so important past the first year
sp_data <- biom_selex_all %>%
  # filter(Code == sp) %>%
  left_join(catch_long, by = c("Year", "Code")) %>%
  # left_join(filter(biom_tot_all, Code == sp), by = c("Time", "Code")) %>%
  mutate(
    mu_old = catch_mt / lag(winter_mt),
    f_old = -log(1 - mu_old),
    mu_max = catch_mt / lag(max_mt),
    f_max = -log(1 - mu_max)
  ) %>%
  filter(Year > (max(Year) - avg_period)) %>%
  group_by(Code) %>%
  summarize(across(winter_mt:f_max, ~mean(.x, na.rm = T))) %>%
  ungroup()


# this is SO MUCH CLOSER to where it should be
# data points MATTER
# there is enough selected biomass to support the catch
# however, the problem of excessive steepness stands, I believe
# bring in the key - this is run 27 for COD
sp_data <- sp_data %>% left_join(key %>% filter(idx == 27), by = c("Code"))

# get input target F (i.e. the F that the input mFC corresponds to)
sp_data <- sp_data %>%
  mutate(targ_f = -365 * log(1 - targ_mfc))

# perceived F is still higher than input f, but soooo much closer to it than before
# still higher because mFC acts every day, and eats into biomass continuously
# so even picking the highest value (e.g. after a year class aged into the fishery), you still lose some
# that's fine

# This looks very promising and could mean we solved at least one problem
# the issue remains that the stock is too productive, which is why the work on beta is important
# let's figure out the fraction of unfished at the highest F as well
biom_spawn_old <- biom_filtered %>%
  left_join(fspb_df, by = c("Code", "Age"="age")) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_spawn = sum(mt * fspb), .groups = 'drop') %>%
  mutate(Year = Time / 365) %>%
  filter(Year %% 1 == 0) %>%
  rename(winter_mt = biom_mt_spawn) %>%
  select(-Time)

biom_spawn_max <- biom_filtered %>%
  left_join(fspb_df, by = c("Code", "Age"="age")) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_spawn = sum(mt * fspb), .groups = 'drop') %>%
  mutate(Year = Time / 365) %>%
  mutate(Year = floor(Year)) %>%
  group_by(Year, Code) %>%
  slice_max(biom_mt_spawn) %>%
  ungroup() %>%
  rename(max_mt = biom_mt_spawn) %>%
  select(-Time)

# join these two
biom_spawn_all <- biom_spawn_old %>%
  left_join(biom_spawn_max, by = c("Year","Code")) %>%
  select(Year, Code, winter_mt, max_mt)

# get SSB0, just for the max
simdir_b0 <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_2246"
biom_b0 <- read.delim(here::here(simdir_b0, "outputGOA02246_testAgeBiomIndx.txt"), sep = " ")
# subset time early on if necessary
biom_b0 <- biom_b0 %>% filter(Time/365 <= yr_end)

# Process biomass data once, outside the loop
# Convert to long format once and use faster string splitting
biom_long_b0 <- biom_b0 %>%
  pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")

biom_long_b0$Code <- sapply(code_age_split, `[`, 1)
biom_long_b0$Age <- as.numeric(sapply(code_age_split, `[`, 2))

# Pre-filter to only relevant species and times
biom_filtered_b0 <- biom_long_b0 %>%
  filter(Code == "COD") %>%#,
  select(-Code.Age)  # Remove original column

ssb0 <- biom_filtered_b0 %>%
  left_join(fspb_df, by = c("Code", "Age"="age")) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_spawn = sum(mt * fspb), .groups = 'drop') %>%
  mutate(Year = Time / 365) %>%
  mutate(Year = floor(Year)) %>%
  group_by(Year, Code) %>%
  slice_max(biom_mt_spawn) %>%
  ungroup() %>%
  rename(max_mt = biom_mt_spawn) %>%
  select(-Time) %>%
  filter(Year > (max(Year) - avg_period)) %>%
  group_by(Code) %>%
  summarize(b0 = mean(max_mt)) %>%
  ungroup() # 248608. mt

biom_spawn_max %>%
  filter(Year > (max(Year) - avg_period)) %>%
  group_by(Code) %>%
  summarize(b0 = mean(max_mt)) %>%
  ungroup() # 61170.


# 68087 / 252097

# using the highest biomass logic, depletion at 4x "FMSY" is 0.246
# This sounds high, HOWEVER, the value of FMSY = 0.253 used here for COD was from the OY paper, which was a different model
# we tuned COD a lot to resemble recent assessments. compared to assessment F, perceived F = 1.25 is more like ~2.5x FOFL
# considering that FMSY intended as proxy of FOFL should occur at B40%... we are off, but not as far off as I originally thought!
# in summary, the stock is too productive but tuning beta may be setting us on the right track

# in the run with original beta: depletion = 0.27
# run with 4x original beta: depletion = 0.246
# so moving beta moves things the right direction I think
