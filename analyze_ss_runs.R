# ALberto Rovellini
# 08/14/2025
# Purpose:
# Read biomass.txt and catch.txt output files from the 85 single-species simulations
# Extract terminal biomass and catch
# Plot yield curves and depletion curves; overlay Rovellini et al. (2025) values, and stock assessment FOFL values 
# write out table with B0 (all species) and FMSY (HCR species)
library(tidyverse)
library(patchwork)
library(here)
library(viridis)

rm(list = ls())

# Data needed -------------------------------------------------------------

# functional groups
grp <- read.csv("data/GOA_Groups.csv")
verts <- grp %>% filter(GroupType %in% c("MAMMAL","SHARK","BIRD","FISH")) %>% pull(Code)
fmsy_grp <- c("POL", "COD", "POP", "SBF", "HAL") # these are groups managed with HCR and halibut; they will need the FMSY proxy and have done the full mfc ramp
other_fmp_grp <- c("ATF", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "RFS", "RFP", "RFD", "THO", "DFS", "DFD", "SCU")
oy_grp <- c(fmsy_grp, other_fmp_grp)
# for debugging only delete when done
# oy_grp <- c("POL", "COD", "POP", "HAL")

# read maturity ogives as presented in the biol.prm file
bio_prm <- "data/GOAbioparam_SS.prm"
bio <- readLines(bio_prm)

fspb_df <- data.frame()
for(i in 1:length(oy_grp)){
  
  sp <- oy_grp[i]
  fspb_line <- bio[grep(paste0("FSPB_", sp), bio) + 2]
  fspb <- as.numeric(strsplit(fspb_line, split = " ")[[1]])
  
  fspb_sp <- data.frame("Code" = rep(sp, length(fspb)), 
                        "age" = 0:(length(fspb)-1), 
                        "fspb" = fspb)
  
  fspb_df <- rbind(fspb_df, fspb_sp)
  
}

# read selectivity start ages from one of the 85 prm files - they are consistent across files
harvest_file <- "output/single_species_runs/mfc_ramps_NEW/GOA_harvest_000.prm"
harvest <- readLines(harvest_file)

selex_df <- data.frame()
for(i in 1:length(oy_grp)) {
  
  sp <- oy_grp[i]
  startage_line <- harvest[grep(paste0(sp, "_mFC_startage"), harvest) + 1]
  startage <- as.numeric(strsplit(startage_line, split = " ")[[1]])[1]
  
  selex_df <- rbind(selex_df, data.frame(
    "Code" = sp,
    "startage" = startage
  ))
  
}

# key of mfc applied to each run
key <- read.csv("output/single_species_runs/mfc_ramps_NEW/reference_table.csv") # fix rthe space
# debugging
#key <- key %>% filter(idx<47)

# fmsy values from Rovellini et al. (2025)
fmsy_oy <- read.csv("data/fmsy_Rovellini_2025.csv")
fmsy_oy <- fmsy_oy %>% left_join(grp %>% select(Code, LongName))

# output file lists
biom_files <- list.files("output/single_species_runs/ss-80yr-october/", pattern = "BiomIndx", full.names = T, recursive = T)
catch_files <- list.files("output/single_species_runs/ss-80yr-october/", pattern = "Catch", full.names = T, recursive = T)

# run information (change for real runs)
yr_end <- 80
burnin <- 30
avg_period <- 10

# Biomass -----------------------------------------------------------------

make_ss_df <- function(this_run, ssb = T){
  
  #this_run <- 3
  
  print(this_run)
  
  this_run_char <- sprintf("%03d", this_run)
  
  # File paths
  biom_file <- biom_files[grep(this_run_char, biom_files)]
  catch_file <- catch_files[grep(this_run_char, catch_files)]
  
  # Read files once
  biom <- read.csv(biom_file, sep = " ", header = T)
  catch <- read.csv(catch_file, sep = " ", header = T)
  
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
    filter(Code %in% oy_grp) %>%#,
    #Time > 0, 
    #Time < max(Time)) %>%
    select(-Code.Age)  # Remove original column
  
  # Process all species at once using vectorized operations
  # Create biomass summaries for all species
  # Albi: changed code so that it uses max biomass the year prior to account for recruitment / aging dynamics
  # NB: right now we have one time step per year so this should be identical to how it was before
  biom_selex_max <- biom_filtered %>%
    left_join(selex_df, by = "Code") %>%
    filter(Age >= startage) %>%
    group_by(Time, Code) %>%
    summarise(biom_mt_selex = sum(mt), .groups = 'drop') %>% # sum across cohort for total biomass for a group
    mutate(Year = Time / 365) %>%
    mutate(Year = floor(Year)) %>%
    group_by(Year, Code) %>%
    slice_max(biom_mt_selex) %>%
    ungroup() %>%
    mutate(Time = Year * 365) %>% # restore the Time column
    select(-Year)
  
  if(ssb){
    
    biom_tot_all <- biom_filtered %>%
      left_join(fspb_df, by = c("Code", "Age"="age")) %>%
      group_by(Time, Code) %>%
      summarise(biom_mt_tot = sum(mt * fspb), .groups = 'drop') %>%
      mutate(Year = Time / 365) %>%
      mutate(Year = floor(Year)) %>%
      group_by(Year, Code) %>%
      slice_max(biom_mt_tot) %>% # get the highest biomass record for the year
      ungroup() %>%
      mutate(Time = Year * 365) %>% # restore the Time column
      select(-Year)
    
  } else {
    
    biom_tot_all <- biom_filtered %>%
      group_by(Time, Code) %>%
      summarise(biom_mt_tot = sum(mt), .groups = 'drop')
    
  }
  
  # Process catch data for all species
  catch_long <- catch %>%
    select(Time, all_of(oy_grp)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt")
  
  # join operation
  # to compute F from exploitation rates we need to stagger catch and biomass
  # catch at t356 / biomass at t0
  # this is not so important past the first year
  result_list <- list()
  for(sp in oy_grp) { # to make sure we don't mess up the staggering, do one sp at a time
    
    # print(sp)
    
    sp_data <- biom_selex_max %>%
      filter(Code == sp) %>%
      left_join(filter(catch_long, Code == sp), by = c("Time", "Code")) %>%
      left_join(filter(biom_tot_all, Code == sp), by = c("Time", "Code")) %>%
      mutate(
        mu = catch_mt / lag(biom_mt_selex),
        f = -log(1 - mu)
      ) %>%
      mutate(year = Time / 365) %>%
      filter(year > (max(year) - avg_period)) %>%
      group_by(Code) %>%
      summarize(across(biom_mt_selex:f, ~mean(.x, na.rm = T))) %>%
      ungroup()
    
    result_list[[sp]] <- sp_data
  }
  
  res_df <- do.call(rbind, result_list)
  
  res_df$run <- this_run
  
  return(res_df)
}

ss_df_tmp <- bind_rows(lapply(key$idx, make_ss_df))

# lots of noise here, need to bring in the run indices
ss_df <- ss_df_tmp %>% left_join(key, by = c("run" = "idx"))

# subset to the species of interest based on Code.y
ss_df <- ss_df %>%
  mutate(tmp = ifelse(Code.x == Code.y, 1, 0)) %>%
  filter(tmp > 0) %>%
  select(-Code.y, -tmp) %>%
  rename(Code = Code.x)

# get input target F (i.e. the F that the input mFC corresponds to)
# there are very large differences between input F and realized F
ss_df <- ss_df %>%
  mutate(targ_f = -365 * log(1 - targ_mfc))

# get b0, b40, bmsy frames
b0_df <- ss_df %>% 
  filter(f == 0) %>% 
  select(Code, biom_mt_tot) %>%
  rename(b0 = biom_mt_tot)

b40_df <- ss_df %>%
  filter(Code %in% fmsy_grp) %>%
  left_join(b0_df, by = "Code") %>%
  mutate(prop_biom = biom_mt_tot / b0) %>%
  group_by(Code) %>%
  slice_min(abs(prop_biom - 0.4), n=  1, with_ties = F) %>%
  select(Code, biom_mt_tot) %>%
  rename(b40 = biom_mt_tot)

fmsy_ss <- ss_df %>% 
  filter(Code %in% fmsy_grp) %>%
  group_by(Code) %>%
  slice_max(catch_mt) %>%
  select(Code, f, targ_f, run) %>%
  rename(fmsy_realized = f,
         fmsy_input = targ_f)

# add OFL estimates from most recent full assessment instead of old Atlantis estimates
ss_ofl_2025 <- data.frame(
  Code = c("POL","COD","SBF","POP","HAL"),
  fofl = c(0.307, 0.57, 0.094, 0.1, NA)
)

# Curves ------------------------------------------------------------------
# this is only for species for which we need FMSY

p_biom <- ss_df %>%
  filter(Code %in% fmsy_grp) %>%
  left_join(b40_df, by = "Code") %>%
  left_join(fmsy_ss, by = "Code") %>%
  left_join(fmsy_oy %>% select(Code, atlantis_fmsy), by = "Code") %>%
  left_join(ss_ofl_2025, by = "Code") %>%
  filter(!is.na(f)) %>%  # for plotting purposes - ideally at some point you should no longer have NaNs
  filter(Code != "HAL") %>% # HAL does not need the HCR so exclude it from this plot
  ggplot(aes(x = f, y = biom_mt_tot / 1000, color = targ_f))+
  geom_point()+
  scale_color_viridis()+
  geom_hline(aes(yintercept = b40 / 1000), linetype = "dashed", color = "darkgrey")+
  geom_vline(aes(xintercept = fmsy_realized), linetype = "dashed", color=  "blue")+
  geom_vline(aes(xintercept = atlantis_fmsy), linetype = "dashed", color=  "orange")+
  geom_vline(aes(xintercept = fofl), linetype = "dashed", color = "red")+
  theme_bw()+
  guides(color = "none")+  # Remove the color legend for points
  scale_y_continuous(limits = c(0,NA))+
  labs(x = "Fishing mortality", y = "Spawning biomass (1000 mt)", color = "Input F")+
  facet_wrap(~Code, scales = "free", ncol = 1)
p_biom

# ggsave

p_catch <- ss_df %>%
  filter(Code %in% fmsy_grp) %>%
  left_join(fmsy_ss, by = "Code") %>%
  left_join(fmsy_oy %>% select(Code, atlantis_fmsy), by = "Code") %>%
  left_join(ss_ofl_2025, by = "Code") %>%
  filter(!is.na(f)) %>%
  filter(Code != "HAL") %>% # HAL does not need the HCR so exclude it from this plot
  ggplot(aes(x = f, y = catch_mt / 1000))+
  geom_point(aes(color = targ_f))+
  scale_color_viridis(name = "Input F")+
  geom_vline(aes(xintercept = fmsy_realized, linetype = "This Atlantis FMSY"), color = "blue")+
  geom_vline(aes(xintercept = atlantis_fmsy, linetype = "Rovellini et al. 2025 FMSY"), color = "orange")+
  geom_vline(aes(xintercept = fofl, linetype = "FOFL latest full assessment"), color = "red")+
  scale_linetype_manual(name = "Reference lines",
                        values = c("This Atlantis FMSY" = "dashed", 
                                   "Rovellini et al. 2025 FMSY" = "dashed", 
                                   "FOFL latest full assessment" = "dashed"))+
  theme_bw()+
  scale_y_continuous(limits = c(0,NA))+
  labs(x = "Fishing mortality", y = "Catch (1000 mt)")+
  facet_wrap(~Code, scales = "free", ncol = 1)
p_catch

# Combine plots side by side
combined_plot <- p_biom | p_catch
combined_plot

ggsave(here("plots", "combined_plots_oct2025.png"), combined_plot, width = 7, height = 5, dpi = 600)

# This looks a lot better than the earlier plots, meaning that tuning steepness was useful
# some comparisons of reference points with the real world

# Write output ------------------------------------------------------------

# produce an output table for b0 and fmsy, which will be used to build the HCR parameters for the ms projections
# NOTE: fref as used in Atlantis is an exploitation rate
# It will have to be the input exploitation rate corresponding to FMSY
# It should be equivalent to mfc*365 for that run
# b0 is in SSB so should be able to use it as-is, if we use the SSB switch for the HCR

out_tab <- b0_df %>%
  left_join(fmsy_ss, by = "Code") %>%
  mutate(fref = 1 - exp(-fmsy_input),
         fref_mfc_for_check = 1 - exp(-fmsy_input/365)) # these track

write.csv(out_tab, "output/ref_points_from_SS_runs_oct2025.csv")

# Diagnostics -------------------------------------------------------------

# check the relationship between input F and realized F
ss_df %>%
  filter(Code %in% fmsy_grp) %>%
  ggplot(aes(x = targ_f, y = f))+
  geom_point()+
  theme_bw()+
  geom_abline(color = "red")+
  labs(x = "Input F", y = "Realized F")+
  facet_grid(~Code)


# Biomass by age class ----------------------------------------------------

get_biomass_by_age <- function(this_run){
  
  #this_run <- 3
  
  print(this_run)
  
  this_run_char <- sprintf("%03d", this_run)
  
  # File paths
  biom_file <- biom_files[grep(this_run_char, biom_files)]
  
  # Read files once
  biom <- read.csv(biom_file, sep = " ", header = T)
  
  # subset time early on if necessary
  biom <- biom %>% filter(Time/365 <= yr_end)
  
  # Process biomass data once, outside the loop
  # Convert to long format once and use faster string splitting
  biom_long <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")
  
  code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
  biom_long$Code <- sapply(code_age_split, `[`, 1)
  biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
  
  # Pre-filter to only relevant species and times
  biom_filtered <- biom_long %>%
    filter(Code %in% oy_grp) %>%#,
    select(-Code.Age)  # Remove original column
  
  # Here you cannot do the slice_max trick because you will double-count biomass as you're keeping the age classes
  # Instead, get the average biomass of an age class during a certain year
  biom_tot_all <- biom_filtered %>%
    mutate(Year = Time / 365) %>%
    mutate(Year = floor(Year)) %>%
    group_by(Year, Code, Age) %>%
    summarize(mean_mt = mean(mt, na.rm = T)) %>% # get the highest biomass record for the year
    ungroup() %>%
    filter(Year > max(Year) - avg_period) %>%
    group_by(Code, Age) %>%
    summarise(mean_mt_terminal = mean(mean_mt, na.rm = T)) %>%
    ungroup()
  
  # turn to proportions
  biom_prop <- biom_tot_all %>%
    group_by(Code) %>%
    mutate(totbiom = sum(mean_mt_terminal)) %>%
    ungroup() %>%
    mutate(prop = mean_mt_terminal / totbiom) %>%
    select(Code, Age, prop) %>%
    mutate(run = this_run)
  
  return(biom_prop)
}

age_df_tmp <- bind_rows(lapply(key$idx, get_biomass_by_age))

# lots of noise here, need to bring in the run indices
age_df <- age_df_tmp %>% left_join(key, by = c("run" = "idx"))

# subset to the species of interest based on Code.y
age_df <- age_df %>%
  mutate(tmp = ifelse(Code.x == Code.y, 1, 0)) %>%
  filter(tmp > 0) %>%
  select(-Code.y, -tmp) %>%
  rename(Code = Code.x)

# get input target F (i.e. the F that the input mFC corresponds to)
# there are very large differences between input F and realized F
age_df <- age_df %>%
  mutate(targ_f = -365 * log(1 - targ_mfc))

p_agebiom <- age_df %>%
  filter(Code %in% fmsy_grp) %>%
  ggplot(aes(x = mult, y = prop, fill = Age))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis()+
  theme_bw()+
  labs(x = "Fishing mortality multiplier", y = "Proportion of total biomass", fill = "Age class")+
  facet_wrap(~Code, scales = "free", ncol = 1)
p_agebiom
