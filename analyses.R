# Alberto Rovellini
# 11/03/2025
# Plots and analyses of runs using the ecosystem cap / OY algorithm

library(tidyverse)
library(ggh4x)
library(patchwork)
library(PNWColors)
library(tidync)
library(ncdf4)
library(data.table)

rm(list = ls())

# General info ------------------------------------------------------------

# set path for reference run
ref_run <- "000" # check that this makes sense
oy_dir <- "AtlantisGOA_MS/"
output_dir <- "../v2"

# source functions script
source("functions_YEAR.R")

# handle time
burnin <- 30
biom_file <- paste0(output_dir, "/outputGOA_", ref_run, "AgeBiomIndx.txt")
biom <- read.csv(biom_file, sep = " ", header = T)
yr_end <- ceiling(max(unique(biom$Time)))/365 # 110

# identify boundary boxes - will need this for NAA extraction
# TODO: this can't be run on the server now, problem with rbgm
# fl <- 'data/GOA_WGS84_V4_final.bgm'
# bgm <- rbgm::read_bgm(fl)
# goa_sf <- rbgm::box_sf(bgm)
# boundary_boxes <- goa_sf %>% sf::st_set_geometry(NULL) %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes
boundary_boxes <- c(0, 2, 8, 11, 30, 53, 58, 62, 69, 81, 89, 95, 97, 108)

# read in species info
grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)

##############################################################
# species and fleets
# pull these from a reference OY run - comparing them makes sense only if all species are consistent
harvest <- readLines(paste(oy_dir, "GOA_harvest_800equal.prm", sep = "/"))
cap_vec <- as.numeric(unlist(strsplit(harvest[grep("FlagSystCapSP", harvest)+1], split = " ")))

# get species that are managed under the OY
oy_species <- codes[which(cap_vec>0)]
oy_names <- grps %>% filter(Code %in% oy_species) %>% pull(Name) # need this for the NC files

oy_fleets <- "background" # manually set this, it will just be bg for the foreseeable future

# get species that are managed with the HCRs
hcr_spp <- c()
for(sp in codes){
  # Look for "tierXX\t" pattern to match the tab-delimited format
  pattern <- paste0("tier", sp, "\t")
  matches <- harvest[grep(pattern, harvest, fixed = TRUE)]
  if(length(matches) > 0) {
    tier <- as.numeric(gsub("[^0-9]", "", matches[1]))
    if(tier > 0){hcr_spp <- c(hcr_spp, sp)}
  }
}

# get "preferred" species
# ALBI: check that we still need this
# TODO: base this on harvest specification data and keep it constant across scenarios
pref <- data.frame("Code" = oy_species) %>%
  rowwise()%>%
  mutate(pref = ifelse(Code %in% c("POL","COD","SBF","POP"),1,0))

##############################################################
# fishery information
# TODO: move estbo info into fishery PRM file
# estbo_files <- list.files("data/estBo/", full.names = T)
# estbo_list <- list()
# for(i in 1:length(estbo_files)){
#   estbo_list[[i]] <- read.csv(estbo_files[i])
# }
# estbo_key <- bind_rows(estbo_list) %>% select(Code, mean_biom) %>% rename(estbo = mean_biom)

ref_points_ss <- read.csv("output/ref_points_from_SS_runs_nov2025.csv")
estbo_key <- ref_points_ss %>% select(Code, b0) %>% rename(estbo = b0)
# estbo_key <- estbo_key %>% mutate(b40 = 0.4*estbo,
#                                   b25 = 0.25*estbo)

# need to get fref
# the value in the prm is an input that may or may not be close
# for HCR species, get from reference points. THOSE ARE MU, so turn to F
fref_frame <- ref_points_ss %>%
  select(Code, fref) %>%
  mutate(fref = -log(1 - fref)) %>%
  filter(Code != "HAL") %>%
  drop_na()

# NB: 11/04/2025
# An error in the parametrization of the first MS batch caused Fref to be off for the HCR species.
# basically I left the mfcchange_mult on for the HCR species, thinking that the HCR would overwrite that
# it did not. mfc, after the hcr, would get multiplied by that factor
# this caused the HCR to be off. way off for COD
# for the purpose of plotting relative to Fref, multiply Fref by the mult factors
# TODO: drop this when we fix that in future runs

# mult_df <- data.frame("Code" = c("POL", "COD", "POP", "SBF"),
#                       mult = c(1.06868, 2.67137, 1.22662, 0.84828))
# fref_frame <- fref_frame %>%
#   left_join(mult_df, by = "Code") %>%
#   mutate(fref = fref * mult) %>%
#   select(-mult)
# ALBI 3/6 we no longer need the bit above, leaving for reference but delete when cleaning up

# NB: this is a hack purely for plotting purpose to see if the HCR works - catch is all off

# for others, need to get it from external analysis of recent F
fref_other <- read.csv("data/fref_from_assessments_proj.csv") %>%
  filter(!Code %in% fref_frame$Code) %>%
  select(Code, `F`) %>%
  rename(fref = `F`)

fref_frame <- rbind(fref_frame, 
                    fref_other,
                    data.frame("Code" = "DFS", # handle DFS but we'll need something else
                               fref = NA))

##############################################################
# biology information
# maturity at age
bio_prm <- paste0("AtlantisGOA_MS/GOAbioparam.prm")
bio <- readLines(bio_prm)

fspb_df <- data.frame()
for(i in 1:length(oy_species)){
  
  sp <- oy_species[i]
  fspb_line <- bio[grep(paste0("FSPB_", sp), bio) + 2]
  fspb <- as.numeric(strsplit(fspb_line, split = " ")[[1]])
  
  fspb_sp <- data.frame("Code" = rep(sp, length(fspb)), 
                        "age" = 0:(length(fspb)-1), 
                        "fspb" = fspb)
  
  fspb_df <- rbind(fspb_df, fspb_sp)
  
}

# diets of pollock's predators
dietfile <-  paste0(output_dir, "/outputGOA_", ref_run, "DietCheck.txt")
diet <- fread(dietfile, 
              select = c("Time", "Predator", "Cohort", "POL"))

# identify POL's predators: say all groups who have >5% POL in their diets at the end of the burn-in period
# these are the species that have at least one cohort eating > 5% POL by the end of the burn-in
POL_predators <- diet %>%
  mutate(Time = ceiling(Time/365)) %>%
  filter(Time == burnin) %>%
  group_by(Predator, Cohort) %>%
  summarise(POL = mean(POL)) %>%
  filter(POL>0.05) %>%
  dplyr::select(Predator, Cohort) %>%
  distinct() %>%
  mutate(isPred = 1)

##############################################################
# make a directory to store plots
# make a plot directory
plotdir <- paste0("plots/oy/",  format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
if(!dir.exists(plotdir)){
  dir.create(plotdir, recursive = T)
} else {
  print("This directory exists")
}

# Run properties ----------------------------------------------------------
# get this from the MS run key
run_combs <- read.csv("AtlantisGOA_MS/atlantis_run_combinations.csv")
run <- run_combs %>% pull(run_id) %>% sprintf("%03d", .)

# rename "ramp" to "attainment-based" which is more informative
run_combs$weights <- gsub("ramp", "attainment-based", run_combs$weights)

# format this consistently with the key used before to limit changes to the code below
key_config <- run_combs %>% 
  mutate(env = gsub (".prm", "", gsub("GOA_force_", "", force_file))) %>% # climate scenario
  mutate(cap = cap * 1000, # transform to tons
         run_id = sprintf("%03d", run_id), # run as character
         env = ifelse(env == "GOA_force", "NoClimate", env)) %>%
  select(run_id, cap, weights, env) %>%
  rename(run = run_id, # naming consistency
         wgts = weights)

# Apply function to process run output ---------------------------------------------------------

# TODO: fix so that we compute F same as in the SS runs
# Probably won't make a difference but we should be consistent

catch_df_tmp <- bind_rows(lapply(run, pull_fishery_info, ssb = T))

# bind to key, and to group names
catch_df <- catch_df_tmp %>%
  left_join(key_config, by = "run") %>%
  left_join(grps %>% select(Code, Name), by = "Code")

# turn caps into factors for better plotting
catch_df$cap <- as.character(catch_df$cap)
catch_df$cap[is.na(catch_df$cap)] <- "No cap"
catch_df$cap <- factor(catch_df$cap, levels = c("8e+05","6e+05","4e+05","2e+05"))

# order weigths
catch_df$wgts <- factor(catch_df$wgts, levels = c("equal","binary","attainment-based"))

# ALBI: delete the below
# get preferred species given a certain weighting scheme
# using > mean(w) here in an attempt to automate the choice
# pref <- catch_df %>%
#   select(Code, w, wgts) %>%
#   distinct() %>%
#   group_by(wgts) %>%
#   mutate(mean_w = mean(w)) %>%
#   rowwise() %>%
#   mutate(pref = ifelse(w>=mean_w,1,0))%>%
#   ungroup()%>%
#   filter(pref>0)%>%
#   select(Code,wgts,pref)

# this does not really work for the equal weight scenario though. The preferred species would need to be the same across scenarios for meaningful comparisons
# as a second option, set the preferred species to the usual suspects, which we did above

# Apply function to produce fishery plots -------------------------------------------------------------------

plot_fishery(catch_df %>% filter(wgts != "binary"))

# Plot OY rescaling only --------------------------------------------------

resc_plot_sp <- c("POL", "COD","ATF")

cap_col <- pnw_palette(name="Sunset2",n=length(unique(catch_df$cap)),type="discrete")

resc_plot <- catch_df %>%
  filter(wgts != "binary") %>%
  filter(env == "NoClimate") %>%
  filter(Time > 0) %>%
  filter(Time >= burnin*365) %>%
  filter(Code %in% resc_plot_sp) %>%
  filter(!is.na(catch_mt)) %>%
  ggplot(aes(x = (Time/365)+1990, y = oy_rescale, color = factor(cap), linetype = factor(wgts))) +
  # annotate("rect", xmin = 0+1990, xmax = burnin+1990, ymin = -Inf, ymax = Inf, 
  #          fill = "grey", alpha = 0.3) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(values = cap_col)+
  scale_linetype_manual(values = c("solid","dashed","dotted"))+
  scale_x_continuous(breaks = c(seq(1990,2100,10))) +
  scale_y_continuous(limits = c(0,NA)) +
  facet_grid(~ Name) +
  labs(x = "", 
       y = "",
       color = "Cap (mt)",
       linetype = "Weight scheme") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(plotdir, "/", "resc_plot.png"), resc_plot, width = 10, height = 3.5, dpi = 300, units = "in")

# Shannon Index -----------------------------------------------------------
h_frame <- data.frame() # initialize df
H <- bind_rows(lapply(run, get_H, do_mature=F)) # apply function to the nc files and get time series for all runs

# join with catch df from above
h_plot <- catch_df %>%
  filter(!is.na(f)) %>%
  mutate(Time = Time/365) %>%
  left_join(H, by = c("Time"="year","Name","run")) 

# make plots
for(i in 1:length(oy_species)){
  
  current_code <- oy_species[i]
  current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
  
  p6 <- h_plot %>%
    filter(Time >= burnin + 5) %>%
    filter(Name == current_name) %>% 
    filter(wgts != "binary") %>%
    ggplot(aes(x = f, y = H, color = Time, shape = wgts))+
    geom_point(aes(shape = factor(wgts)), size = 1)+
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    scale_color_viridis_c(option = "viridis")+
    #scale_y_continuous(limits = c(0,NA))+
    theme_bw()+
    labs(x = "F", 
         y = "H", 
         color = "Year",
         shape = "Weight scheme",
         title = paste("Shannon index for", current_name))+
    facet_grid(factor(env)~cap)
  
  ggsave(paste0(plotdir, "/h/", current_code, "_hcr.png"), p6, 
         width = 10, height = 4.5, 
         units = "in", dpi = 300)
  
}

# Diets -------------------------------------------------------------------
# Proportion of pollock in the diet of its key predators
# do for all runs
diets_pol_pred <- bind_rows(lapply(run, get_polprop))

# join with run info - also need POL f info
diet_key <- catch_df %>%
  filter(!is.na(f), Code == "POL") %>%
  mutate(Time = Time/365) %>%
  select(Time,f,run,cap,wgts,env)

diet_plot <- diets_pol_pred %>%
  left_join(diet_key, by = c("Time","run")) %>%
  filter(!is.na(cap)) # this happens because the last record of f, which you filter for, is for year 99

for(i in 1:length(unique(POL_predators$Predator))){
  
  current_code <- unique(POL_predators$Predator)[i]
  current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
  
  p7 <- diet_plot %>%
    #filter(Time >= burnin) %>%
    filter(Time > 5) %>%
    filter(Time %% 2 == 0) %>% # filter to even years only to make the plot less confusing - for some reason points fluctuate up and down every 2 years
    filter(Predator == current_code) %>%
    filter(wgts != "binary") %>%
    ggplot(aes(x = Time, y = POL, color = f, shape = wgts))+
    annotate("rect", xmin = 0, xmax = burnin, ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.3) +
    geom_point(aes(shape = factor(wgts)), size = 1)+
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    scale_color_viridis_c(option = "cividis")+
    #geom_vline(xintercept = burnin, linetype = "dashed", color = "red")+
    #scale_y_continuous(limits = c(0,NA))+
    theme_bw()+
    labs(x = "Year", 
         y = "Proportion of pollock in diet", 
         color = "F(pol)",
         shape = "Weight scheme",
         title = current_name)+
    facet_grid(factor(env)~cap)
  
  ggsave(paste0(plotdir, "/diet/", current_code, "_diet.png"), p7, 
         width = 10, height = 4.5, 
         units = "in", dpi = 300)
  
}

# This is tied to the climate scenario to an extent, but it's also just POL declining over time in the base model
# Differences overall are very small

# Ecosystem indicators ----------------------------------------------------

# Calculate indicators for all runs (with progress)
indicators_all <- map_df(run, ~{
  cat("Processing run", .x, "\n")
  calc_ecosystem_indicators(.x)
})

# Create plots
plot_ecosystem_indicators(indicators_all)

# biomass of predatory and piscivorous fish
# first of all we need to identify which species are predatory / piscivorous
# piscivorous means it eats >x% fish in the realized diets
# predatory means it eats >x% non-plankton 
# we then use these lists to compute biomass of piscivores and predators 
fish_grps <- grps %>% filter(GroupType %in% c("FISH","SHARK")) %>% pull(Code)
non_plkt_grps <- grps %>% filter(!GroupType %in% 
                                   c("LG_INF",
                                     "SM_INF",
                                     "PHYTOBEN",
                                     "MED_ZOO",
                                     "SM_ZOO",
                                     "LG_PHY",
                                     "SM_PHY",
                                     "SED_BACT",
                                     "PL_BACT",
                                     "CARRION",
                                     "LAB_DET",
                                     "REF_DET")) %>% pull(Code)
  
get_diet_base <- function(this_run){
  
  print(this_run)
  
  # File paths
  #wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  wd <- "../v2"
  dietfile <-  paste0(wd, "/outputGOA_", this_run, "DietCheck.txt")
  diet <- fread((dietfile)) %>% select(-Stock, -Updated) %>% filter(Predator %in% fish_grps) 
  
  # filter to year 30, end of the burn-in
  diet_base <- diet %>%
    pivot_longer(-c(Time:Cohort), names_to = "prey", values_to = "prop") %>%
    mutate(year = ceiling(Time/365)) %>%
    filter(year == burnin) %>%
    group_by(Predator, Cohort, prey) %>%
    summarise(mean_prop = mean(prop))
  
  return(diet_base)
} 

diet_30 <- get_diet_base("000") # avg diet comps in year 30 (end of burn-in)

piscivores <- diet_30 %>%
  group_by(Predator, Cohort) %>%
  # keep only the oldest (max) cohort per predator, with the idea that their diets reflect the adult diets
  group_by(Predator) %>%
  filter(Cohort == max(Cohort)) %>%
  # sum proportions for fish_grps prey items
  summarise(
    fish_prop = sum(mean_prop[prey %in% fish_grps]),
    .groups = "drop"
  ) %>%
  mutate(is_piscivore = ifelse(fish_prop > 0.5, 1, 0)) %>% # arbitrary cutoff of 50% as majority of diet
  filter(is_piscivore > 0) %>%
  pull(Predator) #  [1] "ATF" "BDF" "BSF" "DOG" "DOL" "KWR" "PIN" "SHD" "SHP" "SSL" "WHT"

predators <- diet_30 %>%
  group_by(Predator, Cohort) %>%
  # keep only the oldest (max) cohort per predator
  group_by(Predator) %>%
  filter(Cohort == max(Cohort)) %>%
  # sum proportions for non_plkt_grps prey items
  summarise(
    pred_prop = sum(mean_prop[prey %in% non_plkt_grps]),
    .groups = "drop"
  ) %>%
  mutate(is_predator = ifelse(pred_prop > 0.5, 1, 0))  %>%
  filter(is_predator > 0) %>%
  pull(Predator)

# [1] "ATF" "BDF" "BDI" "BSF" "BSI" "COD" "CRO" "DFD" "DFS" "DOG" "DOL" "EUL" "FFS" "FHS" "FOS" "HAL" "HER" "KIN" "KWR" "KWT" "OCT" "PIN" "POL" "POP" "RFD" "RFP"
# [27] "RFS" "SBF" "SCO" "SCU" "SHD" "SHP" "SKB" "SKL" "SKO" "SSL" "THO" "WHB" "WHH" "WHT"

# the definition of piscivore is simpler than "predator", in some way everything is a predator except for true plant grazers / detritivores
# the proportion is also arbitrary
# check with Isaac / Kerim


# Revenue -----------------------------------------------------------------

# Load price data
price_dat <- read.csv("data/price_safe.csv")

# For all runs
all_runs <- key_config %>% pull(run)

revenue_all <- map_df(all_runs, ~calc_revenue(.x, price_dat))

# turn caps into factors for better plotting
revenue_all$cap <- as.character(revenue_all$cap)
revenue_all$cap[is.na(revenue_all$cap)] <- "No cap"
revenue_all$cap <- factor(revenue_all$cap, levels = c("8e+05","6e+05","4e+05","2e+05"))

# order weigths
revenue_all$wgts <- factor(revenue_all$wgts, levels = c("equal","binary","attainment-based"))

# Create plot
plot_revenue(revenue_all %>% filter(wgts != "binary"))

# Tradeoff plots ----------------------------------------------------------------

plot_delta(catch_df, revenue_df = revenue_all)

plot_ecosystem_delta(catch_df)

# LTL ---------------------------------------------------------------------

ltl_biomass_all <- map_df(run, calc_ltl_biomass, groups = c("EUP", "ZL", "ZM", "PL", "PAN", "PWN"))

# turn caps into factors for better plotting
ltl_biomass_all$cap <- as.character(ltl_biomass_all$cap)
ltl_biomass_all$cap[is.na(ltl_biomass_all$cap)] <- "No cap"
ltl_biomass_all$cap <- factor(ltl_biomass_all$cap, levels = c("8e+05","6e+05","4e+05","2e+05"))

# order weigths
ltl_biomass_all$wgts <- factor(ltl_biomass_all$wgts, levels = c("equal","binary","attainment-based"))

# Create plot
plot_ltl_biomass(ltl_biomass_all %>% 
                   filter(wgts != "binary") %>%
                   filter(Time %% 365 == 0))

# WAA / NAA ---------------------------------------------------------------------

# this takes time

# For all runs with specific species
these_names <- c("Pollock", "Cod", "Arrowtooth_flounder", "Sablefish", "Pacific_ocean_perch" )

waa_core <- map_df(run, ~calc_weight_at_age(.x, sp_names = these_names, boundary_boxes = boundary_boxes))
naa_core <- map_df(run, ~calc_numbers_at_age(.x, sp_names = these_names, boundary_boxes = boundary_boxes))

# drop year 110, because there is only one time step there and it messes up computations
waa_core <- waa_core %>% filter(year < max(year))
naa_core <- naa_core %>% filter(year < max(year))

# plot
# Comparisons plotted against these defaults:
# Env: NoClimate
# Cap: 800K, unless you are comparing wgts scheme
# Wgts: equal and 400K cap (no wgts at 800K) 

plot_age_heatmap(waa_core, naa_core, plot_type = "both", by_env = T, tag = "core", which_decade = 10)
plot_age_heatmap(waa_core, naa_core, plot_type = "both", by_cap = T, tag = "core", which_decade = 10)
plot_age_heatmap(waa_core, naa_core, plot_type = "both", by_wgts = T, tag = "core", which_decade = 10)

# plot_age_heatmap(waa_core, plot_type = "WAA", by_env = T, tag = "core")
# plot_age_heatmap(waa_core, plot_type = "WAA", by_cap = T, tag = "core")
# plot_age_heatmap(waa_core, plot_type = "WAA", by_wgts = T, tag = "core")
# 
# plot_age_heatmap(naa_df = naa_core, plot_type = "NAA", by_env = T, tag = "core")
# plot_age_heatmap(naa_df = naa_core, plot_type = "NAA", by_cap = T, tag = "core")
# plot_age_heatmap(naa_df = naa_core, plot_type = "NAA", by_wgts = T, tag = "core")

# all grps for appendix
waa_all <- map_df(run, ~calc_weight_at_age(.x, sp_names = oy_names, boundary_boxes = boundary_boxes))
naa_all <- map_df(run, ~calc_numbers_at_age(.x, sp_names = oy_names, boundary_boxes = boundary_boxes))

plot_age_heatmap(waa_all, naa_all, plot_type = "both", by_env = T, tag = "all", which_decade = 10)
plot_age_heatmap(waa_all, naa_all, plot_type = "both", by_cap = T, tag = "all", which_decade = 5)
plot_age_heatmap(waa_all, naa_all, plot_type = "both", by_wgts = T, tag = "all", which_decade = 10)
plot_age_heatmap(waa_all, naa_all %>% filter(Name != "Cod"), plot_type = "both", by_env = T, tag = "all_nocod", which_decade = 10)

# do it for top predators
# For all runs with specific species
preds <- c("Steller_sea_lion", "Pinnipeds", "Seabird_dive_fish", "Seabird_surface_fish", "Dolphins")
waa_pred <- map_df(run, ~calc_weight_at_age(.x, sp_names = preds, boundary_boxes = boundary_boxes))

# plot
plot_age_heatmap(waa_pred, plot_type = "WAA", by_env = T, tag = "pred", which_decade = 10)
pw <- plot_age_heatmap(waa_pred, plot_type = "WAA", by_cap = T, tag = "pred", which_decade = 5) # for figure
plot_age_heatmap(waa_pred, plot_type = "WAA", by_wgts = T, tag = "pred", which_decade = 5)

# Top pred diets -----------------------------------------------------------------

pred_diets <- bind_rows(lapply(run, get_dietcomp_preds, which_decade = 5))

# All predators by cap
pd <- plot_diet_heatmap(pred_diets, by_cap = TRUE, tag = "test")

# stack for figure
p_pred <- pw / pd +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

ggsave(paste0(plotdir, "/heatmaps/stacked_preds.png"), p_pred,
       width = 10, height = 5.5, units = "in", dpi = 300)


# Cap computations --------------------------------------------------------
# look at aggregate catch at the end of the run (year 109)
catch_df %>%
  filter(Time == 109*365, run == "000") %>%
  pull(catch_mt) %>%
  sum()

# 363563 mt. This is in the neighborhood of Mueter and Megrey
# This is saying here is the catch the system equilibrates to if HCR stocks get harvested at Ftarget, and other stock keep being exploited at current levels
# save this

# Tcorr tracker -----------------------------------------------------------
# To corroborate the effects of Tcorr on WAA for groundfish, make plots that show where along the curve we are
# ingredients: 
# 1. mean output temperature from the NC files over the time series. maybe summer temp only
# 2. Unimodal response parameters for all fish groups
tcorr_dat <- read.csv("data/tcorr_pars.csv", header = TRUE)

# subset to groundfish as that's what we are showing
tcorr_dat <- tcorr_dat %>%
  left_join(grps %>% select(Code, LongName), by = c("Species"="LongName")) %>%
  filter(Code %in% oy_species) %>%
  select(-Code)

# Build curves across a Tamb range for every species
Tamb_seq <- seq(5, 12, by = 0.1)

curves_df <- tcorr_dat %>%
  expand_grid(Tamb = Tamb_seq) %>%                        # expand here
  mutate(                                                  # fully vectorised, no rowwise
    Y     = log(q10) * (Tmax - Topt + 2),
    Z     = log(q10) * (Tmax - Topt),
    X     = (Z^2 * (1 + (1 + 40/Y)^0.5)^2) / 400,
    V     = (Tmax - Tamb) / (Tmax - Topt),
    Tcorr = V^X * exp(X * (1 - V))
  ) %>%
  select(Species, Tamb, Tcorr)

# extract temperature from the output
# Apply across all runs
# READ box_areas.csv. We need box areas to weigh temperatures by box areas in getting model-wide averages
# as we have trouble installing rbgm on this server, I am importing this as CSV that I produced locally. areas are derived from the BGM file
# the data file will be in the repo so it can be reran
box_areas <- read.csv("data/box_areas.csv")
temp_all <- map_df(run, ~calc_mean_temperature(.x, boundary_boxes = boundary_boxes, layer = 7)) # layer 7 is the sediment, which has the same temp as the botttom water layer

temp_df <- temp_all %>% 
  filter(cap == 800000, wgts == "equal", env != "NoClimate") %>%  # subset to the status quo cap
  mutate(year = floor(t) + 1990) %>%
  select(-run, -cap, -wgts) %>%
  filter(year < 2100, year > 2019) %>% # drop last time step as that's winter of the last year 
  group_by(year, env) %>%
  summarise(mean_temp = mean(mean_temp)) # get mean temperature for each year, understanding that there is seasonal variation

# temp_df %>%
#   ggplot(aes(x = year, y = mean_temp, color = env))+
#   geom_line()

# Compute Tcorr at actual model temperatures for every species × env × year
model_tcorr <- temp_df %>%
  cross_join(tcorr_dat) %>%
  mutate(
    Y     = log(q10) * (Tmax - Topt + 2),
    Z     = log(q10) * (Tmax - Topt),
    X     = (Z^2 * (1 + (1 + 40/Y)^0.5)^2) / 400,
    V     = (Tmax - mean_temp) / (Tmax - Topt),
    Tcorr = V^X * exp(X * (1 - V))
  )

# Plot
p_tcorr <- ggplot() +
  geom_line(
    data = curves_df,
    aes(x = Tamb, y = Tcorr, group = Species),
    color = "grey75", linewidth = 1, alpha = 0.8
  ) +
  geom_point(
    data = model_tcorr %>% filter(env == "ssp585"),
    aes(x = mean_temp, y = Tcorr, color = year),
    size = 1.7
  ) +
  scale_color_viridis_c(name = "Year") +
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  facet_wrap(~ Species, ncol = 4, labeller = labeller(Species = ~gsub(" - ", "\n", .x)))+
  labs(x = "Temperature (°C)", y = expression(T[corr])) +
  theme_bw()
p_tcorr

ggsave(paste0(plotdir, "/tcorr.png"), p_tcorr,
       width = 8.5, height = 8, units = "in", dpi = 600)

