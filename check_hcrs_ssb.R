# Alberto Rovellini
# 98/04/2025
# Analyses of flagSSBforHCR 1
# Done without cap and for 3 climate pathways

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
ref_run <- 2100 # for now this does not really make sense - this run has OY and HCR management starting on day 1
oy_dir <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", ref_run)

# source functions script
source("functions.R")

# handle time
burnin <- 15
biom_file <- paste0("outputGOA0", ref_run, "_testAgeBiomIndx.txt")
biom <- read.csv(paste(oy_dir, biom_file, sep = "/"), sep = " ", header = T)
#yr_end <- ceiling(max(unique(biom$Time)))/365
yr_end <- 90 # if the end of the run should be shorter, for example for runs that had been cut short, specify it here


# identify boundary boxes - will need this for NAA extraction
fl <- 'data/GOA_WGS84_V4_final.bgm'
bgm <- rbgm::read_bgm(fl)
goa_sf <- rbgm::box_sf(bgm)
boundary_boxes <- goa_sf %>% sf::st_set_geometry(NULL) %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes

# read in species info
grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)

##############################################################
# species and fleets
# pull these from a reference OY run - comparing them makes sense only if all species are consistent
harvest_prm <- list.files(oy_dir)[grep("GOA_harvest_.*.prm", list.files(oy_dir))]
harvest <- readLines(paste(oy_dir, harvest_prm, sep = "/"))
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
# TODO: base this on harvest specification data and keep it constant across scenarios
pref <- data.frame("Code" = oy_species) %>%
  rowwise()%>%
  mutate(pref = ifelse(Code %in% c("POL","COD","SBF","POP"),1,0))

##############################################################
# fishery information
# TODO: move estbo info into fishery PRM file
estbo_files <- list.files("data/estBo/", full.names = T)
estbo_list <- list()
for(i in 1:length(estbo_files)){
  estbo_list[[i]] <- read.csv(estbo_files[i])
}
estbo_key <- bind_rows(estbo_list) %>% select(Code, mean_biom) %>% rename(estbo = mean_biom)

##############################################################
# biology information
# maturity at age
bio_prm <- list.files(oy_dir)[grep("GOAbioparam_.*.prm", list.files(oy_dir))]
bio <- readLines(paste(oy_dir, bio_prm, sep = "/"))

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
dietfile <-  paste0("outputGOA0", ref_run, "_testDietCheck.txt")
diet <- fread(paste(oy_dir, dietfile, sep = "/"), 
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
plotdir <- paste0("plots/oy/", ref_run, "_", Sys.Date())
if(!dir.exists(plotdir)){
  dir.create(plotdir, recursive = T)
} else {
  print("This directory exists")
}

# Run properties ----------------------------------------------------------
# model runs
#run <- c(2060, 2061, 2064, 2065, 2066, 2067, 2072, 2073, 2074, 2075)
#run <- c(2072, 2073, 2074, 2075, 2080, 2081, 2082, 2083)
#run <- c(2097:2112, 2124:2143)
run <- c(2100,2104,2127,2222:2224)

# caps
#cap <- c(400, 200, 600, 400, 200, 600, NA, 400, 200, 600) * 1000
#cap <- rep(c(NA, 400, 200, 600), 2) * 1000
#cap <- rep(c(200,400,600,NA), 9) * 1000
cap <- c(rep("Tot",3), rep("SSB",3)) # use this for the HCRs

# weight scheme
#wgts <- c(rep("equal", 3), rep("binary", 4), rep("ramp", 3))
#wgts <- c(rep("ramp", 8))
#wgts <- c(rep("equal", 8), rep("binary", 8), rep("equal", 4), rep("binary", 4), rep("ramp", 12))
wgts <- rep("equal",6)

# climate forcings
# env <- c(rep(1999, 4), rep("2014", 4))
# env <- c(rep(c(rep("ssp126", 4), rep("ssp585", 4)), 2),
#          rep("ssp245", 8),
#          c(rep("ssp126",4), rep("ssp585",4), rep("ssp245", 4)))
env <- c(rep(c("ssp126","ssp585","ssp245"),2))

# combine all run factors in a label key
set_key <- function(run, cap = NULL, wgts = NULL, env = NULL, other = NULL) { # gamma = NULL, 
  # Get the length of run vector
  n <- length(run)
  
  # Set defaults for all arguments except run
  if (is.null(cap)) cap <- rep(NA, n)
  if (is.null(wgts)) wgts <- rep(NA, n)
  #if (is.null(gamma)) gamma <- rep(NA, n) # this is only useful for HCR shape 5
  if (is.null(env)) env <- rep(NA, n)
  if (is.null(other)) other <- rep(NA, n) # this could be rec devs or or off, etc
  
  key_config <- data.frame(run, cap, wgts, env, other) # gamma, 
  return(key_config)
  
}

key_config <- set_key(run, cap, wgts, env)

# Apply function to process run output ---------------------------------------------------------

catch_df <- bind_rows(lapply(run, pull_fishery_info))

# bind to key, and to group names
catch_df <- catch_df %>%
  left_join(key_config, by = "run") %>%
  left_join(grps %>% select(Code, Name), by = "Code")

# turn caps into factors for better plotting
catch_df$cap <- as.character(catch_df$cap)
catch_df$cap[is.na(catch_df$cap)] <- "No cap"

# order weigths
catch_df$wgts <- factor(catch_df$wgts, levels = c("equal","binary","ramp"))

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

plot_fishery(catch_df)

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
    filter(Time >= burnin) %>%
    filter(Name == current_name) %>%
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
# do for all runs
diets_pol_pred <- bind_rows(lapply(run, get_polprop))

# join with run info - also need POL f info
diet_key <- catch_df %>%
  filter(!is.na(f), Code == "POL") %>%
  mutate(Time = Time/365) %>%
  select(Time,f,run,cap,wgts,env,other)

diet_plot <- diets_pol_pred %>%
  left_join(diet_key, by = c("Time","run")) %>%
  filter(!is.na(cap)) # this happens because the last record of f, which you filter for, is for year 99

for(i in 1:length(unique(POL_predators$Predator))){
  
  current_code <- unique(POL_predators$Predator)[i]
  current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
  
  p7 <- diet_plot %>%
    filter(Time >= burnin) %>%
    filter(Predator == current_code) %>%
    ggplot(aes(x = Time, y = POL, color = f, shape = wgts))+
    geom_point(aes(shape = factor(wgts)), size = 1)+
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    scale_color_viridis_c(option = "cividis")+
    #scale_y_continuous(limits = c(0,NA))+
    theme_bw()+
    labs(x = "Year", 
         y = "Proportion of pollock in diet", 
         color = "F(pol)",
         shape = "Weight scheme",
         title = current_name)+
    facet_grid(factor(env)~cap)
  
  ggsave(paste0(plotdir, "/diet/", current_code, "_hcr.png"), p7, 
         width = 10, height = 4.5, 
         units = "in", dpi = 300)
  
}

# NB: at the moment there is a pervasive time series effect that is masking most of these plots
# This is tied to the climate scenario to an extent, but it's also just POL declining over time in the base model
