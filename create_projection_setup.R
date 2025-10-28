# Alberto Rovellini
# 10/23/2025
# This code builds the sets of input files needed to run the projection runs for the OY caps paper
# These runs have several dimensions:

# 4 caps (800, 600, 400, 200). The 800 run SHOULD result in the cap never getting applied, AKA the status quo
# 3 weight systems. If the 800K cap run is not constraining, we do not need to do it with 3 weight systems
# 4 climate scenarios: 3 SSPs + 1 with fixed, historical climate 
# total: (3*3*4)+(1*4) = 40 runs

# these dimensions will be defined by parameters in various files:

# Caps and weights are all done in the harvest.prm file
# change in f and start of HCR management at year 30 are all done in the harvest.prm files ( see older test runs to get this right )
# Climate scenarios will be determined by force.prm and will include salt, temp, hydro, and plankton scalar forcings
# runs will be 110 years: 30 years burn-in, 80 years projections
# The setup will take:
# 10 harvest files (3*3 + 1), but only if the 800K cap is unconstraining
# 4 force.prm files
# 40 sh files

library(tidyverse)

options(scipen = 999)

rm(list = ls())

grp <- read.csv("data/GOA_Groups.csv")
codes <- grp %>% pull(Code)
verts <- grp %>% filter(GroupType %in% c("MAMMAL","SHARK","BIRD","FISH")) %>% pull(Code)
hcr_grp <- c("POL", "COD", "POP", "SBF") # these are groups managed with HCR; they will need the FMSY proxy and have done the full mfc ramp
fmsy_grp <- c(hcr_grp, "HAL") # these are groups managed with HCR and halibut; they will need the FMSY proxy and have done the full mfc ramp
other_fmp_grp <- c("ATF", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "RFS", "RFP", "RFD", "THO", "DFS", "DFD", "SCU")
oy_grp <- c(hcr_grp, other_fmp_grp) # halibut is not a OY group

# begin from the SS model that is used to run the F ramps and determine the reference points
# that model had been built on run 2303, and it uses 1990s mFC values and 1990s climatologies
ss_model <- "AtlantisGOA_SS/"

# copy to a new model
# Copy a folder and all its contents
if(!dir.exists("AtlantisGOA_MS")){dir.create("AtlantisGOA_MS/")}
file.copy(from = list.files("AtlantisGOA_SS/", full.names = TRUE), 
          to = "AtlantisGOA_MS/", 
          recursive = TRUE)

# remove all harvest.prm and sh files
files_to_remove <- list.files("AtlantisGOA_MS/", 
                              pattern = "run_atlantis|GOA_harvest", 
                              full.names = TRUE, 
                              recursive = TRUE)
file.remove(files_to_remove)

# Harvest.prm -------------------------------------------------------------

# create template harvest.prm with correct parameters, mfc modifiers, etc.
# the start point is the same harvest.prm file used to create the set of mfc ramps used in the SS runs
file.copy("AtlantisGOAV0_02303/GOA_harvest_background.prm", "AtlantisGOA_MS/GOA_harvest_2303.prm")
harvest_base <- readLines("AtlantisGOA_MS/GOA_harvest_2303.prm")

# first, change do_sumB_HCR to 2
# all calibration runs have beeen using 0
# if the FlagSystCapSP vector is all set to 0, and if all groups have tierXXX 0, results should be identical (testing in progress)
do_sumB_HCR_line <- grep("do_sumB_HCR", harvest_base)
harvest_base[do_sumB_HCR_line] <- gsub("0", "2", harvest_base[do_sumB_HCR_line])

# turn on management. This will slow the model down even further, but it allows us to have the burn-in
# at year 30, the projections begin
# HCR species (POL, COD, POP, SBF) will use the HCR to determine F each year
# HAL will use FMSY as determined in the SS runs
# all other stocks will use fixed F, as the avg F for the 2015-2019 period
# deal with HCR stocks: activate the HCR at year 30
# This requires the management module to be running - need to set background_flagmanage 1
flagmanage_line <- grep("background_flagmanage", harvest_base)
harvest_base[flagmanage_line] <- gsub("\t0", "\t1", harvest_base[flagmanage_line])

# then set when you want the management for each fleet to kick in - we only have one fleet (background)
startmanage_line <- grep("background_start_manage", harvest_base)
harvest_base[startmanage_line] <- gsub("\t0", "\t10950", harvest_base[startmanage_line])

# need to setup the HCRs 
# this is going to use FMSY values from the SS runs for fref, and SSB0 for estBo
##############
ref_points <- read.csv("output/ref_points_from_SS_runs_oct2025.csv")

# clean up
ref_points <- ref_points %>%
  select(Code, b0, fref, fref_mfc_for_check) %>%
  drop_na() %>%
  rename(estBo = b0,
         Fref = fref)

# set up a data frame with: whichref, estBo, and Fref (if you used estCV and estBias they'd go here as well)
ref_points_hcr <- ref_points %>%
  filter(Code != "HAL") %>%
  mutate(whichref = 1,
         f_restart = 0.1) # this one ends up being Frescale, and as such it will be applied to the base mfc, which is for the 1990s. Nothing we can do about that.

df_hcr <- data.frame(Code = codes)
df_hcr <- df_hcr %>% 
  left_join(ref_points_hcr, by = "Code") 
df_hcr[is.na(df_hcr)] <- 0

##############
# set parameters
targ_refA_line <- grep("targ_refA", harvest_base)
harvest_base[targ_refA_line] <- gsub("\t0", "\t0.4", harvest_base[targ_refA_line])

targ_refB_line <- grep("targ_refB", harvest_base)
harvest_base[targ_refB_line] <- gsub("\t0", "\t0.4", harvest_base[targ_refB_line])

# targ_refE_line <- grep("targ_refE", harvest_base) # this one was already good
# harvest_base[targ_refE_line] <- gsub("\t0", "\t0.2", harvest_base[targ_refE_line])

lim_ref_line <- grep("lim_ref", harvest_base)
harvest_base[lim_ref_line] <- gsub("\t0", "\t0.02", harvest_base[lim_ref_line])

whichref_line <- grep("whichref	78", harvest_base)
harvest_base[whichref_line + 1] <- df_hcr %>% pull(whichref) %>% as.character() %>% paste(collapse = " ")

estBo_line <- grep("estBo	78", harvest_base)
harvest_base[estBo_line + 1] <- df_hcr %>% pull(estBo) %>% as.character() %>% paste(collapse = " ")

Fref_line <- grep("Fref	78", harvest_base) # reminder to self: this is called F but it will have to be an exploitation rate
harvest_base[Fref_line + 1] <- df_hcr %>% pull(Fref) %>% as.character() %>% paste(collapse = " ")

Frestart_scalar_line <- grep("Frestart_scalar	78", harvest_base)
harvest_base[Frestart_scalar_line + 1] <- df_hcr %>% pull(f_restart) %>% as.character() %>% paste(collapse = " ")

# now tiers
# pollock and code will use tier 14, which has the SSL 20% cutoff
# POP and SBF will use tier 1, which is the regular ramp 40/02
tierPOL_line <- grep("tierPOL", harvest_base)
harvest_base[tierPOL_line] <- gsub("\t0", "\t14", harvest_base[tierPOL_line])

tierCOD_line <- grep("tierCOD", harvest_base)
harvest_base[tierCOD_line] <- gsub("\t0", "\t14", harvest_base[tierCOD_line])

tierPOP_line <- grep("tierPOP", harvest_base)
harvest_base[tierPOP_line] <- gsub("\t0", "\t1", harvest_base[tierPOP_line])

tierSBF_line <- grep("tierSBF", harvest_base)
harvest_base[tierSBF_line] <- gsub("\t0", "\t1", harvest_base[tierSBF_line])

# need to do flagmfc change for everything that does not have the HCR
# the burn-in uses tuned mFC values hitting mean F from the 1990s
# the projection needs to use tuned mFC values hitting mean F from 2015-2019
# the run that has those properties is 2306
# the mFC in 2303 (and the PRM currently being modified) is the base value applied in the burn-in
# the mFC in 2306 is the target mFC to apply in proj
# set up a multiplier for the projection
# set all the necessary flags and switches for the mFC change

# TODO: all this will have to be tested with 1095 (3 years before it all kick in)
# I tested all this before and it should work but more testing will not hurt
file.copy("AtlantisGOAV0_02306/GOA_harvest_background.prm", "AtlantisGOA_MS/GOA_harvest_2306.prm")
harvest_target <- readLines("AtlantisGOA_MS/GOA_harvest_2306.prm")

# get mFC from current template PRM (the 1990s mFC to scale to proj values)
burnin_mfc <- list()
for(sp in codes) {
  
  mfc_line <- harvest_base[grep(paste0("mFC_", sp, " "), harvest_base) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  burnin_mfc[[sp]] <- data.frame("Code" = sp, "burnin_mfc"=  mfc)
}
burnin_mfc <- bind_rows(burnin_mfc)

# get mFC from target harvest (the 2015-2019 mFC)
# NB: Halibut is going to need a different approach. The target mFC will be the mFC associated with the FMSY run
# NB2: the HCR species will not be using these values in projection because the HCR handles F dynamically, but let's still set it up in case we need to do runs without the HCR
proj_mfc <- list()
for(sp in codes) {
  
  mfc_line <- harvest_target[grep(paste0("mFC_", sp, " "), harvest_target) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  proj_mfc[[sp]] <- data.frame("Code" = sp, "proj_mfc"=  mfc)
}
proj_mfc <- bind_rows(proj_mfc)

# plug in target value for HAL, which won't be using the HCR
proj_mfc[proj_mfc$Code == "HAL",]$proj_mfc <- ref_points %>% filter(Code == "HAL") %>% pull(fref_mfc_for_check)

# get scalars
scalars <- burnin_mfc %>% 
  left_join(proj_mfc, by = "Code") %>%
  mutate(scalar = proj_mfc / burnin_mfc) %>%
  mutate(scalar = ifelse(is.nan(scalar), 1, scalar))

# One unexpected effect is that even species that should have been using the same mfc in 1990 and 2015 have scalars != 1
# in some cases big ones too, like the usual culprit migrating groups
# this must have happened in the mfc tuning step, which must have gone in different directions in the two runs
# effectively, you retune the whole model so that f matches expectation 
# because salmon and hake are NOT relevant to this paper, I'd rather fix them to their 1990's F values than risk collapses that will then ripple through the model
# in general we'll have to look out for discontinuities after year 30
# many of these F's are small enough that doubling them should hardly have any effect
to_fix <- c("SCH", "SCO", "SCM", "SPI", "SSO", "HAK")
scalars <- scalars %>%
  mutate(scalar = ifelse(Code %in% to_fix, 1, scalar))

# however, I see that the step below only sets the multipliers on mfc for FMP species
# think on what makes the most sense here
# it probably hardly matters and it would be a safer approach

# write scalars to mFC change in the harvest.prm template and make all other relevant changes to activate mFC change
oy_grp_and_hal <- c(oy_grp, "HAL") # we need the mfc change setup for HAL as well

for (i in 1:length(oy_grp_and_hal)){
  
  sp <- oy_grp_and_hal[i]
  
  # what is the mfc_change factor to apply to reach the target mfc?
  mfc_change <- scalars %>% filter(Code == sp) %>% pull(scalar)
  
  # write all relevant lines in the PRM file
  # make changes to PRM lines
  # flagchangeF 1
  old_line <- harvest_base[grep("flagchangeF", harvest_base)]
  new_line <- gsub("\t0", "\t1", old_line)
  harvest_base[grep("flagchangeF", harvest_base)] <- new_line
  
  # flagFchange_XXX 1 for first element (which is BG and needs to go off)
  old_line <- harvest_base[grep(paste0("flagFchange_", sp), harvest_base) + 1]
  new_line <- sub("0", "1", old_line)
  harvest_base[grep(paste0("flagFchange_",sp), harvest_base) + 1] <- new_line
  
  # XXX_mFC_changes  - how many changes? Set to 1 for first element
  old_line <- harvest_base[grep(paste0(sp, "_mFC_changes"), harvest_base) + 1]
  new_line <- sub("0", "1", old_line)
  harvest_base[grep(paste0(sp, "_mFC_changes"), harvest_base) + 1] <- new_line
  
  # mFCchange_start_XXX set to 30*365
  old_line <- harvest_base[grep(paste0("mFCchange_start_", sp), harvest_base) + 1]
  change_start <- 365 * 30
  new_line <- sub("0", as.character(change_start), old_line)
  harvest_base[grep(paste0("mFCchange_start_", sp), harvest_base) + 1] <- new_line
  
  # mFCchange_period_XXX (1 day)
  old_line <- harvest_base[grep(paste0("mFCchange_period_", sp), harvest_base) + 1]
  new_line <- sub("0", "1", old_line)
  harvest_base[grep(paste0("mFCchange_period_", sp), harvest_base) + 1] <- new_line
  
  # mFCchange_mult_XXX set to mfc_change
  old_line <- harvest_base[grep(paste0("mFCchange_mult_", sp), harvest_base) + 1]
  new_line <- as.character(round(mfc_change, digits = 5))
  harvest_base[grep(paste0("mFCchange_mult_", sp), harvest_base) + 1] <- new_line
  
  # done
}

# fix the mEff lines
harvest_base <- gsub("0.25 0.25 0.25 0.25", "0 0 0 0", harvest_base, fixed = T)

# write out
writeLines(harvest_base, "AtlantisGOA_MS/GOA_harvest_template.prm")

# TODO: test template run without caps but with mfc change 

###############################
# Adding caps
# starting from the template, create files for the cap runs
# 3 caps: c(600K, 400K, 200K)
# 3 weigth schemes: Equal, binary, ramp
# for the weight schemes, see C:/Users/Alberto Rovellini/Documents/GOA/Paper3_projections/Atlantis_GOA_projections/oy_weigths_from_harvest_specs.R

# create 9 new PRM files with unique names
cap_df <- data.frame(cap = c(800, rep(c(600, 400, 200), each = 3)),
                     weights = c("equal", rep(c("equal","binary","ramp"), 3))) %>%
  mutate(filename = paste0("GOA_harvest_", cap, weights, ".prm"))

# weights
w_from_specs <- read.csv("data/w_from_harvest_specs.csv")
w_from_specs <- w_from_specs %>% select(Code, w)
# DFD is missing from these weights but part of the OY
# DFD includes grenadiers, but these do not show up in the harvest specs data, nor in the summarized catch accounting data
# according to the latest assessment, they have an estimated exploitation rate corresponding to F = 0.017
# this would insert them just above FHS in terms of exploitation, so in the low value stocks
w_from_specs <- rbind(w_from_specs, data.frame(Code = "DFD", w = 3)) %>%
  mutate(w = ifelse(Code %in% c("FFD","FHS","DFD"), w, w+1)) %>%
  arrange(-w)

# now construct the weight vectors for each scenario
w <- data.frame(Code = codes) %>%
  left_join(w_from_specs %>% select(Code, w), by = "Code") %>%
  mutate(w = replace_na(w, 0)) %>%
  mutate(equal = ifelse(w == 0, 0, 1),
         binary = case_when(
           w == 0 ~ 0,
           w <= 7 ~ 1,
           .default = 5
         )) %>% # 7 corresponds to ATF; the biggest gap in catch/TAC was between SKL and ATF
  rename(ramp = w)

# build PRM files
for(i in 1:nrow(cap_df)){
  
  this_harvest <- harvest_base
  
  filename <- cap_df$filename[i]
  oy_vec <- w %>% pull(equal) %>% as.character() %>% paste(collapse = " ")
  w_vec <- w %>% pull(cap_df$weights[i]) %>% as.character() %>% paste(collapse = " ")
  cap <- cap_df$cap[i] * 1000
  
  # modify lines
  FlagSystCapSP_line <- grep("FlagSystCapSP 78", this_harvest)
  this_harvest[FlagSystCapSP_line + 1] <- oy_vec
  
  SystCapSPpref_line <- grep("SystCapSPpref 78", this_harvest)
  this_harvest[SystCapSPpref_line + 1] <- w_vec
  
  Ecosystm_Cap_tonnes_line <- grep("Ecosystm_Cap_tonnes", this_harvest)
  this_harvest[Ecosystm_Cap_tonnes_line] <- paste("Ecosystm_Cap_tonnes", cap)
  
  flagSSBforHCR_line <- grep("flagSSBforHCR", this_harvest)
  this_harvest[flagSSBforHCR_line] <- gsub("0", "1", this_harvest[flagSSBforHCR_line])
  
  writeLines(this_harvest, paste0("AtlantisGOA_MS/", filename))
}

# TODO: do test that the mgmt and F change kicks in at the right time
# 2 runs: base / template file; and one of these (but move burnin! do like 3 years + 2)

# Force.prm ---------------------------------------------------------------
# 4 climate scenarios
# one with no climate at all and stable conditions - this one is the base force.prm
# use that as template to create the other ones
force_file <- "AtlantisGOA_MS/GOA_force.prm"
force_template <- readLines(force_file)

scenarios <- c(126, 245, 585)
burnin_t <- 30
proj_t <- 80
runtime <- burnin_t + proj_t

for(i in 1:length(scenarios)){
  
  clim <- scenarios[i]
  filename <- paste0("AtlantisGOA_MS/GOA_force_ssp", clim,".prm")
  this_force <- force_template
  
  # hydro
  start_h <- grep("nhdfiles", this_force) - 1
  end_h <- grep("hd0.name", this_force) + 1
  
  head_file <- this_force[1:start_h]
  tail_file <- this_force[end_h:length(this_force)]
  
  # make hydro files section
  h_section <- paste("nhdfiles", runtime)
  yr <- 2020
  for(j in 1:runtime){
    if(j <= burnin_t){
      string <- paste0("hd", j-1, ".name forcings/hydro/goa_hydro_1999.nc")
    } else {
      string <- paste0("hd", j-1, ".name forcings_proj/ssp", clim, "/hydro/goa_hydro_", yr, ".nc")
      yr <- yr + 1
    }
    h_section <- c(h_section, string)
  }
  
  # bind after each step
  this_force <- c(head_file, h_section, tail_file)
  
  # now temperature
  start_t <- grep("ntempfiles", this_force) - 1
  end_t <- grep("Temperature0.name", this_force) + 1
  
  head_file <- this_force[1:start_t]
  tail_file <- this_force[end_t:length(this_force)]
  
  # make hydro files section
  t_section <- paste("ntempfiles", runtime)
  yr <- 2020
  for(j in 1:runtime){
    if(j <= burnin_t){
      string <- paste0("Temperature", j-1, ".name forcings/temp/mean_1990s_temperature.nc")
    } else {
      string <- paste0("Temperature", j-1, ".name forcings_proj/ssp", clim, "/temp/goa_roms_temp_", yr, ".nc")
      yr <- yr + 1
    }
    t_section <- c(t_section, string)
  }
  
  # bind after each step
  this_force <- c(head_file, t_section, tail_file)
  
  # salinity
  start_s <- grep("nsaltfiles", this_force) - 1
  end_s <- grep("Salinity0.name", this_force) + 1
  
  head_file <- this_force[1:start_s]
  tail_file <- this_force[end_s:length(this_force)]
  
  # make hydro files section
  s_section <- paste("nsaltfiles", runtime)
  yr <- 2020
  for(j in 1:runtime){
    if(j <= burnin_t){
      string <- paste0("Salinity", j-1, ".name forcings/salt/mean_1990s_salinity.nc")
    } else {
      string <- paste0("Salinity", j-1, ".name forcings_proj/ssp", clim, "/salt/goa_roms_salt_", yr, ".nc")
      yr <- yr + 1
    }
    s_section <- c(s_section, string)
  }
  
  # bind after each step
  this_force <- c(head_file, s_section, tail_file)
  
  # plankton
  start_p <- grep("use_external_scaling", this_force) - 1
  
  head_file <- this_force[1:start_p]
  tail_file <- c(
    "use_external_scaling 1",
    "",
    "scale_all_mortality 0",
    "mortality_addition 0",
    "",
    paste0("externalBiologyForcingFile scalar_proj_ROMS_ssp", clim, "_110yr.nc"),
    "",
    "externalBiologyForcingFile_rewind 0"
  )
  
  this_force <- c(head_file, tail_file)

  # write force file
  writeLines(this_force, filename)
}

# sh files ----------------------------------------------------------------

# now build 40 sh files that will be used to launch these runs
# Create shell scripts for all combinations of harvest and force files

# Define force files
force_files <- c("GOA_force.prm", "GOA_force_ssp126.prm", "GOA_force_ssp245.prm", 
                 "GOA_force_ssp585.prm")

# Create all combinations
combinations <- expand.grid(
  harvest_file = cap_df$filename,
  force_file = force_files,
  stringsAsFactors = FALSE
)

# Add run numbers (000-039)
combinations$run_id <- sprintf("%03d", 0:(nrow(combinations) - 1))

# Add cap and weights info by merging back with cap_df
combinations <- combinations %>% left_join(cap_df, by = c("harvest_file" = "filename"))

# Create shell scripts for each combination
for (i in 1:nrow(combinations)) {
  run_id <- combinations$run_id[i]
  harvest <- combinations$harvest_file[i]
  force <- combinations$force_file[i]
  
  # Create shell script content matching the example format
  sh_content <- paste0(
    "atlantisMerged -i GOA_cb_summer.nc  0 -o outputGOA_", run_id, ".nc ",
    "-r GOA_run.prm -f ", force, " -p GOA_physics.prm -b GOAbioparam.prm ",
    "-h ", harvest, " -m GOAMigrations.csv -s GOA_Groups.csv ",
    "-q GOA_fisheries.csv -d outputFolder", run_id, "\n"
  )
  
  # Write shell script file
  sh_filename <- paste0("AtlantisGOA_MS/run_atlantis_", run_id, ".sh")
  writeLines(sh_content, sh_filename)
  
  # Make it executable (on Unix-like systems)
  Sys.chmod(sh_filename, mode = "0755")
}

# Write tracking CSV
write.csv(combinations, "AtlantisGOA_MS/atlantis_run_combinations.csv", row.names = FALSE)

# Print summary
cat("Created", nrow(combinations), "shell scripts (GOA_run_000.sh to GOA_run_039.sh)\n")
cat("Tracking file saved as: atlantis_run_combinations.csv\n")
cat("\nFirst few combinations:\n")
print(head(combinations, 10))
