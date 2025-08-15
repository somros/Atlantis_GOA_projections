# Alberto Rovellini
# 08/15/2025
# build sets of PRM files for the MS runs
# write some intro and rationale

options(scipen = 999)

library(tidyverse)
grp <- read.csv("data/GOA_Groups.csv")
codes <- grp %>% pull(Code)
verts <- grp %>% filter(GroupType %in% c("MAMMAL","SHARK","BIRD","FISH")) %>% pull(Code)
hcr_grp <- c("POL", "COD", "POP", "SBF") # these are groups managed with HCR; they will need the FMSY proxy and have done the full mfc ramp
fmsy_grp <- c(hcr_grp, "HAL") # these are groups managed with HCR and halibut; they will need the FMSY proxy and have done the full mfc ramp
other_fmp_grp <- c("ATF", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "RFS", "RFP", "RFD", "THO", "DFS", "DFD", "SCU")
oy_grp <- c(hcr_grp, other_fmp_grp) # halibut is not a OY group

# 3 caps
# 3 weight systems
# 4 climate scenarios 
# no cap, 4 climate scenarios
# total: (3*3*4)+(1*4) = 40 runs

# caps and weight systems are all done in the PRM files
# climate scenarios will be determined by force.prm and will include salt, temp, hydro, and plankton forcings
# runs will be 115 years: 85 years projections, 30 years burn-in
# change in f and start of HCR management at year 30 ( see older test runs to get this right )
# this will be:
# 10 harvest files
# 4 force.prm files
# 40 sh files
# 1 biol and everything else

# begin from the SS model that is used to run the F ramps and determine the reference points
# that model had been built on run 2238, and it uses 1990s mFC values and 1990s climatologies
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
file.copy("AtlantisGOAV0_02238/GOA_harvest_background.prm", "AtlantisGOA_MS/GOA_harvest_2238.prm")
harvest_base <- readLines("AtlantisGOA_MS/GOA_harvest_2238.prm")

# first, change do_sumB_HCR to 2
# all calibration runs have beeen using 0
# if the FlagSystCapSP vector is all set to 0, and if all groups have tierXXX 0, results should be identical (testing in progress)
do_sumB_HCR_line <- grep("do_sumB_HCR", harvest_base)
harvest_base[do_sumB_HCR_line] <- gsub("0", "2", harvest_base[do_sumB_HCR_line])

# turn on management. This will slow the model down even further, but it allows us to have the burn-in
# at year 30, the projections begin
# HCR species (POL, COD, POP, SBF) will use the HCR to determine F each year
# HAL will use FMSY as determined in the SS runs
# all other stocks will use fixed F
# deal with HCR stocks: activate the HCR at year 30
# This requires the management module to be running - need to set background_flagmanage 1
flagmanage_line <- grep("background_flagmanage", harvest_base)
harvest_base[flagmanage_line] <- gsub("\t0", "\t1", harvest_base[flagmanage_line])

# then set when you want the management for each lfeet to kick in - we only have one fleet (background)
startmanage_line <- grep("background_start_manage", harvest_base)
harvest_base[startmanage_line] <- gsub("\t0", "\t10950", harvest_base[startmanage_line])

# need to setup the HCRs 
# this will need to be done using FMSY as Fref and B0 as estbo when the SS runs are done
# for now do a placeholder setup
##############
# set up a data frame with: whichref, estBo, and Fref (if you used estCV and estBias they'd go here as well)
ref_points <- data.frame(Code = c("POL","COD","POP","SBF"),
                         whichref = c(1, 1, 1, 1),
                         estBo = c(2e+06, 2e+05, 1e+05, 1e+05), # placeholders!
                         Fref = c(0.3, 0.25, 0.1, 0.1), # placeholders! remember this has to be a mu
                         f_restart = c(0.1, 0.1, 0.1, 0.1)) 

df_hcr <- data.frame(Code = codes)
df_hcr <- df_hcr %>% 
  left_join(ref_points, by = "Code") 
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
# the projection needs to use tuned mFC values hitting mean F from 2016-2020
# the run that has those properties is 2241
# the mFC in 2238 (and the PRM currently being modified) is the base value applied in the burn-in
# the mFC in 2241 is the target mFC to apply in proj
# set up a multiplier for the projection
# set all the necessary flags and switches for the mFC change

# TODO: all this will have to be tested with 1095 (3 years before it all kick in)
# I tested all this before and it should work but more testing will not hurt
file.copy("AtlantisGOAV0_02241/GOA_harvest_background.prm", "AtlantisGOA_MS/GOA_harvest_2241.prm")
harvest_target <- readLines("AtlantisGOA_MS/GOA_harvest_2241.prm")

# get mFC from current template PRM (the 1990s mFC to scale to proj values)
burnin_mfc <- list()
for(sp in codes) {
  
  mfc_line <- harvest_base[grep(paste0("mFC_", sp, " "), harvest_base) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  burnin_mfc[[sp]] <- data.frame("Code" = sp, "burnin_mfc"=  mfc)
}
burnin_mfc <- bind_rows(burnin_mfc)

# get mFC from target harvest (the 2016-2020 mFC)
# NB: Halibut is going to need a different approach. The target mFC will be the mFC associated with the FMSY run
# TODO: plug in HAL values when the SS runs are done
proj_mfc <- list()
for(sp in codes) {
  
  mfc_line <- harvest_target[grep(paste0("mFC_", sp, " "), harvest_target) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  proj_mfc[[sp]] <- data.frame("Code" = sp, "proj_mfc"=  mfc)
}
proj_mfc <- bind_rows(proj_mfc)

# get scalars
scalars <- burnin_mfc %>% 
  left_join(proj_mfc, by = "Code") %>%
  mutate(scalar = proj_mfc / burnin_mfc) %>%
  mutate(scalar = ifelse(is.nan(scalar), 1, scalar))

# these scalars are slightly different than those one obtains from doing Frecent/F1990s
# That is because after tuning 1990s mFC I did some calibration of some groups' biology (recruitment and mortality)
# that biological calibration changed productivity, so now mFC is slightly off
# most stocks are in the right ballpark still (except SKO) - good enough for projecting fixed F for 85 years

# write scalars to mFC change in the harvest.prm template and make all other relevant changes to activate mFC change
for (i in 1:length(oy_grp)){
  
  sp <- oy_grp[i]
  
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
cap_df <- data.frame(cap = rep(c(600, 400, 200), each = 3),
                     weights = rep(c("equal","binary","ramp"), 3)) %>%
  mutate(filename = paste0("GOA_harvest_", cap, weights, ".prm"))

# weights
w_from_specs <- read.csv("data/w_from_harvest_specs.csv")
w_from_specs <- w_from_specs %>% select(Code, w)
# DFD is missing from these weights but part of the OY
# DFD includes grenadiers, but these do not show up in the harvest specs data, nor in the summarized catch accounting data
# according to the latest assessment, they have an estimated exploiation rate corresponding to F = 0.017
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
  
  syst_cap_calc_method_line <- grep("syst_cap_calc_method", this_harvest)
  this_harvest[syst_cap_calc_method_line] <- gsub("syst_cap_calc_method 1", "syst_cap_calc_method 0", this_harvest[syst_cap_calc_method_line])
  
  Ecosystm_Cap_tonnes_line <- grep("Ecosystm_Cap_tonnes", this_harvest)
  this_harvest[Ecosystm_Cap_tonnes_line] <- paste("Ecosystm_Cap_tonnes", cap)
  
  flagSSBforHCR_line <- grep("flagSSBforHCR", this_harvest)
  this_harvest[flagSSBforHCR_line] <- gsub("0", "1", this_harvest[flagSSBforHCR_line])
  
  writeLines(this_harvest, paste0("AtlantisGOA_MS/", filename))
}

# Force.prm ---------------------------------------------------------------




# sh files ----------------------------------------------------------------


