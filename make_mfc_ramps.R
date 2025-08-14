# Alberto Rovellini
# 08/07/2025
# Code to create mFC ramps for reference point simulations
library(tidyverse)

# read in groups
grps <- read.csv("data/GOA_Groups.csv")

# start from reference run's PRM - using run 2238 as best start
# read in reference PRM file
ref_prm_file <- "AtlantisGOAV0_02238/GOA_harvest_background.prm"
ref_prm <- readLines(ref_prm_file)

# HCR stocks and Pacific halibut ------------------------------------------
# For these stocks we need to run F profiles to identify FMSY
# For the HCR stocks (POL,COD,POP,SBF), FMSY will be used as Ftarget in the HCR
# For HAL, we assume fishing at FMSY as proxy of F46%
# this is akin to what was done in Rovellini et al. (2025), except instead of scaling FOFL we scale the Atlantis FMSY
# We also allow for some more resolution around the peak while the tails are coarser
# the profile must start from F=0 so that we capture B0
# go high, F=4*FMSY, because it is interesting to see the descending limb of the curve

# read in FMSY vecs from Rovellini et al. (2025)
fmsy_atlantis <- read.csv("data/fmsy_Rovellini_2025.csv")

# filter 5 stocks of interest - the stocks for which we need the full ramp to get FMSY
hcr_stocks <- grps %>% filter(Code %in% c("POL","COD","POP","SBF","HAL")) %>% select(Code,LongName)

fmsy_atlantis <- fmsy_atlantis %>%
  filter(LongName %in% hcr_stocks$LongName) %>%
  select(LongName, atlantis_fmsy)

# turn to mFC corresponding (approximately) to FMSY
# mFC =1-exp(-F/365)
fmsy_atlantis <- fmsy_atlantis %>%
  mutate(mfc_at_fmsy = 1 - exp(-atlantis_fmsy / 365)) %>%
  left_join(hcr_stocks)

# make ramp, using multipliers: c(0, 0.5, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2, 3, 4)
ramp <- c(0, 0.25, 0.5, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2, 3, 4)

# create a folder to store the harvest.prm and the run_atlantis.sh files
# if files already exist the file writing below will get corrupted so purge old versions first
outdir <- "output/single_species_runs/mfc_ramps/"
if(dir.exists(outdir)){
  unlink(outdir, recursive = T)
}
dir.create(outdir, recursive = T)

sim_idx <- 0 # start the counter, it will be used to label the runs
ref_table <- data.frame()
for(i in 1:length(hcr_stocks$Code)){
  
  sp <- hcr_stocks$Code[i]
  fmsy <- fmsy_atlantis %>% filter(Code == sp) %>% pull(atlantis_fmsy)
  mfc_at_fmsy <- fmsy_atlantis %>% filter(Code == sp) %>% pull(mfc_at_fmsy)
  
  for(j in 1:length(ramp)){
    
    mult <- ramp[j]
    
    # what is the mfc in the burnin period?
    mfc_line <- ref_prm[grep(paste0("mFC_", sp, " "), ref_prm) + 1]
    burnin_mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
    
    # what is the target mfc?
    target_mfc <- mfc_at_fmsy * mult
    
    # what is the mfc_change factor to apply to reach the target mfc?
    mfc_change <- target_mfc / burnin_mfc
    
    # write all relevant lines in the PRM file
    # create a new PRM file copying the base file
    idx <- sprintf("%03d", sim_idx)
    new_prm_file <- paste0("GOA_harvest_", idx, ".prm")
    file.copy(ref_prm_file, paste(outdir, new_prm_file, sep = "/"))
    
    # open file to modify
    new_prm <- readLines(paste(outdir, new_prm_file, sep = "/"))
    
    # make changes to PRM lines
    # flagchangeF 1
    old_line <- new_prm[grep("flagchangeF", new_prm)]
    new_line <- gsub("\t0", "\t1", old_line)
    new_prm[grep("flagchangeF", new_prm)] <- new_line
    
    # flagFchange_XXX 1 for first element (which is BG and needs to go off)
    old_line <- new_prm[grep(paste0("flagFchange_", sp), new_prm) + 1]
    new_line <- sub("0", "1", old_line)
    new_prm[grep(paste0("flagFchange_",sp), new_prm) + 1] <- new_line
    
    # XXX_mFC_changes  - how many changes? Set to 1 for first element
    old_line <- new_prm[grep(paste0(sp, "_mFC_changes"), new_prm) + 1]
    new_line <- sub("0", "1", old_line)
    new_prm[grep(paste0(sp, "_mFC_changes"), new_prm) + 1] <- new_line
    
    # mFCchange_start_XXX set to 30*365
    old_line <- new_prm[grep(paste0("mFCchange_start_", sp), new_prm) + 1]
    change_start <- 365 * 30
    new_line <- sub("0", as.character(change_start), old_line)
    new_prm[grep(paste0("mFCchange_start_", sp), new_prm) + 1] <- new_line
    
    # mFCchange_period_XXX (1 day)
    old_line <- new_prm[grep(paste0("mFCchange_period_", sp), new_prm) + 1]
    new_line <- sub("0", "1", old_line)
    new_prm[grep(paste0("mFCchange_period_", sp), new_prm) + 1] <- new_line
    
    # mFCchange_mult_XXX set to mfc_change
    old_line <- new_prm[grep(paste0("mFCchange_mult_", sp), new_prm) + 1]
    new_line <- as.character(round(mfc_change, digits = 5))
    new_prm[grep(paste0("mFCchange_mult_", sp), new_prm) + 1] <- new_line
    
    # write out modified prm file
    writeLines(new_prm, paste0(outdir, new_prm_file))
    
    # create a new SH file
    this_sh_file <- paste0("run_atlantis_", idx, ".sh")
    this_sh <- paste0("atlantisMerged -i GOA_cb_summer.nc  0 -o outputGOA_",
                      idx,
                      ".nc -r GOA_run.prm -f GOA_force.prm -p GOA_physics.prm -b GOAbioparam.prm -h GOA_harvest_",
                      idx,
                      ".prm -m GOAMigrations.csv -s GOA_Groups.csv -q GOA_fisheries.csv -d outputFolder",
                      idx)
    writeLines(this_sh, paste(outdir, this_sh_file, sep = "/"))
    
    # fill in the reference table
    ref_table <- rbind(ref_table, data.frame(
      "idx" = idx,
      "Code" = sp,
      "fmsy" = fmsy,
      "fmsy_mfc" = mfc_at_fmsy,
      "mult" = mult,
      "targ_mfc" = target_mfc
    ))
    
    # update sim index for next iteration
    sim_idx <- sim_idx + 1
    
    # done
  }
}

# All other FMP stocks ----------------------------------------------------
# For all other groundfish stocks we only need F=0. We need B0 for plotting purposes (biomass fraction)
# we do not need FMSY because the target F in the projections will be recent F - they have no HCR
# these should be 15 stocks:
# leaving sharks and cephalopods out - seem to be negligible catches
oy_groups <- c("ATF", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "RFS", "RFP", "RFD", "THO", "DFS", "DFD", "SCU")

# keep sim_idx as it was 
for (i in 1:length(oy_groups)){
  
  sp <- oy_groups[i]
  mult <- 0
  
  # what is the target mfc?
  target_mfc <- 0
  
  # what is the mfc_change factor to apply to reach the target mfc?
  mfc_change <- 0
  
  # write all relevant lines in the PRM file
  # create a new PRM file copying the base file
  idx <- sprintf("%03d", sim_idx)
  new_prm_file <- paste0("GOA_harvest_", idx, ".prm")
  file.copy(ref_prm_file, paste(outdir, new_prm_file, sep = "/"))
  
  # open file to modify
  new_prm <- readLines(paste(outdir, new_prm_file, sep = "/"))
  
  # make changes to PRM lines
  # flagchangeF 1
  old_line <- new_prm[grep("flagchangeF", new_prm)]
  new_line <- gsub("\t0", "\t1", old_line)
  new_prm[grep("flagchangeF", new_prm)] <- new_line
  
  # flagFchange_XXX 1 for first element (which is BG and needs to go off)
  old_line <- new_prm[grep(paste0("flagFchange_", sp), new_prm) + 1]
  new_line <- sub("0", "1", old_line)
  new_prm[grep(paste0("flagFchange_",sp), new_prm) + 1] <- new_line
  
  # XXX_mFC_changes  - how many changes? Set to 1 for first element
  old_line <- new_prm[grep(paste0(sp, "_mFC_changes"), new_prm) + 1]
  new_line <- sub("0", "1", old_line)
  new_prm[grep(paste0(sp, "_mFC_changes"), new_prm) + 1] <- new_line
  
  # mFCchange_start_XXX set to 30*365
  old_line <- new_prm[grep(paste0("mFCchange_start_", sp), new_prm) + 1]
  change_start <- 365 * 30
  new_line <- sub("0", as.character(change_start), old_line)
  new_prm[grep(paste0("mFCchange_start_", sp), new_prm) + 1] <- new_line
  
  # mFCchange_period_XXX (1 day)
  old_line <- new_prm[grep(paste0("mFCchange_period_", sp), new_prm) + 1]
  new_line <- sub("0", "1", old_line)
  new_prm[grep(paste0("mFCchange_period_", sp), new_prm) + 1] <- new_line
  
  # mFCchange_mult_XXX set to mfc_change
  old_line <- new_prm[grep(paste0("mFCchange_mult_", sp), new_prm) + 1]
  new_line <- as.character(round(mfc_change, digits = 5))
  new_prm[grep(paste0("mFCchange_mult_", sp), new_prm) + 1] <- new_line
  
  # write out modified prm file
  writeLines(new_prm, paste0(outdir, new_prm_file))
  
  # create a new SH file
  this_sh_file <- paste0("run_atlantis_", idx, ".sh")
  this_sh <- paste0("atlantisMerged -i GOA_cb_summer.nc  0 -o outputGOA_",
                    idx,
                    ".nc -r GOA_run.prm -f GOA_force.prm -p GOA_physics.prm -b GOAbioparam.prm -h GOA_harvest_",
                    idx,
                    ".prm -m GOAMigrations.csv -s GOA_Groups.csv -q GOA_fisheries.csv -d outputFolder",
                    idx)
  writeLines(this_sh, paste(outdir, this_sh_file, sep = "/"))
  
  # fill in the reference table
  ref_table <- rbind(ref_table, data.frame(
    "idx" = idx,
    "Code" = sp,
    "fmsy" = NA,
    "fmsy_mfc" = NA,
    "mult" = mult,
    "targ_mfc" = target_mfc
  ))
  
  # update sim index for next iteration
  sim_idx <- sim_idx + 1
  
  # done
}

# write out the reference table
write.csv(ref_table, paste0(outdir, "reference_table.csv"), row.names = F)

# done
