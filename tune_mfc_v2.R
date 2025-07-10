# Alberto Rovellini
# 6/28/2025
# tune mFC values based on realized catch/selexbiomass
# F values in the assessments are fully selected, so we need to use selected biomass from Atlantis

library(tidyverse)
rm(list = ls())

# read in species info
grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)
fished <- grps %>% filter(IsFished==1) %>% pull(Code)

# Reference run with expected (target) mfc --------------------------------
# reference PRM with the target mFC
ref_run <- 2179
ref_wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", ref_run)
ref_harvest_prm <- list.files(ref_wd)[grep("GOA_harvest_.*.prm", list.files(ref_wd))]
ref_harvest <- readLines(paste(ref_wd, ref_harvest_prm, sep = "/"))

# get harvest parameters outside loop
ref_mfc <- list()
for(sp in fished) {
  
  mfc_line <- ref_harvest[grep(paste0("mFC_", sp, " "), ref_harvest) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  ref_mfc[[sp]] <- list(mfc = mfc)
}

# turn mfc to f
f_expected <- do.call(rbind, lapply(names(ref_mfc), function(x) {
  data.frame(Code = x, mfc = ref_mfc[[x]]$mfc)
})) %>%
  mutate(f_e = -(365)*log(1-mfc)) %>%
  select(Code,f_e,mfc)

# Runt to tune ------------------------------------------------------------
# what run are you tuning
this_run <- 2181

# prepare file paths
wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]

# Read files
biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
harvest <- readLines(paste(wd, harvest_prm, sep = "/"))

# get harvest parameters outside loop
harvest_params <- list()
for(sp in fished) {
  
  startage_line <- harvest[grep(paste0(sp, "_mFC_startage"), harvest) + 1]
  startage <- as.numeric(strsplit(startage_line, split = " ")[[1]])[1]
  
  mfc_line <- harvest[grep(paste0("mFC_", sp, " "), harvest) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  harvest_params[[sp]] <- list(startage = startage, mfc = mfc)
}

# # store the mfc values from the run to tune
mfc_to_tune <- do.call(rbind, lapply(names(harvest_params), function(x) {
  data.frame(Code = x, mfc = harvest_params[[x]]$mfc)
}))

# Process biomass data
biom <- biom %>% filter(Time==0) # biomass at t0

# reshape
biom_long <- biom %>% pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") # Convert to long format

# handle ages and names
code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
biom_long$Code <- sapply(code_age_split, `[`, 1)
biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
biom_long <- biom_long %>% select(-Code.Age)  # Remove original column

# get selected biomass at t0
biom_selex <- biom_long %>%
  left_join(
    data.frame(Code = names(harvest_params),
               startage = sapply(harvest_params, `[[`, "startage")),
    by = "Code"
  ) %>%
  filter(Age >= startage) %>%
  group_by(Time, Code) %>%
  summarise(biom_mt_selex = sum(mt), .groups = 'drop') %>%
  select(-Time)

# process catch data
catch <- catch %>% filter(Time==365) # catch at t1 - note the staggered time step
catch_long <- catch %>%
  dplyr::select(Time:BIV) %>%
  pivot_longer(-Time, names_to = "Code", values_to = "catch") %>%
  select(-Time)

# get annual realized exploitation rate
f_realized <- biom_selex %>%
  left_join(catch_long, by = "Code") %>%
  mutate(exp = catch / biom_mt_selex,
         f_r = -log(1-exp)) %>%
  select(Code,f_r)

# compare realized to expected
correction_factors <- f_realized %>%
  left_join(f_expected, by = "Code") %>%
  mutate(prop = f_e / f_r) %>%
  select(Code, prop)

# if NaN, turn to 1
correction_factors[is.nan(correction_factors$prop),]$prop <- 1

# as noted before, typically those the need the most correction are the migrating salmon bc they leave the model
# everything else is between 0.7 and 1.25 of where it should be so doing OK
# probably just need a couple rounds of tuning
# SPI ends up being very high though...

# bump up mfc values by prop
new_mfc <- mfc_to_tune %>%
  left_join(correction_factors, by = "Code") %>%
  mutate(new_mfc = mfc * prop)

# write new mfc vectors

# turn mFC_frame_new back to section of prm file and write out
newfile <- paste0('mFC_tuning/mFC_from_', this_run, '.prm')
file.create(newfile)

for(i in 1:length(codes)){
  
  sp <- codes[i]
  parname <- paste0("mFC_", sp, " 33")
  
  parval <- new_mfc %>% filter(Code == sp) %>% pull(new_mfc)
  if(length(parval)==0){
    parval <- 0
  }
  
  if(sp == "SPI"){ # do not tune SPI
    parval <- harvest_params[[which(names(harvest_params)==sp)]]$mfc
  }
  
  parline <- paste(as.character(c(parval,rep(0,32))),collapse = " ")
  
  cat(parname, file=newfile, append=TRUE,'\n')
  cat(parline, file=newfile, append=TRUE, '\n')
  
}

