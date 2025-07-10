# Alberto Rovellini
# 6/27/2025
# This script produces mFC values for Atlantis GOA projections based on historical F values from GOA assessments
# F time series were extracted from the most recent FULL assessment for each stock as of May 2025
# https://docs.google.com/spreadsheets/d/1RpWddJOC8l4aPhS85SQ615Qd3qMUhs6v/edit?usp=drive_link&ouid=105160877592505917667&rtpof=true&sd=true

# Some decision points will include:
# What to do with exploitation rates
# What time period we should use (1990's vs last n years before 2020 vs last n years before 2025)
# This last point opens up a whole different conversation
# For the purpose of getting the reference poinbts with the F curves, we should use 1990-2000 rates
# This is because those are the biomasses we calibrated the model to
# The burn-in could use those values, then switch to 2015-2020 values for projections
# Start projecting from 2020 because that is the end of the ROMS hindcast

# maybe the way to go is: 
# 1990's F for reference point calculations
# 30-yr burn-in with 1990's climatology
# tune initial biomasses so that they are closer to 2015-2020 averages than 1990
# start projections from 2020 
# the goal ultimately is to be in the right ballpark in 2020
# this is particularly important for COD, POP, SBF

library(readxl)
library(tidyverse)

# read in species info
grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)

# FMP groundfish ----------------------------------------------------------

dat <- read_xlsx("data/F_from_assessments.xlsx", range = "A1:V49")

# pivot
dat_long <- dat %>%
  pivot_longer(-Year, names_to = "Code_tmp", values_to = "F")

# handle duplicate names:
# New names:
# • `FFS` -> `FFS...9`
# • `FFS` -> `FFS...10`
# • `RFS` -> `RFS...15`
# • `RFS` -> `RFS...16`

key <- data.frame("Code_tmp" = c("FFS...9", "FFS...10", "RFS...15", "RFS...16"),
                  "Species" = c("NRS", "SRS", "Northern", "REBS"),
                  "Code" = c("FFS","FFS","RFS","RFS"),
                  "Biomass_1990" = c(42569,96164,215131,43492))

dat_long1 <- dat_long %>% filter(!Code_tmp %in% key$Code_tmp) %>% rename(Code = Code_tmp)
dat_long2 <- dat_long %>% filter(Code_tmp %in% key$Code_tmp)

dat_long2 <- dat_long2 %>%
  left_join(key, by = "Code_tmp") %>%
  group_by(Year,Code) %>%
  summarise(`F` = weighted.mean(`F`,Biomass_1990))

dat_final <- rbind(dat_long1, dat_long2) %>%
  arrange(Year,Code)

# TODO: handle exploitation rates where needed

dat_1990s <- dat_final %>%
  filter(Year %in% c(1990:1999)) %>%
  group_by(Code) %>%
  summarise(`F` = mean(`F`, na.rm = T))

# turn to mFC, but only where needed
exp_rate_stocks <- c("SKL","SKO","SKB","RFD","THO","DFD","SCU")

dat_1990s <- dat_1990s %>%
  rowwise() %>%
  mutate(mFC = ifelse(Code %in% exp_rate_stocks,`F`/365,1-exp(-`F`/365))) %>%
  ungroup()

# Other groups ------------------------------------------------------------
# Lots of groups to fill in. The only groups without estimates for 1990s but with estimates for later decades are skates
# skates have exploitation rates
# I'd say use old 1/4 M values for consistency...

# Pacific halibut
# the target is SPR 46%
# we need to convert that to biomass as per the assessment and get the equivalent depletion rference point
# we then need to profile over F and get the F that leads us there
# the old model used a low F for HAL, 1/4 M
# F in the new model will be much higher, and as such we will need to calibrate HAL most likely
# to get the model closer to where it will be, use FMSY from the OY paper
# FMSY for HAL was 0.0115, which corresponded to 44% depletion
# this is probably in the ballpark of F46%, so use that for calibration purposes
f_hal <- 0.0115
mfc_hal <- 1-exp(-f_hal/365)

dat_1990s <- rbind(dat_1990s,
                   data.frame("Code"="HAL",
                              "F" = f_hal,
                              "mFC" = mfc_hal))

# For all others use the old values (for now)

# Compare to F values used in the OY runs ---------------------------------

# you can only compare the groups that have avalue
dat_1990s <- dat_1990s %>% filter(!is.nan(`F`))

# read recent harvest.prm
harvest_file <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_2097/GOA_harvest_background.prm"
harvest <- readLines(harvest_file)

# get mFC values
prm_df <- data.frame()
for(i in seq_along(unique(dat_1990s$Code))){
  
  sp <- dat_1990s$Code[i]
  
  mfc_line <- harvest[grep(paste0("mFC_", sp, " "), harvest) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  prm_df <- rbind(prm_df, data.frame(sp, mfc))
  
}

colnames(prm_df) <- c("Code","mFC_prm")

# compare for groups that have new values
dat_1990s %>%
  left_join(prm_df) %>%
  mutate(prop = mFC / mFC_prm) %>%
  arrange(-prop)

# there are some differences here and we should test this ASAP - may take some calibration
# POL and COD can take it based on the HCR? POL is fished at FMSY and barely goes below B40
# worried that SBF will collapse and ATF will boom
# neat that ATF will be even more inconsequential for BC caluculations

# while we are at it and doing a run, let's do another with 2010s values
dat_2010s <- dat_final %>%
  filter(Year %in% c(2010:2019)) %>%
  group_by(Code) %>%
  summarise(`F` = mean(`F`, na.rm = T))

# turn to mFC
dat_2010s <- dat_2010s %>%
  rowwise() %>%
  mutate(mFC = ifelse(Code %in% exp_rate_stocks,`F`/365,1-exp(-`F`/365))) %>%
  ungroup()

# compare for groups that have new values
dat_2010s %>%
  filter(!is.nan(`F`)) %>%
  left_join(prm_df) %>%
  mutate(prop = mFC / mFC_prm) %>%
  arrange(-prop)

# would be very rough for Pcod

# if stocks suffer too much from these you can always ease on mL and mQ (which would be a godsend)
# halibut will be a key example
# for stocks exploited more lightly we may have the opposite problem and have to introduce more mL or mQ

# NB: these values have to be calibrated

# write mfc file

# neet to get a recent mFC for all the values that are NOT being changed to 1990s
oldfile <- "C:/Users/Alberto Rovellini/Documents/GOA/Atlantis_GOA_OY_MS/GOA_harvest_background.prm"
ref_harvest <- readLines(oldfile)

# get harvest parameters outside loop
ref_mfc <- list()
for(sp in codes) {
  
  mfc_line <- ref_harvest[grep(paste0("mFC_", sp, " "), ref_harvest) + 1]
  mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
  
  ref_mfc[[sp]] <- list(mfc = mfc)
}

newfile <- paste0('mFC_tuning/mfc_projections_BEFORE_TUNING.prm')
file.create(newfile)

for(i in 1:length(codes)){
  
  sp <- codes[i]
  parname <- paste0("mFC_", sp, " 33")
  
  print(sp)
  
  if(sp %in% dat_1990s$Code){
    parval <- dat_1990s %>% filter(Code == sp) %>% pull(mFC)
  } else {
    parval <- ref_mfc[[which(names(ref_mfc)==sp)]]$mfc
  }
  
  parval <- signif(parval, 8)
  
  parline <- paste(as.character(c(parval,rep(0,32))),collapse = " ")
  
  cat(parname, file=newfile, append=TRUE,'\n')
  cat(parline, file=newfile, append=TRUE, '\n')
  
}
