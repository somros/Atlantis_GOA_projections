# 06/25/2025
# TODO: major issue at the moment: all these analyses are model-wide. Should we be limiting these to US?
# Arguments for no: 
# 1. it is one stock; estBo is computed model-wide because it is a population property. It could be computed for the relevant boxes only, 
# but the model will apply it across all boxes - HCRs are not spatial, nor is the OY cap calculation - it is based on whole model biomass.
# Settign F=0 in BC boxes won't help - management calcs will still be carried out across biomass for the whole model domain
# Arguments for yes: FMSY / proxies will be compared to Alaska; OY cap is for Alaska

# as always, BC is a big issue...
# we'd need to at least know the proportion of biomass and catch that is coming from BC, and put that in the discussion
# there is no way to separate the HCR / OY management of the two areas so all calculations depending on biomass will always depend on total biomass
# plotting only alaska biomass in the plots may make sense only for catch vs cap. any panel that incorporates reference point B0, which is used in the HCR,
# cannot be rescaled. Biomass for non-HCR groups may be rescaled, but then we become inconsistent
# let it be a separate function / analysis

# plots that could be rescaled / use information from the US boxes only:
# total biomass and catch - will need biomass from spatial txt and catch from netcdf file
# catch-biomass tradeoff - as above
# H plot - use only boxes in AK

# plots that CANNOT be rescaled because they incorporate B0
# pol-atf tradeoff

# plots that are somewhere in between - they have some info that is based on B0
# ecosystem index plot (we barely look at that)
# by species lineplots: biom_frac has B0 baked into it; oy_rescale column refers to the whole model - then again there are no spatial diffs in the rescaling so it applies to AK-only portion
# diets: technically these refer to average across the model domain. Extreme spatial patterns in BC would make this plot spurious
# HCR: maybe least of all, as they have B/B0 info that is used directly in the HCR calculations

# The thing is, if B0 is the problem, we could certainly compute AK-specific B0 for the sake of the plots
# BUT: doing so would be misleading, would it not? The mgmt routines use a model-wide B0 (they can't do anything other than that)
# And it goes beyond B0. The OY rescaling is operated at model-level. So suppose you have a huge biomass of species X in BC. its predicted catch would be large and affect the rescaling. 

# Maybe this is an opportunity to steer away from specific numbers after all
# Step 1 is running a separate analysis that tells us, in every run we do, biomass and catch proportions between Ak and BC for all OY species

# let's do it here.

# this script explores the effect of biomass and catch from the Atlantis boxes in BC on the application of the HCRs and the OY rescaling
# the script will need some objects created in "analysys.R"

run <- c(2097:2112, 2124:2143)

# biomass: use the BoxBiomass.txt files for this
get_AK_biom <- function(this_run){
  
  print(this_run)
  
  # File paths
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testBoxBiomass.txt")
  
  # Read files once
  biom <- data.table::fread(paste(wd, biom_file, sep = "/"), 
                            sep = " ", 
                            header = TRUE,
                            nThread = getDTthreads())
  
  # subset time early on if necessary
  time_filter <- seq(1, yr_end) * 365
  biom <- biom[Time %in% time_filter & Time/365 <= yr_end]
  
  # Process biomass data once, outside the loop
  # Convert to long format once and use faster string splitting
  biom_long <- data.table::melt(biom, 
                                id.vars = c("Time", "Box"), 
                                variable.name = "Code", 
                                value.name = "mt")
  
  # Pre-filter to only relevant species and times
  biom_filtered <- biom_long[Code %in% oy_species]
  
  # add space
  biom_filtered[, region := fifelse(Box < 92, "AK", "BC")]
  
  # summarize
  biom_space <- biom_filtered[, .(mt = sum(mt)), by = .(Time, region, Code)]
  biom_space[, tot_mt := sum(mt), by = .(Time, Code)]
  biom_space[, prop := mt / tot_mt]
  
  # Filter to AK only and select relevant columns
  result <- biom_space[region == "AK", .(Time, Code, prop, tot_mt)]
  result[, run := this_run]
  
  return(result)
  
}

biom_spatial <- bind_rows(lapply(run, get_AK_biom))

# bring in run info
biom_spatial2 <- biom_spatial %>%
  left_join(key_config, by = "run")

# summary of proportions of biomass in AK
props <- biom_spatial2 %>%
  group_by(Code) %>%
  summarise(mean_prop = mean(prop)) %>%
  arrange(-mean_prop)

# summary of total biomass contribution 
biom_spatial2 %>%
  slice_min(Time) %>%
  slice_min(run) %>%
  mutate(all_biom = sum(tot_mt)) %>%
  mutate(totprop = tot_mt / all_biom) %>%
  select(Code, totprop) %>%
  arrange(-totprop)

# plot
biom_spatial2 %>%
  mutate(Code = factor(Code, levels = sort(oy_species))) %>%
  ggplot(aes(x = Time/365, y = prop, color = factor(run)))+
  geom_line()+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  scale_y_continuous(limits = c(0,1))+
  facet_wrap(~Code)

# For biomass:
# The good news is that the proportion of catch in AK is very constant over time
# Also good news is that proportions stay identical across all runs short of some blips - meaning that management indeed acts homogeneously it seems
# Also good news is that POL and COD are almost entirely in AK (how correct that is is a different question)
# The lowest is RFD with 56% in AK. Here are the props:
# A tibble: 19 × 2
# Code  mean_prop
# <fct>     <dbl>
# 1 SCU       0.996
# 2 POL       0.995
# 3 DFD       0.989
# 4 COD       0.966
# 5 SKO       0.964
# 6 SKB       0.961
# 7 FHS       0.957
# 8 RFP       0.927
# 9 FFS       0.918
# 10 ATF       0.878
# 11 SKL       0.859
# 12 RFS       0.827
# 13 REX       0.824
# 14 SBF       0.810
# 15 THO       0.775
# 16 POP       0.772
# 17 FFD       0.747
# 18 DFS       0.724
# 19 RFD       0.565

# So, where does this leave us? 
# Biomass-based mgmt computations are being done model-wide. estBo provided to the model must be model-wide
# This is most relevant to HCR and OY, but how are they affected?

# HCR management is based on total biomass of a species / B0
# if a species has a very high biomass in BC, F setting will account for that
# BUT: F will also be applied equally across the model domain
# And because proportions remain constant over time (for now), I think that the F applied to AK is still correct';
# Suppose that B/B0 was for AK only - if proportions are constant between the two areas, then B/B0 in AK = B/B0 in BC = B/B0 model-wide
# And F is applied equally in AK and BC as proportion of such biomass)
# So, in fact, the HCR management is not affected because it acts proportionally.
# This will all fall apart with spatially-explicit fleets

# OY management: here is where we have a problem. Each year, projected catch is computed based on total biomass
# if a stock has 30% of its biomass in BC, that 30% will go towards computing projected catch model-wide, which will influence total catch against the cap.
# In other words, if AK only was involved, the cap may not have been met and the OY rescaling, as such, may not have acted
# This is particularly problematic for ATF, which is a large biomass and as such a large catch. 12% of its biomass is in BC.
# Other stocks that have a problem with this include POP and SBF, especially if we want to increase their biomass

# # here is the biomass over total biomass in the model by species, at the BEGINNING of the runs:
# Code     totprop
# <fctr>       <num>
# 1:    POL 0.307396645
# 2:    ATF 0.233660383
# 3:    COD 0.126729313
# 4:    DFS 0.074818847
# 5:    FFS 0.051900526
# 6:    RFS 0.039649252
# 7:    FHS 0.035580129
# 8:    SBF 0.030102056
# 9:    DFD 0.019971701
# 10:    POP 0.019869018
# 11:    FFD 0.017187748
# 12:    REX 0.012487293
# 13:    RFD 0.008278762
# 14:    THO 0.006766575
# 15:    RFP 0.005256663
# 16:    SKB 0.003943915
# 17:    SCU 0.003032599
# 18:    SKL 0.001927199
# 19:    SKO 0.001441374

# so, over a run, what is the proportion of total FMP groundfish biomass that ends up in BC?

spillover <- biom_spatial2 %>%
  left_join(props, by = "Code") %>%
  mutate(ak_mt = tot_mt * mean_prop) %>%
  group_by(Time,run) %>%
  summarise(all_biom = sum(tot_mt), # total model-wide biomass
         all_biom_ak = sum(ak_mt)) %>% # total ak biomass
  ungroup() %>%
  mutate(akprop = all_biom_ak / all_biom)

spillover %>%
  ggplot(aes(x = Time, y = akprop, color = factor(run)))+
  geom_line()

# 84% is the minimum proportion of biomass coming from Alaska at the end of a run
# 91% is the maximum
# the reason why this occurs is that some AK-based stocks (POL and COD) decline more than "BC-based" ATF, and SBF in fact increases
# so the make up of AK-based and BC-based stocks change
# I use BC-based liberally here, all stocks have a majority of their biomass in the AK portion
# This can have an effect on whether the cap-based rescaling is engaged or not
# for example, the "excess" catch from the biomass in BC is going to be considered against the cap

# Catch -------------------------------------------------------------------

# look at the other half of the story
# use the TOTCATCH.nc files as we do not need catch by age class

#' Reads bio and catch nc files for a specific run
#'
#' @param catch_nc_tot path to the CATCHTOT.nc NetCDF file
#' @param run number of calibration run
#'
#' @description 
#' 1. Reads CATCHTOT.nc files for a specific run
#' 2. Extracts catch output per cell (in t)
#' @return data frame with catch information in space per species
#' @export
#' 
#' 
#' 

build_catch_output_TOT <- function(this_run){
  
  print(paste(this_run, "NAA"))
  
  # File paths
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  ncfile <- paste0(wd, "/outputGOA0", this_run, "_testTOTCATCH.nc")
  
  this_tidync <- tidync(ncfile)
  this_nc <- ncdf4::nc_open(ncfile)
  
  catch_nc_tot_ls <- list()
  
  for(i in 1:length(oy_names)){
    fg <- oy_names[i] # this needs to use the "Name" to pull from the NC file
    code <- grps %>% filter(Name == fg) %>% pull(Code)
    
    #Extract from the output .nc file the appropriate catch time series variables
    catch_vars <- this_tidync %>%
      activate("D1,D0") %>%
      hyper_vars() %>% # all variables in the .nc file active grid
      filter(grepl("Tot_.*_Catch",name)) %>% # filter for reserve N
      filter(grepl(code,name)) # filter for specific functional group
    
    catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this_nc) 
    
    catch_df <- as.data.frame((catch))
    
    # add box_id
    catch_df <- catch_df %>% mutate(box_id = 0:108)
    
    # reshape
    catch_df_long <- catch_df %>%
      pivot_longer(-box_id, names_to = "ts", values_to = "mt")
    
    # turn ts column to integer
    catch_df_long <- catch_df_long %>%
      mutate(ts = gsub("X","",ts)) %>%
      mutate(ts = as.numeric(ts)) %>%
      mutate(ts = ts - 1) # start numbering ts from 0
    
    catch_box_df <- catch_df_long
    
    # drop ts = 0
    catch_box_df <- catch_box_df %>%
      filter(ts > 0)
    
    # add species
    catch_box_df <- catch_box_df %>%
      mutate(Name = fg,
             Code = code)
    
    catch_nc_tot_ls[[i]] <- catch_box_df
    
  }
  
  catch_nc_TOT <- bind_rows(catch_nc_tot_ls)
  
  catch_nc_TOT <- catch_nc_TOT %>% mutate(run = this_run)
  
  return(catch_nc_TOT)
}

# extract catch by box from all runs
catch_spatial <- bind_rows(lapply(run, build_catch_output_TOT))

# attribute to AK or BC
catch_spatial2 <- catch_spatial %>%
  mutate(region = ifelse(box_id<92,"AK","BC")) %>%
  group_by(ts,Name,Code,run,region) %>% # sum over boxes
  summarise(mt_region = sum(mt)) %>%
  group_by(ts,Name,Code,run) %>% # tot biomass
  mutate(mt_all = sum(mt_region)) %>%
  ungroup() %>%
  mutate(prop = mt_region / mt_all) %>%
  filter(region == "AK")

# summary of total biomass contribution 
tt <- catch_spatial2 %>%
  filter(ts == 15) %>%
  #filter(ts == 100) %>% # you can't do ts=1 because HCR mgmt has not kicked in yet
  group_by(ts, run) %>%
  mutate(all_catch = sum(mt_all)) %>%
  mutate(totprop = mt_all / all_catch) %>%
  select(Code, totprop, run) %>%
  arrange(-totprop) %>%
  filter(Code == "ATF")

# the proportion of ATF's catch is 21%, towards the end of run 2111, a binary ssp585
# a run like this will suffer the most from the issue of having biomass / catch in BC
# the 200K cap is still constraining here, and it would not be if we were using AK only
# at the end of the burn-in, ATF catch is <10% of the total

# plot ratios by species 
catch_spatial2 %>%
  mutate(Code = factor(Code, levels = sort(oy_species))) %>%
  ggplot(aes(x = ts, y = prop, color = factor(run)))+
  geom_line()+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  scale_y_continuous(limits = c(0,1))+
  facet_wrap(~Code)

# this plot is fundamentally identical to the biomass one, which is no surprise because mFC is not spatially explicit
# worth figuring out the contributions to total catch and how much of that is in BC, though


spillover_catch <- biom_spatial2 %>%
  left_join(props, by = "Code") %>%
  mutate(ak_mt = tot_mt * mean_prop) %>%
  group_by(Time,run) %>%
  summarise(all_catch = sum(tot_mt), # total model-wide catch
            all_catch_ak = sum(ak_mt)) %>% # total ak catch
  ungroup() %>%
  mutate(akprop = all_catch_ak / all_catch)

spillover_catch %>%
  ggplot(aes(x = Time, y = akprop, color = factor(run)))+
  geom_line()

# proportion of catch in BC is the same as the proportion of biomass in BC
# does this make sense? Yes, because whichever F rate is used for a group is going to be consistent in space

# in conclusion: 

# As always we have a problem with the BC portion of the model
# The problem is that we think of Alaska when we do the total catch and biomass calculations with respect to the caps, but the model has boxes in BC
# This is not just a visualization problem, which all in all would be of easy solution
# More deeply, the OY algorithm uses model-level catches to rescale catch under the OY cap
# So, a stock that has a lot of biomass and catch coming from BC may be important in its contribution to the cap
# The good news is that POL and COD are mostly in AK, and those are the most important stocks in terms of biomass and catch
# it's also good that the proportions are constant over time and across runs
# also good that this does not really affect the HCR because those calclations are all proportional
# it does affect the OY because those calculations are with raw catch
# so ignoring the smaller stocks, this may become a problem for arrowtooth.

# The threns are not affected from what I can see

# options for when we write this up are:
# 1. simply state that the cap is applied to the whole model domain. I can live with this, 
# it illustrates the concepts nonetheless, it is cleaner to explain. Drawback: there is no OY management in BC, and this will be difficult to go around.
# people will wonder why we would apply an ecosystem cap across jurisditions, it will never happen and it will muddle the waters in terms of this being useful to the NPFMC
# 2. amend the plots so that they only capture the AK portion of catch and biomass, including rescaling B0. Issue with this is
# that the cap hlines will refer to the entire model domain, so your total catch plots will be off


