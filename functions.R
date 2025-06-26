# Functions ---------------------------------------------------------------

#' Extract Fishery Information for OY Species
#'
#' @description
#' Extracts and processes fishery-related data for species managed under the 
#' Optimum Yield (OY) framework from Atlantis model output files. Computes
#' biomass, catch, fishing mortality, and ecosystem cap rescaling factors.
#'
#' @param this_run Integer. The run number/ID for the Atlantis model simulation
#'
#' @return A data frame containing:
#'   \item{Time}{Time step from model output}
#'   \item{Code}{Species code}
#'   \item{w}{Weight/preference value for the species}
#'   \item{biom_mt_selex}{Selectable biomass in metric tons}
#'   \item{catch_mt}{Catch in metric tons}
#'   \item{biom_mt_tot}{Total biomass in metric tons}
#'   \item{mu}{Exploitation rate (catch/biomass)}
#'   \item{f}{Instantaneous fishing mortality rate}
#'   \item{biom_frac}{Biomass as fraction of B0}
#'   \item{fref}{Reference fishing mortality}
#'   \item{run}{Run number}
#'   \item{oy_rescale}{OY rescaling factor from ecosystem capacity}
#'
#' @details
#' Reads biomass, catch, harvest parameters, and ecosystem cap files.
#' Processes data for all OY species, computing fishing mortality from catch
#' and biomass data. For HCR-managed species, uses estimated B0 and Fref
#' from harvest parameters. For non-HCR species, computes Fref from maximum
#' fishing rates from Atlantis (mFC).
#'
#' @examples
#' \dontrun{
#' fishery_data <- pull_fishery_info(2097)
#' }
pull_fishery_info <- function(this_run){
  
  print(this_run)
  
  # File paths
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  ecocap_file <- paste0(wd, "/outputGOA0", this_run, "_test_EcosystemCapResult.txt")
  
  # Read files once
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  harvest <- readLines(paste(wd, harvest_prm, sep = "/"))
  ecocap_report <- read.delim(ecocap_file, sep = " ")
  
  # subset time early on if necessary
  biom <- biom %>% filter(Time/365 <= yr_end)
  catch <- catch %>% filter(Time/365 <= yr_end)
  ecocap_report <- ecocap_report %>% filter(Time/365 <= yr_end)
  
  # Process biomass data once, outside the loop
  # Convert to long format once and use faster string splitting
  biom_long <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")
  
  code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
  biom_long$Code <- sapply(code_age_split, `[`, 1)
  biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
  
  # Pre-filter to only relevant species and times
  biom_filtered <- biom_long %>%
    filter(Code %in% oy_species,
           Time > 0, Time < max(Time)) %>%
    select(-Code.Age)  # Remove original column
  
  # Pre-process harvest parameters outside loop
  harvest_params <- list()
  for(sp in oy_species) {
    startage_line <- harvest[grep(paste0(sp, "_mFC_startage"), harvest) + 1]
    startage <- as.numeric(strsplit(startage_line, split = " ")[[1]])[1]
    
    # get weight
    idx <- grep(sp, grps %>% pull(Code))
    w_line <- harvest[grep("SystCapSPpref", harvest) + 1]
    w <- as.numeric((strsplit(w_line, split = " "))[[1]])[idx]
    
    if(sp %in% hcr_spp) {
      sp_idx <- grep(sp, codes)
      estbo_vec <- as.numeric(strsplit(harvest[grep("estBo\t", harvest) + 1], split = " ")[[1]])
      estbo <- estbo_vec[sp_idx]
      fref_vec <- as.numeric(strsplit(harvest[grep("Fref\t", harvest) + 1], split = " ")[[1]])
      fref <- -log(1 - fref_vec[sp_idx])
    } else {
      estbo <- NA
      mfc_line <- harvest[grep(paste0("mFC_", sp, " "), harvest) + 1]
      mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
      fref <- -(365) * log(1 - mfc)
    }
    
    harvest_params[[sp]] <- list(startage = startage, w = w, estbo = estbo, fref = fref)
  }
  
  # Process all species at once using vectorized operations
  # Create biomass summaries for all species
  biom_selex_all <- biom_filtered %>%
    left_join(
      data.frame(Code = names(harvest_params),
                 startage = sapply(harvest_params, `[[`, "startage"),
                 w = sapply(harvest_params, `[[`, "w")),
      by = "Code"
    ) %>%
    filter(Age >= startage) %>%
    group_by(Time, Code, w) %>%
    summarise(biom_mt_selex = sum(mt), .groups = 'drop')
  
  biom_tot_all <- biom_filtered %>%
    group_by(Time, Code) %>%
    summarise(biom_mt_tot = sum(mt), .groups = 'drop')
  
  # Process catch data for all species
  catch_long <- catch %>%
    select(Time, all_of(oy_species)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt")
  
  # join operation
  result_list <- list()
  for(sp in oy_species) {
    sp_params <- harvest_params[[sp]]
    
    sp_data <- biom_selex_all %>%
      filter(Code == sp) %>%
      left_join(filter(catch_long, Code == sp), by = c("Time", "Code")) %>%
      left_join(filter(biom_tot_all, Code == sp), by = c("Time", "Code")) %>%
      mutate(
        mu = catch_mt / biom_mt_selex,
        f = -log(1 - mu),
        biom_frac = biom_mt_tot / sp_params$estbo,
        fref = sp_params$fref,
        run = this_run
      )
    
    result_list[[sp]] <- sp_data
  }
  
  res_df <- do.call(rbind, result_list)
  
  # Process ecocap report
  ecocap_report <- ecocap_report %>%
    select(Time, SpeciesName, FisheryName, SpBased_Frescale, PostSystCap_Frescale) %>%
    filter(SpeciesName %in% oy_species,
           FisheryName %in% oy_fleets) %>%
    mutate(oy_rescale = PostSystCap_Frescale / SpBased_Frescale) %>%
    select(Time, SpeciesName, oy_rescale) %>%
    rename(Code = SpeciesName)
  
  # Final join
  res_df <- res_df %>%
    left_join(ecocap_report, by = c("Time", "Code"))
  
  # TODO: remove this chunk when we add estBo to the harvest.prm file for all stocks
  res_df <- res_df %>%
    left_join(estbo_key, by = "Code") %>%
    mutate(biom_frac = biom_mt_tot / estbo) %>%
    select(-estbo)
  
  return(res_df)
}

#' Generate Fishery-Related Plots
#'
#' @description
#' Creates a comprehensive set of plots for analyzing fishery performance
#' under different management scenarios, including total catch vs caps,
#' species-specific time series, HCR plots, and ecosystem tradeoffs.
#'
#' @param catch_df Data frame. Output from pull_fishery_info() containing
#'   fishery data across multiple runs and scenarios
#'
#' @return NULL (function creates and saves plots to plotdir)
#'
#' @details
#' Generates the following plot types:
#' \itemize{
#'   \item Total catch and biomass vs cap constraints
#'   \item Species-specific time series (catch, F/Ftarg, B/B0, OY rescaling)
#'   \item Harvest Control Rule (HCR) phase plots for managed species
#'   \item Catch vs biomass tradeoffs for preferred species
#'   \item Pollock-Arrowtooth flounder ecosystem tradeoffs
#'   \item Ecosystem indicators summary (dot plots)
#' }
#'
#' Plots are saved to subdirectories: by_spp/, hcr/, and main plotdir.
#' Uses PNW color palette and distinguishes scenarios by cap levels,
#' weight schemes, and climate scenarios.
#'
#' @examples
#' \dontrun{
#' plot_fishery(catch_df)
#' }
plot_fishery <- function(catch_df){
  
  # make a palette
  cap_col <- pnw_palette(name="Sunset2",n=length(unique(catch_df$cap)),type="discrete")
  
  # total catch against cap
  p1 <- catch_df %>%
    filter(Time >= burnin*365) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time,run,cap,wgts,env,other) %>%
    summarise(catch_tot = sum(catch_mt),
              biom_tot = sum(biom_mt_selex)) %>%
    ungroup() %>%
    pivot_longer(-c(Time:other)) %>%
    ggplot(aes(x = Time/365, y = value, color = factor(cap), linetype = factor(wgts)))+
    geom_line(linewidth = 0.5)+
    scale_color_manual(values = cap_col)+
    #scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    scale_linetype_manual(values = c("solid","dashed","dotted"))+
    geom_hline(aes(yintercept = as.numeric(cap)), linetype = "dashed", linewidth = 0.25)+
    theme_bw()+
    scale_y_continuous(limits = c(0,NA))+
    labs(x = "Year", 
         y = "mt",
         color = "Cap (mt)",
         linetype = "Weight scheme",
         title = "Total biomass and catch of OY species")+
    facet_grid(name~env, scales = "free_y")
  
  ggsave(paste0(plotdir, "/", "oy_tot.png"), p1, 
         width = 9, height = 5,  
         units = "in", dpi = 300)
  
  # by species, biom fraction (need B0), f fraction, catch, biomass (selected), exploitation rate, ...?
  # maybe it would be best to produce this plot one species at a time
  # Loop through each code and create a plot
  
  for(i in seq_along(oy_species)) {
    current_code <- oy_species[i]
    current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
    
    p_tmp <- catch_df %>%
      filter(Time >= burnin*365) %>%
      filter(Code == current_code) %>%
      filter(!is.na(catch_mt)) %>%
      mutate(f_frac = f / fref) %>%
      select(-fref, -f, -biom_mt_tot, -biom_mt_selex, -mu, -w) %>%
      pivot_longer(-c(Time, Code, Name, run, cap, wgts, env, other)) %>%
      ggplot(aes(x = Time/365, y = value, color = factor(cap), linetype = factor(wgts))) +
      geom_line(linewidth = 0.5) +
      scale_color_manual(values = cap_col)+
      scale_linetype_manual(values = c("solid","dashed","dotted"))+
      scale_y_continuous(limits = c(0,NA)) +
      facet_grid2(env ~ name, scales = "free_y", independent = "y") +
      labs(title = current_name,
           x = "Year", 
           y = "",
           color = "Cap (mt)",
           linetype = "Weight scheme") +
      theme_bw()
    
    ggsave(paste0(plotdir, "/by_spp/", current_code, ".png"), p_tmp, 
           width = 10, height = 4.5,  
           units = "in", dpi = 300)
    
  }
  
  # for HCR species only: HCR plots
  # TODO: if we include estbo for all stocks in the PRM, we could plot these for all stocks to see how the OY makes the HCR look (even though there is no HCR with ramp)
  
  for(i in seq_along(oy_species)){
    
    current_code <- oy_species[i]
    current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
    
    # TODO: sort out faceting - you want the interesting feature to be in the same panel (e.g., facet by weight scheme or by climate scenario?)
    p2 <- catch_df %>%
      filter(Time >= burnin*365) %>%
      filter(Code == current_code, !is.na(catch_mt)) %>%
      ggplot(aes(x = biom_frac, y = f/fref, color = Time/365))+
      geom_point(aes(shape = factor(wgts)), size = 1)+
      scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
      scale_color_viridis_c(option = "viridis")+
      geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue", linewidth = 0.5)+
      geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", linewidth = 0.5)+
      geom_vline(xintercept = 0.2, linetype = "dotted", color = "orange", linewidth = 0.5)+
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "steelblue")+
      scale_y_continuous(limits = c(0,NA))+
      theme_bw()+
      labs(x = "B/B0", 
           y = "F/Ftarg", 
           color = "Year",
           shape = "Weight scheme",
           title = current_name)+
      facet_grid(factor(env)~cap)
    
    ggsave(paste0(plotdir, "/hcr/", current_code, "_hcr.png"), p2, 
           width = 10, height = 4.5, 
           units = "in", dpi = 300)
    
  }
  
  
  # catch of preferred species vs biomass of all species
  # preferred species should not be based on the we
  pref_catch <- catch_df %>%
    filter(!is.na(catch_mt)) %>%
    #left_join(pref, by = c("Code","wgts")) %>%
    left_join(pref, by = "Code") %>%
    filter(pref==1) %>%
    group_by(Time,run,cap,wgts,env)%>%
    summarise(catch_mt = sum(catch_mt))
  
  # then biom tot, join, plot
  p3 <- catch_df %>%
    filter(Time >= burnin*365) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time,run,cap,wgts,env)%>%
    summarise(biom_mt_tot = sum(biom_mt_tot)) %>%
    left_join(pref_catch) %>%
    ggplot(aes(x=biom_mt_tot, y = catch_mt, color = factor(cap), shape = factor(wgts)))+
    geom_point(size = 1)+
    scale_color_manual(values = cap_col)+
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    scale_x_continuous(limits = c(0,NA))+
    scale_y_continuous(limits = c(0,NA))+
    labs(x = "Biomass of all OY species (mt)",
         y = "Catch of pollock, cod, POP, sablefish (mt)",
         color = "Cap (mt)",
         shape = "Weight scheme",
         title = "Catch-biomass tradeoff")+
    theme_bw()+
    facet_wrap(~env)
  
  ggsave(paste0(plotdir, "/", "catch_vs_biomass.png"), p3, 
         width = 10, height = 4.5, 
         units = "in", dpi = 300)
  
  # Tradeoff plot between ATF biomass and POL biomass and their respective F
  pol <- catch_df %>%
    filter(Code == "POL") %>%
    select(Time, biom_frac, f, run, cap:env) %>%
    rename(biom_pol = biom_frac,
           f_pol = f)
  
  atf <- catch_df %>%
    filter(Code == "ATF") %>%
    select(Time, biom_frac, f, run, cap:env) %>%
    rename(biom_atf = biom_frac,
           f_atf = f)
  
  pol_vs_atf <- pol %>%
    left_join(atf) %>%
    mutate(f_ratio = f_atf / f_pol)
  
  p4 <- pol_vs_atf %>%
    filter(Time >= burnin*365) %>%
    filter(!is.na(f_ratio)) %>%
    filter(Time %in% (seq(0,yr_end,10)*365)) %>% # thin out
    ggplot(aes(x = biom_atf, y = biom_pol, color = f_ratio, shape = factor(wgts)))+
    geom_line(aes(group = interaction(Time, cap)), size = 0.5) +
    geom_point(aes(shape = factor(wgts)), size = 1.5)+
    geom_text(data = . %>% filter(wgts == "equal"), 
              aes(label = Time / 365), 
              nudge_x = -0.025, nudge_y = 0.025, 
              size = 2.5, color = "black", 
              check_overlap = F) +
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    viridis::scale_color_viridis(option = "cividis", begin = 0.05, end = 0.95)+
    theme_bw()+
    # scale_x_continuous(limits = c(0,NA))+
    # scale_y_continuous(limits = c(0,NA))+
    labs(x = "Arrowtooth B/B0",
         y = "Pollock B/B0",
         shape = "Weight scheme",
         color = "F(atf) / F(pol)",
         title = "Pollock-arrowtooth tradeoffs")+
    facet_grid(env~factor(cap))
  
  ggsave(paste0(plotdir, "/pol_vs_atf.png"), p4, 
         width = 10, height = 5, 
         units = "in", dpi = 300)
  
  # plot quantities relative to reference points and ecosystem indicators
  # some of these are repeated from the time series graphs so probably can do away with this plot, but it makes you compare species
  # summarize for early and late period to simplify the time dimension
  ecoind_df_tmp <- catch_df %>%
    filter(!is.na(catch_mt)) %>%
    rowwise() %>%
    mutate(period = ifelse(between(Time/365,burnin,burnin+10), "early",
                           ifelse(Time/365>(yr_end-10), "late", NA))) %>%
    ungroup() %>%
    filter(!is.na(period)) %>%
    group_by(run,cap,wgts,env,other,period,Code,Name,w,fref) %>%
    summarize(biom_mt_tot = mean(biom_mt_tot),
              catch_mt = mean(catch_mt),
              f = mean(f),
              oy_rescale = mean(oy_rescale)) %>%
    ungroup()
  
  # need data sets of total biomass and total catch
  ecoind_df_tot <- ecoind_df_tmp %>%
    group_by(run,cap,wgts,env,other,period) %>%
    summarise(tot_biom = sum(biom_mt_tot),
              tot_catch = sum(catch_mt)) %>%
    ungroup()
  
  ecoind_df <- ecoind_df_tmp %>%
    left_join(estbo_key) %>%
    left_join(ecoind_df_tot) %>%
    mutate(b_over_b0 = biom_mt_tot / estbo,
           f_over_ftarg = f/fref,
           biom_over_btot = biom_mt_tot / tot_biom,
           catch_over_ctot = catch_mt / tot_catch) %>%
    select(run:Name, oy_rescale, b_over_b0:catch_over_ctot)
  
  # pivot longer
  ecoind_df_long <- ecoind_df %>%
    pivot_longer(-c(run:Name)) %>%
    ungroup()
  
  # dot plot
  p5 <- ecoind_df_long %>%
    filter(env == "ssp585", Code %in% c("POL","COD","ATF","SBF","POP")) %>%
    ggplot(aes(x = name, y = value, color = factor(cap), shape = wgts))+
    geom_point(position = position_dodge(width = .75), size = 2)+
    scale_color_manual(values = cap_col)+
    theme_bw()+
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "", y = "", color = "Cap (mt)", shape = "Weight scheme", title = "ssp585")+
    facet_grid(Code~period)
  
  ggsave(paste0(plotdir, "/ecoind.png"), p5, 
         width = 12, height = 6.5, 
         units = "in", dpi = 300)
  
}

#' Set Boundary Box Values to NA
#'
#' @description
#' Helper function to set values in boundary boxes to NA for Atlantis spatial
#' analysis. Handles both 2D and 3D arrays/matrices.
#'
#' @param mat Numeric matrix or array. Spatial data from Atlantis model with
#'   boundary boxes that need to be excluded from analysis
#'
#' @return Matrix or array of same dimensions with boundary box values set to NA
#'
#' @details
#' Uses the global variable boundary_boxes to identify which spatial boxes
#' represent model boundaries. For 3D arrays, sets entire columns to NA.
#' For 2D matrices, sets entire rows to NA. Boundary boxes are typically
#' excluded from ecological analysis as they represent model edges.
#'
#' @examples
#' \dontrun{
#' clean_matrix <- setNA(spatial_data_matrix)
#' }
setNA <- function(mat) {
  mat2 <- mat
  if(length(dim(mat2))==3) mat2[,(boundary_boxes+1),]<-NA
  if(length(dim(mat2))==2) mat2[(boundary_boxes+1),] <- NA
  mat2
}

#' Calculate Shannon-Wiener Diversity Index
#'
#' @description
#' Computes the Shannon-Wiener diversity index for age structure of OY species
#' from Atlantis NetCDF output files. Can optionally weight by maturity.
#'
#' @param this_run Integer. The run number/ID for the Atlantis model simulation
#' @param do_mature Logical. If TRUE, weights abundance by maturity-at-age
#'   before calculating diversity index
#'
#' @return Data frame containing:
#'   \item{year}{Year of simulation}
#'   \item{Name}{Species name}
#'   \item{H}{Shannon-Wiener diversity index}
#'   \item{run}{Run number}
#'
#' @details
#' Reads numbers-at-age from NetCDF output files, sums across spatial boxes
#' (excluding boundaries), and calculates Shannon diversity index:
#' H = -sum(p_i * log(p_i)) where p_i is proportion in age class i.
#' 
#' When do_mature=TRUE, numbers are weighted by spawning potential (fspb)
#' before calculating proportions, giving a diversity index for spawning
#' population age structure.
#'
#' @examples
#' \dontrun{
#' diversity_data <- get_H(2097, do_mature = FALSE)
#' mature_diversity <- get_H(2097, do_mature = TRUE)
#' }
get_H <- function(this_run, do_mature){
  
  print(paste(this_run, "NAA"))
  
  # File paths
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  ncfile <- paste0(wd, "/outputGOA0", this_run, "_test.nc")
  
  # Open file
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # Pre-allocate list for results (much faster than rbind)
  h_list <- vector("list", length(oy_names))
  
  # do one fg at a time, then bring them back together
  for (i in seq_along(oy_names)){
    
    fg <- oy_names[i]
    sp <- grps %>% filter(Name == fg) %>% pull(Code) # need this for the fspb frame
    
    # Get numbers by box
    all_var_names <- names(this_ncdata$var)
    abun_var_names <- all_var_names[grepl("_Nums", all_var_names) & grepl(fg, all_var_names)]
    
    # Read all abundance data at once
    abun_data <- map(abun_var_names, ~ {
      data <- ncdf4::ncvar_get(this_ncdata, varid = .x)
      # Avoid setNA() - use direct indexing and na.rm instead
      # This eliminates the expensive lapply(setNA) bottleneck
      apply(data, MARGIN = 3, FUN = sum, na.rm = TRUE)
    })
    
    # Create abundance dataframe
    abun_df <- data.frame(
      t = tyrs,
      do.call(cbind, abun_data)
    )
    names(abun_df)[-1] <- abun_var_names
    
    # Process abundance data
    abun_processed <- abun_df %>%
      pivot_longer(cols = -t, names_to = 'age_group', values_to = 'abun') %>%
      mutate(
        age = parse_number(age_group) - 1,  # Atlantis indexing from 0
        year = ceiling(t),
        Name = fg
      ) %>%
      select(year, Name, age, abun)
    
    # Apply maturity adjustment if needed
    if(do_mature && exists("fspb_df")) {
      abun_processed <- abun_processed %>%
        left_join(fspb_df %>% filter(Code == sp), by = "age") %>%
        mutate(abun = abun * fspb)
    }
    
    # Calculate Shannon index
    h_result <- abun_processed %>%
      group_by(year, Name) %>%
      summarise(
        H = {
          p <- abun / sum(abun)
          p <- p[p > 0]  # Remove zeros to avoid log(0)
          -sum(p * log(p))
        },
        .groups = 'drop'
      )
    
    h_list[[i]] <- h_result
  }
  
  # Close the netCDF file
  nc_close(this_ncdata)
  
  # Combine results
  h_frame <- bind_rows(h_list) %>%
    mutate(run = this_run)
  
  return(h_frame)
  
}


#' Calculate Pollock Consumption Proportions by Predators
#'
#' @description
#' Extracts diet composition data to calculate the proportion of pollock
#' in the diets of identified pollock predators over time.
#'
#' @param this_run Integer. The run number/ID for the Atlantis model simulation
#'
#' @return Data frame containing:
#'   \item{Time}{Year of simulation}
#'   \item{Predator}{Predator species code}
#'   \item{POL}{Proportion of pollock in predator's diet}
#'   \item{run}{Run number}
#'
#' @details
#' Reads diet check output files and filters for predators identified in
#' POL_predators (species with >5% pollock in diet at the end of the burn-in).
#' Averages diet proportions across time steps within each year and across
#' age classes for each predator species. Used to analyze ecosystem effects
#' of pollock fishing on predator-prey dynamics.
#'
#' @examples
#' \dontrun{
#' pollock_diet_data <- get_polprop(2097)
#' }
get_polprop <- function(this_run){
  
  print(this_run)
  
  # File paths
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  dietfile <-  paste0("outputGOA0", this_run, "_testDietCheck.txt")
  diet <- fread(paste(wd, dietfile, sep = "/"), 
                select = c("Time", "Predator", "Cohort", "POL"))
  
  diet1 <- diet %>%
    mutate(Time = ceiling(Time/365)) %>%
    left_join(POL_predators) %>%
    filter(!is.na(isPred)) %>%
    group_by(Time, Predator) %>% # take the mean across the multiple time steps within each year AND across age classes for a predator
    summarise(POL = mean(POL)) 
  
  diet1 <- diet1 %>% mutate(run = this_run)
  return(diet1)
}
