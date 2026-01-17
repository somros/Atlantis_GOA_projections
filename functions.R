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
pull_fishery_info <- function(this_run, ecocap = T, ssb = T){
  
  print(this_run)
  
  # File paths
  # Old paths for prototype analyses on local
  # wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  # biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  # catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  # harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  # ecocap_file <- paste0(wd, "/outputGOA0", this_run, "_test_EcosystemCapResult.txt")
  
  # paths for new runs
  wd <- "../v2"
  biom_file <- paste0("outputGOA_", this_run, "AgeBiomIndx.txt")
  catch_file <- paste0("outputGOA_", this_run, "Catch.txt")
  harvest_prm <- run_combs %>% filter(run_id == as.numeric(this_run)) %>% pull(harvest_file)
  ecocap_file <- paste0(wd, "/outputGOA_", this_run, "_EcosystemCapResult.txt")
  
  
  # Read files once
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  harvest <- readLines(paste(oy_dir, harvest_prm, sep = "/"))
  
  # Force convert all columns except Time to numeric
  biom <- biom %>%
    mutate(across(-Time, ~as.numeric(as.character(.))))
  catch <- catch %>%
    mutate(across(-Time, ~as.numeric(as.character(.))))
  
  # Drop incomplete rows at the end of biom file
  # Check which rows have all numeric values (excluding Time column)
  complete_rows_biom <- apply(biom[, -1], 1, function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
  # Find the last complete row
  last_complete_biom <- max(which(complete_rows_biom))
  # Keep only up to the last complete row
  biom <- biom[1:last_complete_biom, ]
  
  # Drop incomplete rows at the end of catch file
  complete_rows_catch <- apply(catch[, -1], 1, function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
  last_complete_catch <- max(which(complete_rows_catch))
  catch <- catch[1:last_complete_catch, ]
  
  if(ecocap){
    ecocap_report <- read.delim(ecocap_file, sep = " ")
    # Remove columns that are all NA (ghost columns)
    ecocap_report <- ecocap_report[, colSums(is.na(ecocap_report)) < nrow(ecocap_report)]
    # Drop incomplete rows at the end of ecocap file
    complete_rows_ecocap <- apply(ecocap_report, 1, function(x) !any(is.na(x)))
    last_complete_ecocap <- max(which(complete_rows_ecocap))
    ecocap_report <- ecocap_report[1:last_complete_ecocap, ]
    
    ecocap_report <- ecocap_report %>%
      mutate(across(-c(Time, SpeciesName, FisheryName), ~as.numeric(as.character(.))))
  }
  
  # subset time early on if necessary
  biom <- biom %>% filter(Time/365 <= yr_end)
  catch <- catch %>% filter(Time/365 <= yr_end)
  
  if(ecocap){
    ecocap_report <- ecocap_report %>% filter(Time/365 <= yr_end)
  }
  
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
           #Time > 0, 
           Time < max(Time)) %>%
    select(-Code.Age)  # Remove original column
  
  # Pre-process harvest parameters outside loop
  harvest_params <- list()
  for(sp in oy_species) {
    startage_line <- harvest[grep(paste0(sp, "_mFC_startage"), harvest) + 1]
    startage <- as.numeric(strsplit(startage_line, split = " ")[[1]])[1]
    
    # get weight
    idx <- grep(sp, grps %>% pull(Code))
    w_line <- harvest[grep("SystCapSPpref", harvest) + 1]
    w_id <- as.numeric((strsplit(w_line, split = " "))[[1]])[idx]
    
    if(sp %in% hcr_spp) {
      sp_idx <- grep(sp, codes)
      estbo_vec <- as.numeric(strsplit(harvest[grep("estBo\t", harvest) + 1], split = " ")[[1]])
      estbo <- estbo_vec[sp_idx]
      fref_vec <- as.numeric(strsplit(harvest[grep("Fref\t", harvest) + 1], split = " ")[[1]])
      #fref <- -log(1 - fref_vec[sp_idx]) # TODO: pull this one from the reference points
    } else { # TODO: you can't take fref from the mfc file because those values are tuned so that the realized f is right
      estbo <- NA
      mfc_line <- harvest[grep(paste0("mFC_", sp, " "), harvest) + 1]
      mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
      #fref <- -(365) * log(1 - mfc) # TODO: pull this one from the F recent data frame? That will show us how good we're doing
    }
    
    fref <- fref_frame %>% filter(Code == sp) %>% pull(fref) # now pull fref from external data frame
    
    harvest_params[[sp]] <- list(startage = startage, w_id = w_id, estbo = estbo, fref = fref)
  }
  
  # Process all species at once using vectorized operations
  # Create biomass summaries for all species
  # biom_selex_all <- biom_filtered %>%
  #   left_join(
  #     data.frame(Code = names(harvest_params),
  #                startage = sapply(harvest_params, `[[`, "startage"),
  #                w = sapply(harvest_params, `[[`, "w")),
  #     by = "Code"
  #   ) %>%
  #   filter(Age >= startage) %>%
  #   group_by(Time, Code, w) %>%
  #   summarise(biom_mt_selex = sum(mt), .groups = 'drop')
  # 
  # if(ssb){
  #   
  #   biom_tot_all <- biom_filtered %>%
  #     left_join(fspb_df, by = c("Code", "Age"="age")) %>%
  #     group_by(Time, Code) %>%
  #     summarise(biom_mt_tot = sum(mt * fspb), .groups = 'drop')
  #   
  # } else {
  #   
  #   biom_tot_all <- biom_filtered %>%
  #     group_by(Time, Code) %>%
  #     summarise(biom_mt_tot = sum(mt), .groups = 'drop')
  #   
  # }
  
  
  # Process all species at once using vectorized operations
  # Create biomass summaries for all species
  # Albi: changed code so that it uses max biomass the year prior to account for recruitment / aging dynamics
  # NB: right now we have one time step per year so this should be identical to how it was before
  biom_selex_max <- biom_filtered %>%
    left_join(
      data.frame(Code = names(harvest_params),
                 startage = sapply(harvest_params, `[[`, "startage"),
                 w_id = sapply(harvest_params, `[[`, "w_id")),
      by = "Code"
    ) %>%
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
      summarise(biom_mt_tot = sum(mt), .groups = 'drop') %>%
      mutate(Year = Time / 365) %>%
      mutate(Year = floor(Year)) %>%
      group_by(Year, Code) %>%
      slice_max(biom_mt_tot) %>% # get the highest biomass record for the year
      ungroup() %>%
      mutate(Time = Year * 365) %>% # restore the Time column
      select(-Year)
    
  }
  
  # Process catch data for all species
  catch_long <- catch %>%
    select(Time, all_of(oy_species)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt")
  
  # join operation
  # to compute F from exploitation rates we need to stagger catch and biomass
  # catch at t356 / biomass at t0
  # this is not so important past the first year
  result_list <- list()
  for(sp in oy_species) {
    sp_params <- harvest_params[[sp]]
    
    sp_data <- biom_selex_max %>%
      filter(Code == sp) %>%
      left_join(filter(catch_long, Code == sp), by = c("Time", "Code")) %>%
      left_join(filter(biom_tot_all, Code == sp), by = c("Time", "Code")) %>%
      mutate(
        mu = catch_mt / lag(biom_mt_selex),
        f = -log(1 - mu),
        biom_frac = biom_mt_tot / sp_params$estbo,
        fref = sp_params$fref,
        run = this_run,
        w_id = sp_params$w_id
      )
    
    result_list[[sp]] <- sp_data
  }
  
  res_df <- do.call(rbind, result_list)
  
  # Process ecocap report
  if(ecocap){
    
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
    
  }
  
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
    # filter(Time > 0) %>%
    filter(Time >= burnin*365+1) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time,run,cap,wgts,env) %>%
    summarise(catch_tot = sum(catch_mt),
              biom_tot = sum(biom_mt_selex)) %>%
    ungroup() %>%
    pivot_longer(-c(Time:env)) %>%
    ggplot(aes(x = (Time/365)+1990, y = value, color = factor(cap), linetype = factor(wgts)))+
    # annotate("rect", xmin = 0+1990, xmax = burnin+1990, ymin = -Inf, ymax = Inf, 
    #          fill = "grey", alpha = 0.3) +
    geom_line(linewidth = 0.5)+
    scale_color_manual(values = cap_col)+
    #scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    scale_linetype_manual(values = c("solid","dashed","dotted"))+
    geom_hline(data = ~ filter(., name == "catch_tot"),
               aes(yintercept = as.numeric(as.character(cap))), 
               linetype = "dashed", 
               linewidth = 0.25)+    #geom_vline(xintercept = burnin, linetype = "dashed", color = "red")+
    theme_bw()+
    # scale_x_continuous(breaks = c(seq(1990,2100,10)))+
    scale_x_continuous(breaks = c(seq(2020,2100,10)))+
    scale_y_continuous(limits = c(0,NA))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "", 
         y = "Metric tons",
         color = "Cap (mt)",
         linetype = "Weight scheme")+
    facet_grid(name~env, scales = "free_y",
               labeller = labeller(name = as_labeller(
                 c("biom_tot" = "SSB", 
                   "catch_tot" = "Catch")
               )))
  
  ggsave(paste0(plotdir, "/", "oy_tot.png"), p1, 
         width = 10, height = 4.5,  
         units = "in", dpi = 300)
  
  # by species, biom fraction (need B0), f fraction, catch, biomass (selected), exploitation rate, ...?
  # maybe it would be best to produce this plot one species at a time
  # Loop through each code and create a plot
  
  for(i in seq_along(oy_species)) {
    current_code <- oy_species[i]
    current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
    
    p_df <- catch_df %>%
      filter(Time > 0) %>%
      #filter(Time >= burnin*365) %>%
      filter(Code == current_code) %>%
      filter(!is.na(catch_mt)) %>%
      mutate(f_frac = f / fref) %>%
      dplyr::select(-c(fref, f, biom_mt_tot, biom_mt_selex, mu, w_id)) %>%
      pivot_longer(-c(Time, Code, Name, run, cap, wgts, env))
    
    # drop the f_frac panel, redundant with the OY rescale panel mostly and we have the HCR information in the dedicated plots
    p_df <- p_df %>% filter(name != "f_frac")
    
    p_tmp <- ggplot(data = p_df %>% filter(Time > burnin*365+1), 
                    aes(x = (Time/365)+1990, y = value, color = factor(cap), linetype = factor(wgts))) +
      # annotate("rect", xmin = 0+1990, xmax = burnin+1990, ymin = -Inf, ymax = Inf, 
      #          fill = "grey", alpha = 0.3) +
      geom_line(linewidth = 0.5) +
      scale_color_manual(values = cap_col)+
      scale_linetype_manual(values = c("solid","dashed","dotted"))+
      # scale_x_continuous(breaks = c(seq(1990,2100,10))) +
      scale_x_continuous(breaks = c(seq(2020,2100,10))) +
      scale_y_continuous(limits = c(0,NA)) +
      geom_hline(data = data.frame(name = "biom_frac", 
                                   env = unique(catch_df$env)),
                 aes(yintercept = 1), linetype = "dotted") +
      {if(oy_species[i] %in% hcr_spp) 
        geom_hline(data = data.frame(name = "biom_frac", 
                                     env = unique(catch_df$env)),
                   aes(yintercept = 0.4), linetype = "dotted", color = "blue")} + 
      {if(oy_species[i] %in% hcr_spp) 
        geom_hline(data = data.frame(name = "biom_frac", 
                                     env = unique(catch_df$env)),
                   aes(yintercept = 0.2), linetype = "dotted", color = "red")} + 
      # geom_vline(xintercept = burnin, linetype = "dashed", color = "red")+
      facet_grid2(env ~ name, scales = "free_y", independent = "y") +
      labs(title = current_name,
           x = "", 
           y = "",
           color = "Cap (mt)",
           linetype = "Weight scheme") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
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
      filter(Time >= (burnin + 5)*365) %>%
      filter(Code == current_code, !is.na(catch_mt)) %>%
      ggplot(aes(x = biom_frac, y = f/fref, color = (Time/365)+1990))+
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
  # First, create pref_catch with decadal aggregation
  pref_catch <- catch_df %>%
    filter(!is.na(catch_mt)) %>%
    left_join(pref, by = "Code") %>%
    filter(pref==1) %>%
    group_by(Time, run, cap, wgts, env) %>%
    summarise(catch_mt = sum(catch_mt), .groups = 'drop') %>%
    mutate(decade = floor(Time/365/10)) %>%  # Create decade variable
    filter(decade %in% c(3, 10)) %>%  # Filter to decades 3, 10 (beginning and end of run)
    group_by(decade, run, cap, wgts, env) %>%
    summarise(
      catch_mt_mean = mean(catch_mt),
      catch_mt_sd = sd(catch_mt),
      .groups = 'drop'
    )
  
  # Then compute biomass totals and join
  p3 <- catch_df %>%
    filter(Time >= burnin*365) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time, run, cap, wgts, env) %>%
    summarise(biom_mt_tot = sum(biom_mt_tot), .groups = 'drop') %>%
    mutate(decade = floor(Time/365/10)) %>%  # Create decade variable
    filter(decade %in% c(3,10)) %>%  # Filter to decades 3, 10 (beginning and end of run)
    group_by(decade, run, cap, wgts, env) %>%
    summarise(
      biom_mt_tot_mean = mean(biom_mt_tot),
      biom_mt_tot_sd = sd(biom_mt_tot),
      .groups = 'drop'
    ) %>%
    left_join(pref_catch, by = c("decade", "run", "cap", "wgts", "env")) %>%
    ggplot(aes(x = biom_mt_tot_mean, y = catch_mt_mean, 
               color = factor(cap), shape = factor(wgts))) +
    geom_errorbar(aes(ymin = catch_mt_mean - catch_mt_sd, 
                      ymax = catch_mt_mean + catch_mt_sd),
                  width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = biom_mt_tot_mean - biom_mt_tot_sd,
                       xmax = biom_mt_tot_mean + biom_mt_tot_sd),
                   height = 0, alpha = 0.3) +
    geom_point(size = 2, position = position_dodge(width = 1)) +
    scale_color_manual(values = cap_col) +
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts)))) +
    # scale_x_continuous(limits = c(0, NA)) +
    # scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Biomass of all OY species (mt)",
         y = "Catch of pollock, cod, POP, sablefish (mt)",
         color = "Cap (mt)",
         shape = "Weight scheme",
         title = "Catch-biomass tradeoff") +
    theme_bw() +
    facet_grid(env ~ decade, 
               labeller = labeller(decade = as_labeller(
                 c("3" = "2020-2029", 
                   "10" = "2090-2099")
               )))
  
  ggsave(paste0(plotdir, "/", "catch_vs_biomass.png"), p3, 
         width = 8, height = 4.5, 
         units = "in", dpi = 300)
  
  # Tradeoff plot between ATF biomass and POL biomass and their respective F
  # TODO: turn this to mean per decade
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
  
  # p4 <- pol_vs_atf %>%
  #   filter(Time >= burnin*365) %>%
  #   filter(!is.na(f_ratio)) %>%
  #   filter(Time %in% (seq(0,yr_end,10)*365)) %>% # thin out
  #   ggplot(aes(x = biom_atf, y = biom_pol, color = f_ratio, shape = factor(wgts)))+
  #   geom_line(aes(group = interaction(Time, cap)), size = 0.5) +
  #   geom_point(aes(shape = factor(wgts)), size = 1.5)+
  #   geom_text(data = . %>% filter(wgts == "equal"), 
  #             aes(label = Time / 365), 
  #             nudge_x = -0.025, nudge_y = 0.025, 
  #             size = 2.5, color = "black", 
  #             check_overlap = F) +
  #   scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
  #   viridis::scale_color_viridis(option = "cividis", begin = 0.05, end = 0.95)+
  #   theme_bw()+
  #   # scale_x_continuous(limits = c(0,NA))+
  #   # scale_y_continuous(limits = c(0,NA))+
  #   labs(x = "Arrowtooth B/B0",
  #        y = "Pollock B/B0",
  #        shape = "Weight scheme",
  #        color = "F(atf) / F(pol)",
  #        title = "Pollock-arrowtooth tradeoffs")+
  #   facet_grid(env~factor(cap))
  
  # alternative version
  p4 <- pol_vs_atf %>%
    filter(Time >= burnin*365, !is.na(f_ratio)) %>%
    mutate(decade = floor(Time/365/10)) %>%
    # filter(decade %in% c(3, 10)) %>%
    filter(decade == 10, wgts != "binary") %>%
    group_by(decade, run, cap, wgts, env) %>%
    summarise(across(c(biom_pol, f_pol, biom_atf, f_atf, f_ratio),
                     list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"),
              .groups = 'drop') %>%
    filter(cap == "2e+05") %>% # because it is the only cap where catch streams from different wgts schemes still diverge at the end (other caps have become unconstraining)
    ggplot(aes(x = biom_atf_mean, y = biom_pol_mean, color = f_atf_mean, shape = factor(wgts)))+
    geom_errorbar(aes(ymin = biom_pol_mean - biom_pol_sd, 
                      ymax = biom_pol_mean + biom_pol_sd), 
                  width = 0, alpha = 0.9) +
    geom_errorbarh(aes(xmin = biom_atf_mean - biom_atf_sd, 
                       xmax = biom_atf_mean + biom_atf_sd), 
                   height = 0, alpha = 0.9) +
    geom_line(aes(group = interaction(decade, cap))) +
    geom_point(aes(shape = factor(wgts)), size = 1.5)+
    scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
    viridis::scale_color_viridis(option = "cividis", begin = 0.2, end = 0.8)+
    theme_bw()+
    labs(x = "Arrowtooth B/B0",
         y = "Pollock B/B0",
         shape = "Weight scheme",
         color = "F(atf)",
         title = "Pollock-arrowtooth tradeoffs")+
    facet_grid(~ env)#, 
               # labeller = labeller(decade = as_labeller(
               #   c("3" = "2020-2029", 
               #     "10" = "2090-2099")
               # )))

  
  ggsave(paste0(plotdir, "/pol_vs_atf.png"), p4, 
         width = 10, height = 3, 
         units = "in", dpi = 300)
  
  # plot quantities relative to reference points and ecosystem indicators
  # some of these are repeated from the time series graphs so probably can do away with this plot, but it makes you compare species
  # summarize for early and late period to simplify the time dimension
  # ecoind_df_tmp <- catch_df %>%
  #   filter(!is.na(catch_mt)) %>%
  #   rowwise() %>%
  #   mutate(period = ifelse(between(Time/365,burnin,burnin+10), "early",
  #                          ifelse(Time/365>(yr_end-10), "late", NA))) %>%
  #   ungroup() %>%
  #   filter(!is.na(period)) %>%
  #   group_by(run,cap,wgts,env,other,period,Code,Name,w_id,fref) %>%
  #   summarize(biom_mt_tot = mean(biom_mt_tot),
  #             catch_mt = mean(catch_mt),
  #             f = mean(f),
  #             oy_rescale = mean(oy_rescale)) %>%
  #   ungroup()
  # 
  # # need data sets of total biomass and total catch
  # ecoind_df_tot <- ecoind_df_tmp %>%
  #   group_by(run,cap,wgts,env,other,period) %>%
  #   summarise(tot_biom = sum(biom_mt_tot),
  #             tot_catch = sum(catch_mt)) %>%
  #   ungroup()
  # 
  # ecoind_df <- ecoind_df_tmp %>%
  #   left_join(estbo_key) %>%
  #   left_join(ecoind_df_tot) %>%
  #   mutate(b_over_b0 = biom_mt_tot / estbo,
  #          f_over_ftarg = f/fref,
  #          biom_over_btot = biom_mt_tot / tot_biom,
  #          catch_over_ctot = catch_mt / tot_catch) %>%
  #   select(run:Name, oy_rescale, b_over_b0:catch_over_ctot)
  # 
  # # pivot longer
  # ecoind_df_long <- ecoind_df %>%
  #   pivot_longer(-c(run:Name)) %>%
  #   ungroup()
  # 
  # # dot plot
  # p5 <- ecoind_df_long %>%
  #   filter(env == "ssp585", Code %in% c("POL","COD","ATF","SBF","POP")) %>%
  #   ggplot(aes(x = name, y = value, color = factor(cap), shape = wgts))+
  #   geom_point(position = position_dodge(width = .75), size = 2)+
  #   scale_color_manual(values = cap_col)+
  #   theme_bw()+
  #   geom_hline(yintercept = 1, linetype = "dashed")+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  #   labs(x = "", y = "", color = "Cap (mt)", shape = "Weight scheme", title = "ssp585")+
  #   facet_grid(Code~period)
  # 
  # ggsave(paste0(plotdir, "/ecoind.png"), p5, 
  #        width = 12, height = 6.5, 
  #        units = "in", dpi = 300)
  
}

#' Plot Impact of Ecosystem Caps on Total Catch
#'
#' @description
#' Creates bar plots showing the difference in total catch between capped
#' scenarios and the no-cap baseline (800,000 mt cap), using equal weight scheme.
#'
#' @param catch_df Data frame. Output from pull_fishery_info() containing
#'   fishery data across multiple runs and scenarios
#'
#' @return NULL (function creates and saves plot to plotdir)
#'
#' @details
#' Calculates total catch across all OY species for each time step, then
#' computes the difference between each capped scenario and the no-cap
#' baseline. Only shows results for equal weight scheme to avoid overcrowding.
#' Positive values indicate higher catch than baseline, negative values
#' indicate reduced catch due to cap constraints.
#' 
#' Plot shows:
#' \itemize{
#'   \item X-axis: Year
#'   \item Y-axis: Difference in catch (mt) from no-cap baseline
#'   \item Fill color: Cap level
#'   \item Facet by climate scenario
#'   \item Only post-burn-in period shown
#'   \item Reference line at y=0
#' }
#'
#' @examples
#' \dontrun{
#' plot_cap_impact(catch_df)
#' }
plot_catch_delta <- function(catch_df){
  
  # Calculate total catch and biomass by time and scenario
  totals <- catch_df %>%
    filter(Time > 0) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time, run, cap, wgts, env) %>%
    summarise(
      catch_tot = sum(catch_mt),
      biom_tot = sum(biom_mt_tot),
      .groups = 'drop'
    )
  
  # Get no-cap baseline (800,000 mt cap)
  baseline <- totals %>%
    filter(cap == 8e+05) %>%
    rename(
      catch_baseline = catch_tot,
      biom_baseline = biom_tot
    ) %>%
    select(Time, wgts, env, catch_baseline, biom_baseline)
  
  # Calculate difference from baseline
  delta <- totals %>%
    filter(cap != 8e+05, wgts == "equal") %>%
    left_join(baseline, by = c("Time", "wgts", "env")) %>%
    mutate(
      catch_delta = catch_tot - catch_baseline,
      biom_delta = biom_tot - biom_baseline
    ) %>%
    select(Time, run, cap, wgts, env, catch_delta, biom_delta) %>%
    pivot_longer(
      cols = c(catch_delta, biom_delta),
      names_to = "variable",
      values_to = "delta"
    ) %>%
    mutate(variable = recode(variable,
                             "catch_delta" = "Catch",
                             "biom_delta" = "Biomass"))
  
  # Make palette for capped scenarios only (excluding 8e+05)
  cap_col <- pnw_palette(name="Sunset2", 
                         n=length(unique(delta$cap)), 
                         type="discrete")
  
  # Create plot
  p <- delta %>%
    filter(Time > burnin * 365) %>%
    ggplot(aes(x = (Time/365) + 1990, y = delta, fill = factor(cap))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = cap_col) +
    scale_x_continuous(breaks = seq(1990, 2100, 10)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Year",
      y = "Difference from no cap (mt)",
      fill = "Cap (mt)",
      title = "Impact of ecosystem caps on total OY catch and biomass (relative to no cap)"
    ) +
    facet_grid(variable ~ env)
  
  ggsave(paste0(plotdir, "/catch_biom_cap_impact.png"), p, 
         width = 12, height = 6,  
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
  
  #TODO: turn this to decadal avg
  
  print(paste(this_run, "NAA"))
  
  # File paths
  #wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  wd <- "../v2"
  ncfile <- paste0(wd, "/outputGOA_", this_run, ".nc")
  
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
    abun_data <- map(abun_var_names, function(var_name) {
      raw_data <- ncdf4::ncvar_get(this_ncdata, varid = var_name)
      raw_data[, (boundary_boxes + 1), ] <- NA  # Same as setNA logic
      apply(raw_data, MARGIN = 3, FUN = sum, na.rm = TRUE)
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
  
  #TODO: avg by decade instead of full time series
  
  print(this_run)
  
  # File paths
  #wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  wd <- "../v2"
  dietfile <-  paste0(wd, "/outputGOA_", this_run, "DietCheck.txt")
  diet <- fread((dietfile), 
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

#' Calculate diet composition of OY groundfish for top predators
#'
#' Extracts and processes diet composition data from Atlantis model output to 
#' determine what proportion of top predator diets (Steller sea lions, pinnipeds, 
#' dolphins, and select demersal fish) consists of overfished/overfishing (OY) 
#' groundfish species. Results are averaged over the last 10 years of the simulation.
#'
#' @param this_run Character string specifying the Atlantis run identifier used to 
#'   locate the appropriate DietCheck.txt output file
#'
#' @return A tibble with the following columns:
#'   \itemize{
#'     \item \code{Predator}: Predator functional group code
#'     \item \code{Cohort}: Age cohort (0-indexed)
#'     \item \code{is_oy}: Binary indicator (1 = OY species)
#'     \item \code{comp_agg}: Aggregated diet composition (proportion)
#'     \item \code{run}: Run identifier
#'     \item \code{Age class}: Age class (1-indexed, Cohort + 1)
#'     \item \code{LongName}: Full name of predator group
#'     \item Additional columns from \code{key_config} joined by run
#'   }
#'
#' @details 
#' The function:
#' \itemize{
#'   \item Reads DietCheck.txt output from the specified Atlantis run
#'   \item Filters for the last 10 years of simulation (Time > 100)
#'   \item Focuses on five top predator groups: SSL, PIN, DOL, BDF, BSF
#'   \item Averages diet composition across years and age cohorts
#'   \item Identifies prey items in the OY groundfish category (defined in \code{oy_species})
#'   \item Returns only rows where OY species are present in the diet
#' }
#' 
#' @note Requires the following objects to be available in the environment:
#'   \code{oy_species}, \code{grps}, and \code{key_config}
#'
#' @examples
#' \dontrun{
#' diet_comp <- get_dietcomp_preds("1043")
#' }
#'
get_dietcomp_preds <- function(this_run, which_decade){
  
  print(this_run)
  
  # File paths
  #wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  wd <- "../v2"
  dietfile <-  paste0(wd, "/outputGOA_", this_run, "DietCheck.txt")
  diet <- fread((dietfile))
  
  diet1 <- diet %>%
    mutate(decade = floor(Time/365/10)) %>%
    filter(decade == which_decade) %>% # keep the last 10 years only
    select(-c(Stock, Updated)) %>%
    pivot_longer(-c(Time, Predator, Cohort,decade), names_to = "Prey", values_to = "comp") %>%
    filter(Predator %in% c("SSL", "PIN", "DOL", "BDF", "BSF")) %>%
    group_by(Predator, Cohort, Prey) %>% # average over last 10 years and across age classes for a predator
    summarise(comp = mean(comp), .groups = "drop") 
  
  # check that props add up to 1
  # diet1 %>% group_by(Predator, Cohort) %>% summarise(check = sum(comp)) %>% pull(check) %>% summary() # OK
  
  diet2 <- diet1 %>%
    mutate(is_oy = ifelse(Prey %in% oy_species, 1, 0)) %>%
    group_by(Predator, Cohort, is_oy) %>%
    summarise(comp_agg = sum(comp))
  
  diet2 <- diet2 %>% mutate(run = this_run) %>%
    filter(is_oy > 0) %>%
    mutate(`Age class` = Cohort + 1) %>%
    left_join(grps %>% select(Code, LongName), by = c("Predator"="Code")) %>%
    left_join(key_config, by = "run")
  
  return(diet2)
}


#' Calculate Ecosystem Indicators from Biomass Output
#'
#' @description
#' Calculates ecosystem indicators based on biomass ratios between functional
#' groups from Atlantis model output. Summarizes by decade for all decades
#' after the burn-in period.
#'
#' @param this_run Character. The run number/ID for the Atlantis model simulation
#'
#' @return Data frame containing:
#'   \item{decade}{Decade number}
#'   \item{indicator}{Indicator name}
#'   \item{value_mean}{Mean indicator value}
#'   \item{value_sd}{Standard deviation of indicator value}
#'   \item{run}{Run number}
#'   \item{cap}{Cap level}
#'   \item{wgts}{Weight scheme}
#'   \item{env}{Climate scenario}
#'
#' @details
#' Reads AgeBiomIndx.txt files and calculates four ecosystem indicators:
#' \itemize{
#'   \item gadids_flatfish: Gadids biomass / Flatfish biomass
#'   \item roundfish_flatfish: Roundfish biomass / Flatfish biomass
#'   \item shrimp_oy: Shrimp biomass / OY species biomass
#'   \item forage_oy: Forage fish biomass / OY species biomass
#' }
#'
#' @examples
#' \dontrun{
#' indicators <- calc_ecosystem_indicators("2097")
#' }
calc_ecosystem_indicators <- function(this_run){
  
  print(this_run)
  
  # Define functional groups
  gadids <- c("POL", "COD")
  flatfish <- c("ATF", "REX", "FFS", "FFD", "FHS")
  roundfish <- c("POL", "COD", "SBF")
  rockfish <- c("POP", "RFS", "RFD", "RFP", "THO")
  shrimp <- c("PAN", "PWN")
  forage <- c("HER", "SAN", "CAP", "FOS", "EUL")
  
  # File paths
  wd <- "../v2"
  biom_file <- paste0("outputGOA_", this_run, "AgeBiomIndx.txt")
  
  # Read biomass file
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  
  # Force convert all columns except Time to numeric
  biom <- biom %>%
    mutate(across(-Time, ~as.numeric(as.character(.))))
  
  # Drop incomplete rows at the end
  complete_rows_biom <- apply(biom[, -1], 1, function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
  last_complete_biom <- max(which(complete_rows_biom))
  biom <- biom[1:last_complete_biom, ]
  
  # Subset to time range
  biom <- biom %>% filter(Time/365 <= yr_end)
  
  # Convert to long format
  biom_long <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")
  
  # Split Code and Age
  code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
  biom_long$Code <- sapply(code_age_split, `[`, 1)
  biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
  
  # Calculate total biomass by species and time
  biom_by_species <- biom_long %>%
    group_by(Time, Code) %>%
    summarise(biom_mt = sum(mt, na.rm = TRUE), .groups = 'drop')
  
  # Calculate functional group biomasses
  biom_functional <- biom_by_species %>%
    mutate(
      gadids = ifelse(Code %in% gadids, biom_mt, 0),
      flatfish = ifelse(Code %in% flatfish, biom_mt, 0),
      roundfish = ifelse(Code %in% roundfish, biom_mt, 0),
      shrimp = ifelse(Code %in% shrimp, biom_mt, 0),
      forage = ifelse(Code %in% forage, biom_mt, 0),
      oy = ifelse(Code %in% oy_species, biom_mt, 0)
    ) %>%
    group_by(Time) %>%
    summarise(
      gadids_biom = sum(gadids),
      flatfish_biom = sum(flatfish),
      roundfish_biom = sum(roundfish),
      shrimp_biom = sum(shrimp),
      forage_biom = sum(forage),
      oy_biom = sum(oy),
      .groups = 'drop'
    )
  
  # Calculate indicators
  indicators <- biom_functional %>%
    mutate(
      gadids_to_flatfish = gadids_biom / flatfish_biom,
      roundfish_to_flatfish = roundfish_biom / flatfish_biom,
      shrimp_to_oy_spp = shrimp_biom / oy_biom,
      forage_to_oy_spp = forage_biom / oy_biom
    ) %>%
    select(Time, gadids_to_flatfish:forage_to_oy_spp)
  
  # Add decade, summarize
  indicators_decadal <- indicators %>%
    filter(Time >= burnin*365) %>%
    mutate(decade = floor(Time/365/10)) %>%
    pivot_longer(-c(Time, decade), names_to = "indicator", values_to = "value") %>%
    group_by(decade, indicator) %>%
    summarise(
      value_mean = mean(value, na.rm = TRUE),
      value_sd = sd(value, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(run = this_run)
  
  # Add metadata from key_config
  indicators_decadal <- indicators_decadal %>%
    left_join(key_config, by = "run")
  
  return(indicators_decadal)
}

#' Plot Ecosystem Indicators as Time Series
#'
#' @description
#' Creates time series plot of ecosystem indicators with decades on the x-axis,
#' faceted by indicator and climate scenario.
#'
#' @param indicators_df Data frame. Output from calc_ecosystem_indicators()
#'   containing ecosystem indicator values across runs and scenarios
#'
#' @return NULL (function creates and saves plot to plotdir)
#'
#' @details
#' Creates a plot with:
#' \itemize{
#'   \item X-axis: Decade
#'   \item Y-axis: Indicator value
#'   \item Points with error bars (no lines)
#'   \item Color by cap level
#'   \item Shape by weight scheme
#'   \item Facets: indicator (rows) × env (columns)
#'   \item Reference line at y=1
#' }
#'
#' @examples
#' \dontrun{
#' plot_ecosystem_indicators(indicators_df)
#' }
plot_ecosystem_indicators <- function(indicators_df){
  
  # turn caps into factors for better plotting
  indicators_df$cap <- as.character(indicators_df$cap)
  indicators_df$cap[is.na(indicators_df$cap)] <- "No cap"
  indicators_df$cap <- factor(indicators_df$cap, levels = c("8e+05","6e+05","4e+05","2e+05"))
  
  # order weigths
  indicators_df$wgts <- factor(indicators_df$wgts, levels = c("equal","binary","attainment-based"))
  
  # Make palette
  cap_col <- pnw_palette(name="Sunset2", n=length(unique(indicators_df$cap)), type="discrete")
  
  # filter out what we don't need for this paper
  indicators_df <- indicators_df %>%
    filter(wgts != "binary", decade %in% c(6,10), indicator != "roundfish_to_flatfish", env %in% c("NoClimate", "ssp585"))
  
  # add decade label
  indicators_df <- indicators_df %>%
    mutate(decade_lab = ifelse(decade == 6, "2050-2059", "2090-2099"))
  
  # Time series plot with points and error bars
  p <- indicators_df %>%
    ggplot(aes(x = factor(decade_lab), y = value_mean, 
               color = factor(cap), 
               shape = factor(wgts))) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = value_mean - value_sd, 
                      ymax = value_mean + value_sd),
                  width = 0.2, alpha = 0.9,
                  position = position_dodge(width = 0.3)) +
    scale_color_manual(values = cap_col) +
    scale_shape_manual(values = c(1:length(unique(indicators_df$wgts)))) +
    scale_y_continuous(limits = c(0, NA)) +
    # scale_x_continuous(breaks = unique(indicators_df$decade)) +
    theme_bw() +
    # theme(
    #   strip.text = element_text(size = 10),
    #   legend.position = "bottom"
    # ) +
    labs(
      x = "",
      y = "Indicator value (mean ± SD)",
      color = "Cap (mt)",
      shape = "Weight scheme"
    ) +
    facet_grid(env ~ indicator, 
               labeller = labeller(indicator = as_labeller(
                 c("forage_to_oy_spp" = "Forage fish / groundfish", 
                   "gadids_to_flatfish" = "Gadids / flatfish",
                   "shrimp_to_oy_spp" = "Shrimp / groundfish")
               )))
  
  ggsave(paste0(plotdir, "/ecosystem_indicators.png"), p,
         width = 9, height = 4,
         units = "in", dpi = 300)
  
}

#' #' Plot Ecosystem Indicators as Radar Charts
#' #'
#' #' @description
#' #' Creates radar (spider) plots of ecosystem indicators using raw ratio values,
#' #' with color representing weight schemes. Allows filtering by cap, env, and decade.
#' #'
#' #' @param indicators_df Data frame. Output from calc_ecosystem_indicators()
#' #'   containing ecosystem indicator values across runs and scenarios
#' #' @param cap_filter Numeric or NULL. Cap value to filter to (default NULL for all)
#' #' @param env_filter Character or NULL. Environment to filter to (default NULL for all)
#' #' @param decade_filter Numeric or NULL. Decade to filter to (default NULL for all)
#' #'
#' #' @return NULL (function creates and saves plots to plotdir)
#' #'
#' #' @details
#' #' Creates radar plots showing four ecosystem indicators using raw ratio values.
#' #' Uses ggradar package with color representing weight schemes.
#' #' Filter to specific combinations to avoid overcrowding.
#' #'
#' #' @examples
#' #' \dontrun{
#' #' # Single specific scenario
#' #' plot_ecosystem_indicators_radar(indicators_all, cap_filter = 4e+05, 
#' #'                                  env_filter = "ssp585", decade_filter = 9)
#' #' 
#' #' # All decades for one cap and env
#' #' plot_ecosystem_indicators_radar(indicators_all, cap_filter = 4e+05, 
#' #'                                  env_filter = "ssp585", decade_filter = NULL)
#' #' }
#' plot_ecosystem_indicators_radar <- function(indicators_df, 
#'                                             cap_filter = NULL, 
#'                                             env_filter = NULL,
#'                                             decade_filter = NULL){
#'   
#'   library(ggradar)
#'   
#'   # Color palette for weight schemes
#'   wgts_col <- c("equal" = "#2E5266", "binary" = "#D65108", "attainment-based" = "#6B9080")
#'   
#'   # Filter data based on arguments
#'   radar_data <- indicators_df
#'   
#'   if(!is.null(cap_filter)) {
#'     radar_data <- radar_data %>% filter(cap == cap_filter)
#'   }
#'   
#'   if(!is.null(env_filter)) {
#'     radar_data <- radar_data %>% filter(env == env_filter)
#'   }
#'   
#'   if(!is.null(decade_filter)) {
#'     radar_data <- radar_data %>% filter(decade == decade_filter)
#'   }
#'   
#'   # Check if data remains
#'   if(nrow(radar_data) == 0) {
#'     stop("No data remaining after filtering. Check filter values.")
#'   }
#'   
#'   # Determine min and max across all indicators for consistent scaling
#'   global_min <- min(radar_data$value_mean, na.rm = TRUE)
#'   global_max <- max(radar_data$value_mean, na.rm = TRUE)
#'   
#'   # Prepare data for ggradar (needs wide format)
#'   radar_prep <- radar_data %>%
#'     select(decade, env, cap, wgts, indicator, value_mean) %>%
#'     mutate(group = as.character(wgts)) %>%
#'     select(decade, env, cap, group, indicator, value_mean) %>%
#'     pivot_wider(names_from = indicator, values_from = value_mean)
#'   
#'   # Create separate plots for remaining combinations
#'   plot_list <- list()
#'   idx <- 1
#'   
#'   remaining_caps <- unique(radar_prep$cap)
#'   remaining_envs <- unique(radar_prep$env)
#'   remaining_decades <- unique(radar_prep$decade)
#'   
#'   for(c in remaining_caps) {
#'     for(e in remaining_envs) {
#'       for(d in remaining_decades) {
#'         plot_data <- radar_prep %>%
#'           filter(env == e, decade == d, cap == c) %>%
#'           select(-env, -decade, -cap)
#'         
#'         if(nrow(plot_data) > 0) {
#'           # Match colors to weight schemes
#'           wgts_levels <- unique(plot_data$group)
#'           plot_colors <- wgts_col[match(wgts_levels, names(wgts_col))]
#'           
#'           # Create title based on what's being shown
#'           title_parts <- c()
#'           if(length(remaining_caps) > 1) title_parts <- c(title_parts, paste("Cap:", c))
#'           if(length(remaining_envs) > 1) title_parts <- c(title_parts, e)
#'           if(length(remaining_decades) > 1) title_parts <- c(title_parts, paste("Decade", d))
#'           
#'           plot_title <- if(length(title_parts) > 0) {
#'             paste(title_parts, collapse = " | ")
#'           } else {
#'             paste("Cap:", c, "|", e, "| Decade", d)
#'           }
#'           
#'           p_tmp <- ggradar(plot_data,
#'                            grid.min = floor(global_min * 10) / 10,
#'                            grid.mid = (floor(global_min * 10) / 10 + ceiling(global_max * 10) / 10) / 2,
#'                            grid.max = ceiling(global_max * 10) / 10,
#'                            values.radar = c(
#'                              as.character(floor(global_min * 10) / 10),
#'                              as.character(round((floor(global_min * 10) / 10 + ceiling(global_max * 10) / 10) / 2, 1)),
#'                              as.character(ceiling(global_max * 10) / 10)
#'                            ),
#'                            group.colours = plot_colors,
#'                            group.point.size = 3,
#'                            group.line.width = 1.2,
#'                            gridline.mid.colour = "grey70",
#'                            gridline.min.colour = "grey90",
#'                            gridline.max.colour = "grey90",
#'                            background.circle.colour = "white",
#'                            axis.label.size = 4,
#'                            grid.label.size = 3.5,
#'                            legend.position = "bottom",
#'                            legend.text.size = 10,
#'                            legend.title = "Weight scheme") +
#'             labs(title = plot_title) +
#'             theme(plot.title = element_text(size = 11, hjust = 0.5))
#'           
#'           plot_list[[idx]] <- p_tmp
#'           idx <- idx + 1
#'         }
#'       }
#'     }
#'   }
#'   
#'   # Combine plots if multiple
#'   if(length(plot_list) > 1) {
#'     library(patchwork)
#'     
#'     # Determine layout
#'     n_plots <- length(plot_list)
#'     if(length(remaining_decades) > 1) {
#'       ncol_val <- length(remaining_decades)
#'     } else if(length(remaining_envs) > 1) {
#'       ncol_val <- length(remaining_envs)
#'     } else {
#'       ncol_val <- min(3, n_plots)
#'     }
#'     
#'     combined_plot <- wrap_plots(plot_list, ncol = ncol_val)
#'     
#'     # Create filename suffix based on filters
#'     suffix_parts <- c()
#'     if(!is.null(cap_filter)) suffix_parts <- c(suffix_parts, paste0("cap", cap_filter))
#'     if(!is.null(env_filter)) suffix_parts <- c(suffix_parts, env_filter)
#'     if(!is.null(decade_filter)) suffix_parts <- c(suffix_parts, paste0("dec", decade_filter))
#'     
#'     filename_suffix <- if(length(suffix_parts) > 0) {
#'       paste0("_", paste(suffix_parts, collapse = "_"))
#'     } else {
#'       ""
#'     }
#'     
#'     filename_suffix <- gsub("\\+", "", filename_suffix)
#'     
#'     ggsave(paste0(plotdir, "/ecosystem_indicators_radar", filename_suffix, ".png"), 
#'            combined_plot,
#'            width = 5 * ncol_val, 
#'            height = 5 * ceiling(n_plots / ncol_val), 
#'            units = "in", dpi = 300)
#'   } else {
#'     # Single plot
#'     p <- plot_list[[1]]
#'     
#'     # Create filename suffix
#'     suffix_parts <- c()
#'     if(!is.null(cap_filter)) suffix_parts <- c(suffix_parts, paste0("cap", cap_filter))
#'     if(!is.null(env_filter)) suffix_parts <- c(suffix_parts, env_filter)
#'     if(!is.null(decade_filter)) suffix_parts <- c(suffix_parts, paste0("dec", decade_filter))
#'     
#'     filename_suffix <- if(length(suffix_parts) > 0) {
#'       paste0("_", paste(suffix_parts, collapse = "_"))
#'     } else {
#'       ""
#'     }
#'     
#'     filename_suffix <- gsub("\\+", "", filename_suffix)
#'     
#'     ggsave(paste0(plotdir, "/ecosystem_indicators_radar", filename_suffix, ".png"), 
#'            p, width = 6, height = 6, units = "in", dpi = 300)
#'   }
#'   
#' }
#' 
#' #' Calculate Total Revenue from Catch
#'
#' @description
#' Calculates total revenue from catch of OY species using species-specific
#' prices and catch data from Atlantis model output.
#'
#' @param this_run Character. The run number/ID for the Atlantis model simulation
#' @param price_dat Data frame. Contains Code and mean_price columns ($/lb)
#'
#' @return Data frame containing:
#'   \item{Time}{Time step from model output}
#'   \item{revenue}{Total revenue in dollars}
#'   \item{run}{Run number}
#'   \item{cap}{Cap level}
#'   \item{wgts}{Weight scheme}
#'   \item{env}{Climate scenario}
#'
#' @details
#' Reads Catch.txt files, filters to OY species, converts catch from metric
#' tons to pounds (1 mt = 2204.62 lbs), multiplies by species-specific prices,
#' and sums across all OY species to get total revenue.
#'
#' @examples
#' \dontrun{
#' price_dat <- read.csv("data/price.csv")
#' revenue <- calc_revenue("2097", price_dat)
#' }
calc_revenue <- function(this_run, price_dat){
  
  print(this_run)
  
  # Conversion factor
  mt_to_lbs <- 2204.62
  
  # File paths
  wd <- "../v2"
  catch_file <- paste0("outputGOA_", this_run, "Catch.txt")
  
  # Read catch file
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  
  # Force convert all columns except Time to numeric
  catch <- catch %>%
    mutate(across(-Time, ~as.numeric(as.character(.))))
  
  # Drop incomplete rows at the end
  complete_rows_catch <- apply(catch[, -1], 1, function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
  last_complete_catch <- max(which(complete_rows_catch))
  catch <- catch[1:last_complete_catch, ]
  
  # Subset to time range
  catch <- catch %>% filter(Time/365 <= yr_end)
  
  # Convert to long format and filter to OY species
  catch_long <- catch %>%
    select(Time, all_of(oy_species)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt")
  
  # Join with prices, convert to lbs, calculate revenue
  revenue_data <- catch_long %>%
    left_join(price_dat, by = "Code") %>%
    mutate(
      catch_lbs = catch_mt * mt_to_lbs,
      revenue = catch_lbs * mean_price
    ) %>%
    group_by(Time) %>%
    summarise(revenue = sum(revenue, na.rm = TRUE), .groups = 'drop') %>%
    mutate(run = this_run)
  
  # Add metadata from key_config
  revenue_data <- revenue_data %>%
    left_join(key_config, by = "run")
  
  return(revenue_data)
}

#' Plot Total Revenue from OY Species
#'
#' @description
#' Creates time series plots of total revenue from OY species catch,
#' distinguished by cap levels and weight schemes.
#'
#' @param revenue_df Data frame. Output from calc_revenue()
#'   containing revenue data across runs and scenarios
#'
#' @return NULL (function creates and saves plot to plotdir)
#'
#' @details
#' Creates a line plot with:
#' \itemize{
#'   \item X-axis: Year
#'   \item Y-axis: Total revenue ($)
#'   \item Color by cap level
#'   \item Line type by weight scheme
#'   \item Facet by climate scenario
#'   \item Shaded burn-in period
#' }
#'
#' @examples
#' \dontrun{
#' plot_revenue(revenue_df)
#' }
plot_revenue <- function(revenue_df){
  
  # Make palette
  cap_col <- pnw_palette(name="Sunset2", n=length(unique(revenue_df$cap)), type="discrete")
  
  # Revenue time series plot
  p <- revenue_df %>%
    filter(Time > 0) %>%
    ggplot(aes(x = Time/365, y = revenue, 
               color = factor(cap), 
               linetype = factor(wgts))) +
    annotate("rect", xmin = 0, xmax = burnin, ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.3) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = cap_col) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(limits = c(0, NA), labels = scales::dollar_format()) +
    theme_bw() +
    labs(
      x = "Year", 
      y = "Total revenue",
      color = "Cap (mt)",
      linetype = "Weight scheme",
      title = "Total revenue from OY species"
    ) +
    facet_grid(~env)
  
  ggsave(paste0(plotdir, "/oy_revenue.png"), p, 
         width = 10, height = 3,  
         units = "in", dpi = 300)
}

#' Plot Catch, Biomass, and Revenue Deltas Across Weight Schemes (Decadal Averages)
#'
#' @description
#' Creates stacked bar plots showing decadal average differences in catch, biomass, 
#' and revenue between capped scenarios and the no-cap baseline (800,000 mt cap) 
#' across different weight schemes. Also produces a summary table of cumulative differences.
#'
#' @param catch_df Data frame. Output from pull_fishery_info() containing
#'   fishery data across multiple runs and scenarios
#' @param revenue_df Data frame. Output from calc_revenue() containing
#'   revenue data across multiple runs and scenarios
#'
#' @return NULL (function creates and saves stacked plots and summary table to plotdir)
#'
#' @details
#' Filters to NoClimate scenario only to focus on cap effects without climate
#' confounding. The 800,000 mt cap baseline (which only has "equal" weights)
#' is expanded to all weight schemes since caps don't apply at that level.
#' 
#' Aggregates to decadal means with standard deviations shown as error bars.
#' 
#' Creates three stacked panels:
#' \itemize{
#'   \item Catch differences (mt)
#'   \item Biomass differences (mt)
#'   \item Revenue differences ($)
#' }
#' Each panel is faceted by weight scheme.
#'
#' @examples
#' \dontrun{
#' price_dat <- read.csv("data/price.csv")
#' all_runs <- key_config %>% pull(run)
#' revenue_all <- map_df(all_runs, ~calc_revenue(.x, price_dat))
#' plot_delta(catch_df, revenue_all)
#' }
#' Plot Catch, Biomass, and Revenue Deltas Across Weight Schemes (Decadal Averages)
#'
#' @description
#' Creates stacked bar plots showing decadal average differences in catch, biomass, 
#' and revenue between capped scenarios and the no-cap baseline (800,000 mt cap) 
#' across different weight schemes. Also produces a summary table of cumulative differences.
#'
#' @param catch_df Data frame. Output from pull_fishery_info() containing
#'   fishery data across multiple runs and scenarios
#' @param revenue_df Data frame. Output from calc_revenue() containing
#'   revenue data across multiple runs and scenarios
#'
#' @return NULL (function creates and saves stacked plots and summary table to plotdir)
#'
#' @details
#' Filters to NoClimate scenario only to focus on cap effects without climate
#' confounding. The 800,000 mt cap baseline (which only has "equal" weights)
#' is expanded to all weight schemes since caps don't apply at that level.
#' 
#' Aggregates to decadal means with standard deviations shown as error bars.
#' 
#' Creates three stacked panels:
#' \itemize{
#'   \item Catch differences (mt)
#'   \item Biomass differences (mt)
#'   \item Revenue differences ($)
#' }
#' Each panel is faceted by weight scheme.
#'
#' @examples
#' \dontrun{
#' price_dat <- read.csv("data/price.csv")
#' all_runs <- key_config %>% pull(run)
#' revenue_all <- map_df(all_runs, ~calc_revenue(.x, price_dat))
#' plot_delta(catch_df, revenue_all)
#' }
plot_delta <- function(catch_df, revenue_df){
  
  # Filter to NoClimate only
  catch_df <- catch_df %>% filter(env == "NoClimate")
  revenue_df <- revenue_df %>% filter(env == "NoClimate")
  
  # Ensure cap is numeric in both dataframes for joining
  catch_df <- catch_df %>% mutate(cap = as.numeric(as.character(cap)))
  revenue_df <- revenue_df %>% mutate(cap = as.numeric(cap))
  
  # Calculate total catch and biomass by time and scenario
  totals <- catch_df %>%
    filter(Time > 0) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time, run, cap, wgts, env) %>%
    summarise(
      catch_tot = sum(catch_mt),
      biom_tot = sum(biom_mt_tot),
      .groups = 'drop'
    )
  
  # Process revenue data
  revenue_totals <- revenue_df %>%
    filter(Time > 0) %>%
    select(Time, run, cap, wgts, env, revenue)
  
  # Get no-cap baseline (800,000 mt cap) - only has "equal" wgts
  baseline_catch_biom <- totals %>%
    filter(cap == 8e+05) %>%
    rename(
      catch_baseline = catch_tot,
      biom_baseline = biom_tot
    ) %>%
    select(Time, env, catch_baseline, biom_baseline)
  
  baseline_revenue <- revenue_totals %>%
    filter(cap == 8e+05) %>%
    rename(revenue_baseline = revenue) %>%
    select(Time, env, revenue_baseline)
  
  # Expand baseline to all wgts levels (since it doesn't matter for 800K cap)
  all_wgts <- unique(totals$wgts)
  baseline_catch_biom <- baseline_catch_biom %>%
    crossing(wgts = all_wgts)
  
  baseline_revenue <- baseline_revenue %>%
    crossing(wgts = all_wgts)
  
  # Calculate differences from baseline
  delta_catch_biom <- totals %>%
    filter(cap != 8e+05) %>%
    left_join(baseline_catch_biom, by = c("Time", "wgts", "env")) %>%
    mutate(
      catch_delta = catch_tot - catch_baseline,
      biom_delta = biom_tot - biom_baseline
    )
  
  delta_revenue <- revenue_totals %>%
    filter(cap != 8e+05) %>%
    left_join(baseline_revenue, by = c("Time", "wgts", "env")) %>%
    mutate(revenue_delta = revenue - revenue_baseline)
  
  # Join catch/biom and revenue deltas
  delta_all <- delta_catch_biom %>%
    left_join(delta_revenue %>% select(Time, run, cap, wgts, env, revenue_delta),
              by = c("Time", "run", "cap", "wgts", "env"))
  
  # Filter to post-burnin
  delta_filtered <- delta_all %>%
    filter(Time > burnin * 365)
  
  # Filter out binary
  delta_filtered <- delta_filtered %>%
    filter(wgts != "binary")
  
  # Add decade variable
  delta_filtered <- delta_filtered %>%
    mutate(decade = floor((Time/365 - burnin) / 10))
  
  # Order wgts and cap for plotting
  delta_filtered$wgts <- factor(delta_filtered$wgts, 
                                levels = c("equal", "binary", "attainment-based"))
  delta_filtered$cap <- factor(delta_filtered$cap, levels = c("8e+05","6e+05","4e+05","2e+05"))
  
  # Make palette consistent with plot_fishery (4 colors, then subset to 3)
  all_cap_col <- pnw_palette(name="Sunset2", n=4, type="discrete")
  cap_col <- all_cap_col[2:4]
  
  # Prepare catch and biomass for combined plot
  delta_catch_biom_long <- delta_filtered %>%
    select(Time, decade, run, cap, wgts, env, catch_delta, biom_delta) %>%
    pivot_longer(
      cols = c(catch_delta, biom_delta),
      names_to = "variable",
      values_to = "delta"
    ) %>%
    mutate(variable = recode(variable,
                             "catch_delta" = "Catch",
                             "biom_delta" = "Biomass"))
  
  # Aggregate to decadal averages for catch and biomass
  catch_biom_decadal <- delta_catch_biom_long %>%
    group_by(decade, cap, wgts, env, variable) %>%
    summarise(
      delta_mean = mean(delta, na.rm = TRUE),
      delta_sd = sd(delta, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Aggregate to decadal averages for revenue
  revenue_decadal <- delta_filtered %>%
    group_by(decade, cap, wgts, env) %>%
    summarise(
      revenue_delta_mean = mean(revenue_delta, na.rm = TRUE),
      revenue_delta_sd = sd(revenue_delta, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Order factors
  catch_biom_decadal$wgts <- factor(catch_biom_decadal$wgts, 
                                    levels = c("equal", "binary", "attainment-based"))
  catch_biom_decadal$cap <- factor(catch_biom_decadal$cap, 
                                   levels = c("8e+05","6e+05","4e+05","2e+05"))
  
  revenue_decadal$wgts <- factor(revenue_decadal$wgts, 
                                 levels = c("equal", "binary", "attainment-based"))
  revenue_decadal$cap <- factor(revenue_decadal$cap, 
                                levels = c("8e+05","6e+05","4e+05","2e+05"))
  
  # Create combined catch/biomass plot
  p_catch_biom <- catch_biom_decadal %>%
    mutate(variable = factor(variable, levels = c("Catch", "Biomass"))) %>%
    ggplot(aes(x = decade, y = delta_mean, fill = factor(cap))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_bar(stat = "identity", position = position_dodge(),
             color = "darkgrey", linewidth = 0.3) +
    geom_errorbar(aes(ymin = delta_mean - delta_sd, 
                      ymax = delta_mean + delta_sd),
                  width = 0.2,
                  position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cap_col) +
    scale_x_continuous(breaks = unique(catch_biom_decadal$decade),
                       labels = seq(2020, 2090, 10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 10)) +
    labs(
      y = "Difference from no cap (mt)",
      fill = "Cap (mt)"#,
      #title = "Impact of ecosystem caps on OY catch, biomass, and revenue (NoClimate scenario)"
    ) +
    facet_grid(variable ~ wgts)
  
  # Create revenue plot
  p_revenue <- revenue_decadal %>%
    ggplot(aes(x = decade, y = revenue_delta_mean, fill = factor(cap))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_bar(stat = "identity", position = position_dodge(),
             color = "darkgrey", linewidth = 0.3) +
    geom_errorbar(aes(ymin = revenue_delta_mean - revenue_delta_sd, 
                      ymax = revenue_delta_mean + revenue_delta_sd),
                  width = 0.2,
                  position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cap_col) +
    scale_x_continuous(breaks = unique(revenue_decadal$decade),
                       labels = seq(2020, 2090, 10)) +
    scale_y_continuous(labels = scales::dollar_format()) +
    theme_bw() +
    theme(strip.text = element_text(size = 10)) +
    labs(
      x = "Decade",
      y = "Revenue difference",
      fill = "Cap (mt)"
    ) +
    facet_wrap(~wgts, nrow = 1)
  
  # Stack plots using patchwork
  library(patchwork)
  p_combined <- p_catch_biom / p_revenue +
    plot_layout(heights = c(2, 1), guides = "collect") +
    plot_annotation(tag_levels = 'A') &
    theme(legend.position = "bottom")
  
  ggsave(paste0(plotdir, "/catch_biom_revenue_delta.png"), p_combined,
         width = 10, height = 7, units = "in", dpi = 300)
  
  # Create summary table
  summary_table <- delta_filtered %>%
    group_by(cap, wgts, env) %>%
    summarise(
      total_catch_loss = sum(catch_delta),
      total_biom_gain = sum(biom_delta),
      total_revenue_loss = sum(revenue_delta),
      .groups = 'drop'
    ) %>%
    arrange(wgts, cap)
  
  write.csv(summary_table, 
            paste0(plotdir, "/summary_catch_biom_revenue_tradeoff.csv"), 
            row.names = FALSE)
  
  # Print to console
  cat("\n=== Summary of Catch-Biomass-Revenue Tradeoffs (NoClimate) ===\n\n")
  print(summary_table, n = Inf)
  cat("\n")
  
}

#' Plot Catch Loss vs Ecosystem Biomass Gains (Decadal Averages)
#'
#' @description
#' Creates stacked bar plots showing decadal average catch losses from OY species 
#' alongside biomass gains in non-target ecosystem components (forage fish, 
#' top predators, and Pacific halibut total biomass) under different cap scenarios.
#'
#' @param catch_df Data frame. Output from pull_fishery_info() containing
#'   fishery data across multiple runs and scenarios
#'
#' @return NULL (function creates and saves plots to plotdir)
#'
#' @details
#' Filters to NoClimate scenario only. Reads biomass data directly from
#' AgeBiomIndx.txt files for:
#' \itemize{
#'   \item Forage fish: CAP, EUL, SAN, FOS, HER (total biomass)
#'   \item Top predators: SSL, PIN, BDF, BSF (total biomass)
#'   \item Pacific halibut: HAL (total biomass)
#' }
#' 
#' Aggregates to decadal means with standard deviations shown as error bars.
#'
#' @examples
#' \dontrun{
#' plot_ecosystem_delta(catch_df)
#' }
plot_ecosystem_delta <- function(catch_df){
  
  # Filter to NoClimate only
  catch_df <- catch_df %>% filter(env == "NoClimate")
  
  # Ensure cap is numeric
  catch_df <- catch_df %>% mutate(cap = as.numeric(as.character(cap)))
  
  # Define functional groups
  forage <- c("CAP", "EUL", "SAN", "FOS", "HER")
  predators <- c("SSL", "PIN", "BDF", "BSF")
  halibut <- "HAL"
  
  # Get all runs
  all_runs <- unique(catch_df$run)
  
  # Function to read and process biomass for one run
  get_ecosystem_biomass <- function(this_run){
    
    print(paste("Processing ecosystem biomass for run", this_run))
    
    # File paths
    wd <- "../v2"
    biom_file <- paste0("outputGOA_", this_run, "AgeBiomIndx.txt")
    
    # Read biomass file
    biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
    
    # Force convert all columns except Time to numeric
    biom <- biom %>%
      mutate(across(-Time, ~as.numeric(as.character(.))))
    
    # Drop incomplete rows
    complete_rows_biom <- apply(biom[, -1], 1, function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
    last_complete_biom <- max(which(complete_rows_biom))
    biom <- biom[1:last_complete_biom, ]
    
    # Subset to time range
    biom <- biom %>% filter(Time/365 <= yr_end)
    
    # Convert to long format
    biom_long <- biom %>%
      pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")
    
    # Split Code and Age
    code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
    biom_long$Code <- sapply(code_age_split, `[`, 1)
    biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
    
    # Calculate forage biomass
    forage_biom <- biom_long %>%
      filter(Code %in% forage) %>%
      group_by(Time) %>%
      summarise(forage_biom = sum(mt, na.rm = TRUE), .groups = 'drop')
    
    # Calculate predator biomass
    predator_biom <- biom_long %>%
      filter(Code %in% predators) %>%
      group_by(Time) %>%
      summarise(predator_biom = sum(mt, na.rm = TRUE), .groups = 'drop')
    
    # Calculate halibut total biomass
    hal_biom <- biom_long %>%
      filter(Code %in% halibut) %>%
      group_by(Time) %>%
      summarise(hal_biom = sum(mt, na.rm = TRUE), .groups = 'drop')
    
    # Combine all metrics
    ecosystem_biom <- forage_biom %>%
      left_join(predator_biom, by = "Time") %>%
      left_join(hal_biom, by = "Time") %>%
      mutate(run = this_run)
    
    return(ecosystem_biom)
  }
  
  # Get ecosystem biomass for all runs
  ecosystem_all <- map_df(all_runs, get_ecosystem_biomass)
  
  # Add metadata from key_config
  ecosystem_all <- ecosystem_all %>%
    left_join(key_config, by = "run") %>%
    filter(env == "NoClimate")
  
  # Ensure cap is numeric
  ecosystem_all <- ecosystem_all %>% mutate(cap = as.numeric(cap))
  
  # Get no-cap baseline (800,000 mt cap)
  baseline_ecosystem <- ecosystem_all %>%
    filter(cap == 8e+05) %>%
    rename(
      forage_baseline = forage_biom,
      predator_baseline = predator_biom,
      hal_baseline = hal_biom
    ) %>%
    select(Time, env, forage_baseline, predator_baseline, hal_baseline)
  
  # Expand baseline to all wgts levels
  all_wgts <- unique(ecosystem_all$wgts)
  baseline_ecosystem <- baseline_ecosystem %>%
    crossing(wgts = all_wgts)
  
  # Calculate differences from baseline
  delta_ecosystem <- ecosystem_all %>%
    filter(cap != 8e+05) %>%
    left_join(baseline_ecosystem, by = c("Time", "wgts", "env")) %>%
    mutate(
      forage_delta = forage_biom - forage_baseline,
      predator_delta = predator_biom - predator_baseline,
      hal_delta = hal_biom - hal_baseline
    )
  
  # Calculate catch deltas
  totals <- catch_df %>%
    filter(Time > 0) %>%
    filter(!is.na(catch_mt)) %>%
    group_by(Time, run, cap, wgts, env) %>%
    summarise(catch_tot = sum(catch_mt), .groups = 'drop')
  
  baseline_catch <- totals %>%
    filter(cap == 8e+05) %>%
    rename(catch_baseline = catch_tot) %>%
    select(Time, env, catch_baseline) %>%
    crossing(wgts = all_wgts)
  
  delta_catch <- totals %>%
    filter(cap != 8e+05) %>%
    left_join(baseline_catch, by = c("Time", "wgts", "env")) %>%
    mutate(catch_delta = catch_tot - catch_baseline)
  
  # Join catch and ecosystem deltas
  delta_all <- delta_catch %>%
    left_join(delta_ecosystem %>% select(Time, run, cap, wgts, env, 
                                         forage_delta, predator_delta, hal_delta),
              by = c("Time", "run", "cap", "wgts", "env"))
  
  # Filter to post-burnin and exclude binary
  delta_filtered <- delta_all %>%
    filter(Time > burnin * 365) %>%
    filter(wgts != "binary")
  
  # Add decade variable
  delta_filtered <- delta_filtered %>%
    mutate(decade = floor((Time/365 - burnin) / 10))
  
  # Aggregate to decadal averages for catch
  catch_decadal <- delta_filtered %>%
    group_by(decade, cap, wgts, env) %>%
    summarise(
      catch_delta_mean = mean(catch_delta, na.rm = TRUE),
      catch_delta_sd = sd(catch_delta, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Prepare ecosystem metrics for plotting
  delta_ecosystem_long <- delta_filtered %>%
    select(Time, decade, run, cap, wgts, env, forage_delta, predator_delta, hal_delta) %>%
    pivot_longer(
      cols = c(forage_delta, predator_delta, hal_delta),
      names_to = "variable",
      values_to = "delta"
    ) %>%
    mutate(variable = recode(variable,
                             "forage_delta" = "Forage fish",
                             "predator_delta" = "Top predators",
                             "hal_delta" = "Halibut"))
  
  # Aggregate to decadal averages for ecosystem metrics
  ecosystem_decadal <- delta_ecosystem_long %>%
    group_by(decade, cap, wgts, env, variable) %>%
    summarise(
      delta_mean = mean(delta, na.rm = TRUE),
      delta_sd = sd(delta, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Order factors
  catch_decadal$wgts <- factor(catch_decadal$wgts, 
                               levels = c("equal", "binary", "attainment-based"))
  catch_decadal$cap <- factor(catch_decadal$cap, 
                              levels = c("8e+05","6e+05","4e+05","2e+05"))
  
  ecosystem_decadal$wgts <- factor(ecosystem_decadal$wgts, 
                                   levels = c("equal", "binary", "attainment-based"))
  ecosystem_decadal$cap <- factor(ecosystem_decadal$cap, 
                                  levels = c("8e+05","6e+05","4e+05","2e+05"))
  
  # Make palette
  all_cap_col <- pnw_palette(name="Sunset2", n=4, type="discrete")
  cap_col <- all_cap_col[2:4]
  
  # Create catch plot
  p_catch <- catch_decadal %>%
    filter(wgts == "equal") %>% # do not need wgts info in this plot
    ggplot(aes(x = decade, y = catch_delta_mean, fill = factor(cap))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_bar(stat = "identity", position = position_dodge(), 
             color = "darkgrey", linewidth = 0.3) +
    geom_errorbar(aes(ymin = catch_delta_mean - catch_delta_sd, 
                      ymax = catch_delta_mean + catch_delta_sd),
                  width = 0.2,
                  position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cap_col) +
    scale_x_continuous(breaks = unique(catch_decadal$decade),
                       labels = seq(2020, 2090, 10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 10)) +
    labs(
      y = "Catch difference (mt)",
      fill = "Cap (mt)"#,
      #title = "Impact of ecosystem caps on OY catch and ecosystem biomass (NoClimate scenario)"
    )
  
  # Create ecosystem metrics plot
  p_ecosystem <- ecosystem_decadal %>%
    filter(wgts == "equal") %>% # do not need wgts info in this plot
    ggplot(aes(x = decade, y = delta_mean, fill = factor(cap))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_bar(stat = "identity", position = position_dodge(),
             color = "darkgrey", linewidth = 0.3) +
    geom_errorbar(aes(ymin = delta_mean - delta_sd, 
                      ymax = delta_mean + delta_sd),
                  width = 0.2,
                  position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cap_col) +
    scale_x_continuous(breaks = unique(ecosystem_decadal$decade),
                       labels = seq(2020, 2090, 10)) +
    theme_bw() +
    theme(strip.text = element_text(size = 10)) +
    labs(
      x = "Decade",
      y = "Biomass difference (mt)",
      fill = "Cap (mt)"
    ) +
    facet_wrap(~variable, nrow = 3, scales = "free_y")
  
  # Stack plots using patchwork
  library(patchwork)
  p_combined <- p_catch / p_ecosystem +
    plot_layout(heights = c(1, 3), guides = "collect") +
    plot_annotation(tag_levels = 'A') &
    theme(legend.position = "bottom")
  
  ggsave(paste0(plotdir, "/catch_ecosystem_delta.png"), p_combined,
         width = 5, height = 8, units = "in", dpi = 300)
  
  # Create summary table
  summary_table <- delta_filtered %>%
    group_by(cap, wgts, env) %>%
    summarise(
      total_catch_loss = sum(catch_delta),
      total_forage_gain = sum(forage_delta),
      total_predator_gain = sum(predator_delta),
      total_hal_gain = sum(hal_delta),
      .groups = 'drop'
    ) %>%
    arrange(wgts, cap)
  
  write.csv(summary_table, 
            paste0(plotdir, "/summary_catch_ecosystem_tradeoff.csv"), 
            row.names = FALSE)
  
  # Print to console
  cat("\n=== Summary of Catch-Ecosystem Tradeoffs (NoClimate) ===\n\n")
  print(summary_table, n = Inf)
  cat("\n")
  
}

#' Calculate Total Biomass for Lower Trophic Level Groups
#'
#' @description
#' Calculates total biomass for specified lower trophic level groups
#' (zooplankton, phytoplankton) from Atlantis model output.
#'
#' @param this_run Character. The run number/ID for the Atlantis model simulation
#' @param groups Character vector. Codes of groups to extract 
#'   (default c("EUP", "ZL", "ZM", "ZS", "PL", "PS"))
#'
#' @return Data frame containing:
#'   \item{Time}{Time step from model output}
#'   \item{Code}{Species/group code}
#'   \item{biom_mt}{Total biomass in metric tons}
#'   \item{run}{Run number}
#'   \item{cap}{Cap level}
#'   \item{wgts}{Weight scheme}
#'   \item{env}{Climate scenario}
#'
#' @details
#' Reads AgeBiomIndx.txt files, sums biomass across all age classes
#' for specified groups. These are typically lower trophic level groups
#' that don't require age-specific analysis.
#'
#' @examples
#' \dontrun{
#' ltl_biomass <- calc_ltl_biomass("2097")
#' }
calc_ltl_biomass <- function(this_run, groups = c("EUP", "ZL", "ZM", "ZS", "PL", "PS")){
  
  print(this_run)
  
  # File paths
  wd <- "../v2"
  biom_file <- paste0("outputGOA_", this_run, "AgeBiomIndx.txt")
  
  # Read biomass file
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  
  # Force convert all columns except Time to numeric
  biom <- biom %>%
    mutate(across(-Time, ~as.numeric(as.character(.))))
  
  # Drop incomplete rows at the end
  complete_rows_biom <- apply(biom[, -1], 1, function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x))))))
  last_complete_biom <- max(which(complete_rows_biom))
  biom <- biom[1:last_complete_biom, ]
  
  # Subset to time range
  biom <- biom %>% filter(Time/365 <= yr_end)
  
  # Convert to long format
  biom_long <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")
  
  # Split Code and Age
  code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
  biom_long$Code <- sapply(code_age_split, `[`, 1)
  biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
  
  # Filter to specified groups and sum across ages
  biom_ltl <- biom_long %>%
    filter(Code %in% groups) %>%
    group_by(Time, Code) %>%
    summarise(biom_mt = sum(mt, na.rm = TRUE), .groups = 'drop') %>%
    mutate(run = this_run)
  
  # Add metadata from key_config
  biom_ltl <- biom_ltl %>%
    left_join(key_config, by = "run")
  
  return(biom_ltl)
}

#' Plot Total Biomass for Lower Trophic Level Groups
#'
#' @description
#' Creates time series plots of total biomass for lower trophic level groups,
#' faceted by group and climate scenario.
#'
#' @param ltl_df Data frame. Output from calc_ltl_biomass()
#'   containing biomass data across runs and scenarios
#'
#' @return NULL (function creates and saves plot to plotdir)
#'
#' @details
#' Creates a line plot with:
#' \itemize{
#'   \item X-axis: Year
#'   \item Y-axis: Total biomass (mt)
#'   \item Color by cap level
#'   \item Line type by weight scheme
#'   \item Facet by env (rows) × Code (columns)
#'   \item Free y-axis scales (requires ggh4x package)
#'   \item Shaded burn-in period
#' }
#'
#' @examples
#' \dontrun{
#' plot_ltl_biomass(ltl_df)
#' }
plot_ltl_biomass <- function(ltl_df){
  
  # Make palette
  cap_col <- pnw_palette(name="Sunset2", n=length(unique(ltl_df$cap)), type="discrete")
  
  # Get group names for better labels
  ltl_df <- ltl_df %>%
    left_join(grps %>% select(Code, Name), by = "Code")
  
  # Biomass time series plot
  p <- ltl_df %>%
    filter(Time > 0) %>%
    ggplot(aes(x = Time/365, y = biom_mt, 
               color = factor(cap), 
               linetype = factor(wgts))) +
    annotate("rect", xmin = 0, xmax = burnin, ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.3) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = cap_col) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 9)
    ) +
    labs(
      x = "Year", 
      y = "Biomass (mt)",
      color = "Cap (mt)",
      linetype = "Weight scheme",
      title = "Lower trophic level biomass"
    ) +
    ggh4x::facet_grid2(env ~ Name, scales = "free_y", independent = "y")
  
  ggsave(paste0(plotdir, "/ltl_biomass.png"), p, 
         width = 14, height = 6,  
         units = "in", dpi = 300)
}

#' Extract Weight-at-Age from NetCDF Output
#'
#' @description
#' Extracts abundance-weighted mean weight-at-age for specified functional groups
#' from Atlantis NetCDF output files. Calculates weight from reserve and
#' structural nitrogen pools weighted by spatial distribution of abundance.
#'
#' @param this_run Character. The run number/ID for the Atlantis model simulation
#' @param sp_names Character vector. Species names to extract (e.g., c("Pollock", "Cod")).
#'   Default is oy_names (all OY species names)
#' @param boundary_boxes Numeric vector. Box IDs that are boundaries (0-indexed).
#'   If provided, values in these boxes will be set to NA and excluded from calculations.
#'   Default is NULL (no boundary box filtering)
#'
#' @return Data frame containing:
#'   \item{year}{Year of simulation}
#'   \item{age_group}{Age group identifier (e.g., "Pollock1_ResN")}
#'   \item{age}{Numeric age}
#'   \item{weight}{Weight in kg}
#'   \item{Name}{Short species name}
#'   \item{LongName}{Full species name}
#'   \item{run}{Run number}
#'   \item{cap}{Cap level}
#'   \item{wgts}{Weight scheme}
#'   \item{env}{Climate scenario}
#'
#' @details
#' Reads reserve nitrogen (_ResN), structural nitrogen (_StructN), and
#' numbers (_Nums) from NetCDF files. Calculates spatial abundance-weighted
#' mean nitrogen content, then converts to weight in kg using the conversion:
#' weight (kg) = totN * 20 * 5.7 / 1,000,000
#' 
#' If boundary_boxes is provided, sets values in those boxes to NA before
#' calculations to exclude boundary effects from weight-at-age estimates.
#'
#' @examples
#' \dontrun{
#' # Without boundary box filtering
#' waa_data <- calc_weight_at_age("2097")
#' 
#' # With boundary box filtering
#' boundary_boxes <- c(0, 1, 2)  # Example boundary boxes (0-indexed)
#' waa_data <- calc_weight_at_age("2097", boundary_boxes = boundary_boxes)
#' }
calc_weight_at_age <- function(this_run, sp_names = oy_names, boundary_boxes = NULL){
  
  print(paste(this_run, "weight-at-age"))
  
  # File paths
  wd <- "../v2"
  ncfile <- paste0(wd, "/outputGOA_", this_run, ".nc")
  
  # Open file
  this_nc <- nc_open(ncfile)
  
  # Get time
  ts <- ncdf4::ncvar_get(this_nc, varid = "t") %>% as.numeric()
  tyrs <- ts/(60*60*24*365)
  
  # Get all variable names
  all_var_names <- names(this_nc$var)
  
  # Filter for reserve N, structural N, and numbers for specified species
  # Pattern is: SpeciesName + digit(s) + _ResN (or _StructN or _Nums)
  resN_vars <- all_var_names[sapply(sp_names, function(x) grepl(paste0("^", x, "[0-9]+_ResN$"), all_var_names)) %>% 
                               apply(1, any)]
  strucN_vars <- all_var_names[sapply(sp_names, function(x) grepl(paste0("^", x, "[0-9]+_StructN$"), all_var_names)) %>% 
                                 apply(1, any)]
  abun_vars <- all_var_names[sapply(sp_names, function(x) grepl(paste0("^", x, "[0-9]+_Nums$"), all_var_names)) %>% 
                               apply(1, any)]
  
  if(length(resN_vars) == 0) {
    nc_close(this_nc)
    message(paste("No ResN data found for run", this_run, "with specified species"))
    return(NULL)
  }
  
  # Actually pull the data from the .nc
  resN <- purrr::map(resN_vars, ncdf4::ncvar_get, nc = this_nc) 
  strucN <- purrr::map(strucN_vars, ncdf4::ncvar_get, nc = this_nc)
  nums <- purrr::map(abun_vars, ncdf4::ncvar_get, nc = this_nc) # numbers by age group, box, layer, time
  
  # Set boundary boxes to NA if provided
  # Arrays are 3D with dimensions [layer, box, time]
  if(!is.null(boundary_boxes)) {
    # Convert 0-indexed to 1-indexed for R
    boundary_idx <- boundary_boxes + 1
    
    # Set boundary boxes to NA in all arrays
    resN <- purrr::map(resN, function(x) {
      x[, boundary_idx, ] <- NA  # Box is the second dimension
      return(x)
    })
    
    strucN <- purrr::map(strucN, function(x) {
      x[, boundary_idx, ] <- NA  # Box is the second dimension
      return(x)
    })
    
    nums <- purrr::map(nums, function(x) {
      x[, boundary_idx, ] <- NA  # Box is the second dimension
      return(x)
    })
  }
  
  # Calculate total numbers and relative numbers
  totnums <- nums %>% purrr::map(apply, MARGIN = 3, FUN = sum, na.rm = TRUE) # total numbers by age group, time
  relnums <- purrr::map2(nums, totnums, sweep, MARGIN = 3, FUN = `/`) # divide nums by totnums along the time axis
  
  # Add the two matrices to get total nitrogen weight
  rnsn <- purrr::map2(resN, strucN, `+`)
  
  # Multiply and sum to get abundance-weighted mean weight at age
  rnsn_summ <- purrr::map2(rnsn, relnums, `*`) %>% 
    purrr::map(apply, MARGIN = 3, FUN = sum, na.rm = TRUE) %>% # mean total N by time (na.rm added)
    bind_cols() %>% # bind age groups elements together
    suppressMessages() %>% 
    set_names(resN_vars) %>% # TODO: change this, this is ResN + StrucN
    mutate(t = tyrs) %>%
    # pivot to long form
    pivot_longer(cols = -t, names_to = 'age_group', values_to = 'totN') %>%
    mutate(age = parse_number(age_group)) %>% 
    mutate(weight = totN * 20 * 5.7 / 1000000) %>%   # convert totN to weight/individual in kg
    dplyr::filter(t > 0) %>%
    mutate(year = t) %>%
    group_by(year, age_group, age) %>%
    summarise(weight = mean(weight), .groups = 'drop') %>%
    ungroup() %>%
    mutate(Name = gsub('[0-9]+_ResN', '', age_group)) %>%
    left_join(grps %>% dplyr::select(Name, LongName), by = 'Name') %>%
    mutate(run = this_run)
  
  # Close the netCDF file
  nc_close(this_nc)
  
  # Add metadata from key_config
  rnsn_summ <- rnsn_summ %>%
    left_join(key_config, by = "run")
  
  return(rnsn_summ)
}

#' Extract Numbers-at-Age from NetCDF Output
#'
#' @description
#' Extracts total numbers-at-age for specified functional groups
#' from Atlantis NetCDF output files by summing across spatial boxes.
#'
#' @param this_run Character. The run number/ID for the Atlantis model simulation
#' @param sp_names Character vector. Species names to extract (e.g., c("Pollock", "Cod")).
#'   Default is oy_names (all OY species names)
#' @param boundary_boxes Numeric vector. Box IDs that are boundaries (0-indexed).
#'   If provided, values in these boxes will be set to NA and excluded from calculations.
#'   Default is NULL (no boundary box filtering)
#'
#' @return Data frame containing:
#'   \item{year}{Year of simulation}
#'   \item{age_group}{Age group identifier (e.g., "Pollock1_Nums")}
#'   \item{age}{Numeric age}
#'   \item{abun}{Total abundance (numbers)}
#'   \item{Name}{Short species name}
#'   \item{LongName}{Full species name}
#'   \item{run}{Run number}
#'   \item{cap}{Cap level}
#'   \item{wgts}{Weight scheme}
#'   \item{env}{Climate scenario}
#'
#' @details
#' Reads numbers (_Nums) from NetCDF files and sums across spatial boxes
#' (excluding boundary boxes if specified) to get total abundance-at-age
#' over time.
#'
#' @examples
#' \dontrun{
#' # Without boundary box filtering
#' naa_data <- calc_numbers_at_age("2097")
#' 
#' # With boundary box filtering
#' boundary_boxes <- c(0, 1, 2)  # Example boundary boxes (0-indexed)
#' naa_data <- calc_numbers_at_age("2097", boundary_boxes = boundary_boxes)
#' }
calc_numbers_at_age <- function(this_run, sp_names = oy_names, boundary_boxes = NULL){
  
  print(paste(this_run, "numbers-at-age"))
  
  # File paths
  wd <- "../v2"
  ncfile <- paste0(wd, "/outputGOA_", this_run, ".nc")
  
  # Open file
  this_nc <- nc_open(ncfile)
  
  # Get time
  ts <- ncdf4::ncvar_get(this_nc, varid = "t") %>% as.numeric()
  tyrs <- ts/(60*60*24*365)
  
  # Get all variable names
  all_var_names <- names(this_nc$var)
  
  # Filter for numbers for specified species
  # Pattern is: SpeciesName + digit(s) + _Nums
  abun_vars <- all_var_names[sapply(sp_names, function(x) grepl(paste0("^", x, "[0-9]+_Nums$"), all_var_names)) %>% 
                               apply(1, any)]
  
  if(length(abun_vars) == 0) {
    nc_close(this_nc)
    message(paste("No Nums data found for run", this_run, "with specified species"))
    return(NULL)
  }
  
  # Actually pull the data from the .nc
  nums <- purrr::map(abun_vars, ncdf4::ncvar_get, nc = this_nc) # numbers by layer, box, time
  
  # Set boundary boxes to NA if provided
  # Arrays are 3D with dimensions [layer, box, time]
  if(!is.null(boundary_boxes)) {
    # Convert 0-indexed to 1-indexed for R
    boundary_idx <- boundary_boxes + 1
    
    # Set boundary boxes to NA in all arrays
    nums <- purrr::map(nums, function(x) {
      x[, boundary_idx, ] <- NA  # Box is the second dimension
      return(x)
    })
  }
  
  # Sum across space (layers and boxes) for each time point
  abun_summ <- nums %>% 
    purrr::map(apply, MARGIN = 3, FUN = sum, na.rm = TRUE) %>% 
    bind_cols() %>% 
    suppressMessages() %>% 
    set_names(abun_vars) %>% 
    mutate(t = tyrs) %>%
    # pivot to long form
    pivot_longer(cols = -t, names_to = 'age_group', values_to = 'abun') %>%
    mutate(age = parse_number(age_group)) %>%
    mutate(year = t) %>%
    group_by(year, age_group, age) %>%
    summarise(abun = mean(abun), .groups = 'drop') %>%
    ungroup() %>%
    mutate(Name = gsub('[0-9]+_Nums', '', age_group)) %>%
    left_join(grps %>% dplyr::select(Name, LongName), by = 'Name') %>%
    mutate(run = this_run)
  
  # Close the netCDF file
  nc_close(this_nc)
  
  # Add metadata from key_config
  abun_summ <- abun_summ %>%
    left_join(key_config, by = "run")
  
  return(abun_summ)
}

#' Create a heatmap of weight-at-age and/or numbers-at-age changes relative to baseline
#'
#' This function creates heatmap visualizations showing percentage changes in
#' weight-at-age and/or numbers-at-age data compared to a baseline scenario 
#' (NoClimate with equal weights). The comparison can be made across environmental 
#' scenarios, capacity levels, or weight schemes. When plotting both metrics, they
#' are stacked vertically using faceting.
#'
#' @param waa_df Data frame containing weight-at-age data with columns: year, Name,
#'   LongName, age, env, cap, wgts, and weight. Required if plot_type is "WAA" or "both"
#' @param naa_df Data frame containing numbers-at-age data with columns: year, Name,
#'   LongName, age, env, cap, wgts, and abun. Required if plot_type is "NAA" or "both"
#' @param plot_type Character. One of "WAA", "NAA", or "both". Determines which
#'   metric(s) to plot. Default is "both"
#' @param by_env Logical. If TRUE, plot comparison across environmental scenarios.
#'   Default is FALSE
#' @param by_cap Logical. If TRUE, plot comparison across capacity levels.
#'   Default is FALSE
#' @param by_wgts Logical. If TRUE, plot comparison across weight schemes.
#'   Default is FALSE
#' @param n_years Integer. Number of most recent years to include in the analysis.
#'   Default is 5
#' @param tag Character string. Optional tag to append to the output filename.
#'   Default is NULL
#'
#' @return Invisibly returns a combined data frame containing ratio data with a 
#'   'metric' column indicating WAA or NAA. If no data remains after filtering, 
#'   returns NULL with a message
#'
#' @details
#' Exactly one of by_env, by_cap, or by_wgts must be TRUE. The function:
#' \itemize{
#'   \item Filters data to integer year values (January 1st)
#'   \item Calculates mean values over the last n_years
#'   \item Compares to baseline: NoClimate scenario, cap=8e+05 (or 4e+05 for by_wgts),
#'         and equal weights
#'   \item Creates diverging color heatmap(s) (blue-white-red) showing percent change
#'   \item Saves plot to plotdir as PNG (10x height, 300 dpi)
#' }
#'
#' @note Requires ggplot2, dplyr, and colorspace packages. Assumes plotdir 
#'   is defined in the global environment
#'
#' @examples
#' \dontrun{
#' # Plot both WAA and NAA by environmental scenario
#' plot_age_heatmap(waa_data, naa_data, plot_type = "both", 
#'                  by_env = TRUE, n_years = 5, tag = "climate")
#'
#' # Plot only WAA by capacity level
#' plot_age_heatmap(waa_df = waa_data, plot_type = "WAA",
#'                  by_cap = TRUE, n_years = 10, tag = "capacity")
#' }
#'
#' @export
plot_age_heatmap <- function(waa_df = NULL, 
                             naa_df = NULL,
                             plot_type = c("both", "WAA", "NAA"),
                             by_env = FALSE,
                             by_cap = FALSE,
                             by_wgts = FALSE,
                             which_decade = NULL,
                             tag = NULL) {
  
  # Validate inputs
  plot_type <- match.arg(plot_type)
  
  if(is.null(which_decade)){
    stop("Please enter a decade as input")
  }
  
  if(sum(by_env, by_cap, by_wgts) != 1) {
    stop("You need to plot either by env, cap, or wgts")
  }
  
  if(plot_type %in% c("WAA", "both") && is.null(waa_df)) {
    stop("waa_df is required when plot_type is 'WAA' or 'both'")
  }
  
  if(plot_type %in% c("NAA", "both") && is.null(naa_df)) {
    stop("naa_df is required when plot_type is 'NAA' or 'both'")
  }
  
  # Determine label for faceting and filename
  if(by_env) {
    lab <- "env"
    facet_var <- "env"
  } else if(by_cap) {
    lab <- "cap"
    facet_var <- "cap"
  } else {
    lab <- "wgts"
    facet_var <- "wgts"
  }
  
  # Helper function to process data
  process_age_data <- function(df, value_col, metric_name) {
    
    # Filter to desired decade
    df_filtered <- df %>%
      mutate(decade = floor(year/10)) %>%
      filter(decade == which_decade)
    
    # calculate means
    df_ny <- df_filtered %>%
      group_by(Name, LongName, age, env, cap, wgts) %>%
      summarise(value_mean = mean(.data[[value_col]], na.rm = TRUE), .groups = 'drop')
    
    # Create baseline
    df_baseline <- df_ny %>% 
      filter(env == "NoClimate", wgts == "equal") %>%
      rename(value_baseline = value_mean)
    
    if(by_wgts) {
      df_baseline <- df_baseline %>% filter(cap == 4e+05)
    } else {
      df_baseline <- df_baseline %>% filter(cap == 8e+05)
    }
    
    df_baseline <- df_baseline %>%
      select(-c(env, cap, wgts))
    
    # Calculate ratios based on comparison type
    if(by_env) {
      df_ratios <- df_ny %>%
        filter(env != "NoClimate", cap == 8e+05, wgts == "equal")
    } else if(by_cap) {
      df_ratios <- df_ny %>%
        filter(env == "NoClimate", cap != 8e+05, wgts == "equal")
    } else {
      df_ratios <- df_ny %>%
        filter(env == "NoClimate", cap == 4e+05, wgts != "equal")
    }
    
    df_ratios <- df_ratios %>%
      left_join(df_baseline, by = c("Name", "LongName", "age")) %>%
      mutate(ratio = value_mean / value_baseline,
             percent_change = ((value_mean - value_baseline) / value_baseline) * 100,
             metric = metric_name)
    
    return(df_ratios)
  }
  
  # Process data and combine
  ratio_list <- list()
  
  if(plot_type %in% c("WAA", "both")) {
    waa_ratios <- process_age_data(waa_df, "weight", "Weight-at-Age")
    if(nrow(waa_ratios) > 0) {
      ratio_list$waa <- waa_ratios
    } else {
      message("No WAA data to plot after filtering")
    }
  }
  
  if(plot_type %in% c("NAA", "both")) {
    naa_ratios <- process_age_data(naa_df, "abun", "Numbers-at-Age")
    if(nrow(naa_ratios) > 0) {
      ratio_list$naa <- naa_ratios
    } else {
      message("No NAA data to plot after filtering")
    }
  }
  
  # Check if we have any data to plot
  if(length(ratio_list) == 0) {
    message("No data to plot after filtering")
    return(NULL)
  }
  
  # Combine the data frames
  combined_ratios <- bind_rows(ratio_list)
  
  # order caps
  combined_ratios$cap <- factor(combined_ratios$cap, levels = c("6e+05", "4e+05", "2e+05"))
  
  # Create the plot
  p <- combined_ratios %>%
    ggplot() +
    geom_tile(aes(x = age, y = LongName, fill = percent_change), color = 'darkgrey') +
    colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                                 oob = scales::squish,
                                                 breaks = scales::pretty_breaks(n = 6)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(1, max(combined_ratios$age, na.rm = TRUE))) +
    labs(x = 'Age class', 
         y = '', 
         fill = '% change\nfrom\nbaseline') +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10))
  
  # Add faceting based on whether we're plotting one or both metrics
  if(plot_type == "both" && length(ratio_list) == 2) {
    # Facet by metric (rows) and comparison variable (columns)
    p <- p + facet_grid(reformulate(facet_var, "metric"))
    plot_height <- 4
    filename <- paste0(plotdir, "/heatmaps/age_heatmap_", lab, "_", tag, ".png")
  } else {
    # Facet only by comparison variable
    p <- p + facet_grid(reformulate(facet_var, "."))
    plot_height <- 3
    
    # Determine filename based on which metric we actually plotted
    if(!is.null(ratio_list$waa)) {
      filename <- paste0(plotdir, "/heatmaps/waa_heatmap_", lab, "_", tag, ".png")
    } else {
      filename <- paste0(plotdir, "/heatmaps/naa_heatmap_", lab, "_", tag, ".png")
    }
  }
  
  # Save the plot
  ggsave(filename, 
         p, 
         width = 10, 
         height = plot_height, 
         units = "in", 
         dpi = 300)
  
  # Return the combined data for inspection if needed
  return(p)
  invisible(combined_ratios)
}

#' Create a heatmap of predator diet composition changes relative to baseline
#'
#' This function creates heatmap visualizations showing percentage changes in
#' diet composition compared to a baseline scenario (NoClimate with equal weights).
#' The comparison can be made across environmental scenarios, capacity levels, or
#' weight schemes.
#'
#' @param diet_df Data frame containing diet composition data with columns: Predator,
#'   Cohort, is_oy, comp_agg, Name, env, cap, wgts
#' @param predator Character. Optional predator code to filter to a single predator.
#'   If NULL (default), plots all predators with predator as a facet row
#' @param by_env Logical. If TRUE, plot comparison across environmental scenarios.
#'   Default is FALSE
#' @param by_cap Logical. If TRUE, plot comparison across capacity levels.
#'   Default is FALSE
#' @param by_wgts Logical. If TRUE, plot comparison across weight schemes.
#'   Default is FALSE
#' @param tag Character string. Optional tag to append to the output filename.
#'   Default is NULL
#'
#' @return Invisibly returns a data frame of diet composition ratios with columns
#'   including percent_change and ratio. If no data remains after filtering,
#'   returns NULL with a message
#'
#' @details
#' Exactly one of by_env, by_cap, or by_wgts must be TRUE. The function:
#' \itemize{
#'   \item Filters data to is_oy == 1
#'   \item Compares to baseline: NoClimate scenario, cap=8e+05 (or 4e+05 for by_wgts),
#'         and equal weights
#'   \item Creates a diverging color heatmap (blue-white-red) showing percent change
#'   \item Saves plot to plotdir as PNG (10x height, 300 dpi)
#' }
#'
#' @note Requires ggplot2, dplyr, and colorspace packages. Assumes plotdir is
#'   defined in the global environment
#'
#' @examples
#' \dontrun{
#' # Plot all predators by environmental scenario
#' plot_diet_heatmap(pred_diets, by_env = TRUE, tag = "climate")
#'
#' # Plot single predator by capacity level
#' plot_diet_heatmap(pred_diets, predator = "BDF", by_cap = TRUE, tag = "capacity")
#' }
#'
#' @export
plot_diet_heatmap <- function(diet_df,
                              by_env = FALSE,
                              by_cap = FALSE,
                              by_wgts = FALSE,
                              tag = NULL) {
  
  # Validate inputs
  if(sum(by_env, by_cap, by_wgts) != 1) {
    stop("You need to plot either by env, cap, or wgts")
  }
  
  # Determine label for faceting and filename
  if(by_env) {
    lab <- "env"
    facet_var <- "env"
  } else if(by_cap) {
    lab <- "cap"
    facet_var <- "cap"
  } else {
    lab <- "wgts"
    facet_var <- "wgts"
  }
  
  # Filter to is_oy == 1
  diet_filtered <- diet_df %>%
    filter(is_oy == 1)
  
  # Create baseline
  # Base conditions: env = NoClimate, cap = 8e+05, wgts = equal
  diet_baseline <- diet_filtered %>% 
    filter(env == "NoClimate", wgts == "equal") %>%
    rename(comp_baseline = comp_agg)
  
  if(by_wgts) {
    diet_baseline <- diet_baseline %>% filter(cap == 4e+05)
  } else {
    diet_baseline <- diet_baseline %>% filter(cap == 8e+05)
  }
  
  # Drop unnecessary cols
  diet_baseline <- diet_baseline %>%
    select(-c(env, cap, wgts, is_oy, run))
  
  # Calculate ratios based on comparison type
  if(by_env) {
    diet_ratios <- diet_filtered %>%
      filter(env != "NoClimate", cap == 8e+05, wgts == "equal")
  } else if(by_cap) {
    diet_ratios <- diet_filtered %>%
      filter(env == "NoClimate", cap != 8e+05, wgts == "equal")
  } else {
    diet_ratios <- diet_filtered %>%
      filter(env == "NoClimate", cap == 4e+05, wgts != "equal")
  }
  
  diet_ratios <- diet_ratios %>%
    left_join(diet_baseline, by = c("Predator", "Age class", "LongName")) %>%
    mutate(ratio = comp_agg / comp_baseline,
           percent_change = ((comp_agg - comp_baseline) / comp_baseline) * 100)
  
  # Check if we have data
  if(nrow(diet_ratios) == 0) {
    message("No data to plot after filtering")
    return(NULL)
  }
  
  # order caps
  diet_ratios$cap <- factor(diet_ratios$cap, levels = c("6e+05", "4e+05", "2e+05"))
  
  # Create heatmap
  p <- diet_ratios %>%
    ggplot() +
    geom_tile(aes(x = `Age class`, y = LongName, fill = percent_change), color = 'darkgrey') +
    colorspace::scale_fill_continuous_divergingx(palette = 'PRGn', mid = 0, rev = F,
                                                 oob = scales::squish,
                                                 breaks = scales::pretty_breaks(n = 6)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(1, max(diet_ratios$`Age class`, na.rm = TRUE))) +
    labs(x = 'Age class', 
         y = '', 
         fill = '% change\nfrom\nbaseline\ngroundfish\nin diet') +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10))
  
  p <- p + facet_grid(reformulate(facet_var, "."))
  
  # Determine filename
  pred_tag <- if(!is.null(predator)) paste0("_", predator) else ""
  filename <- paste0(plotdir, "/heatmaps/diet_heatmap_", lab, pred_tag, "_", tag, ".png")
  
  # Save the plot
  ggsave(filename, 
         p, 
         width = 10, 
         height = 3, 
         units = "in", 
         dpi = 300)
  
  # Return the data for inspection
  return(p)
  invisible(diet_ratios)
}
