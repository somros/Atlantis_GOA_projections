# Alberto Rovellini
# 12/19/2024
# Script to explore TAC/ABC ratio in the EBS and GOA
# Based on harvest specification data from AKFIN
library(tidyverse)
library(readxl)
library(viridis)

# Gulf of Alaska ----------------------------------------------------------
# read in GOA
goa_specs <- read_xlsx("data/GOA_harvest specs_1986-2024.xlsx", 
                       sheet = 1, 
                       range = "A3:DO131",
                       col_names = F) 

# set names for columns
yrs <- 2024:1986
tags <- paste(c("OFL","ABC","TAC"),rep(yrs,each=3),sep="_")
colnames(goa_specs) <- c("Stock","Area",tags)

# need to clean up all the asterisk and comma garbage first
goa_specs <- goa_specs %>%
  # Convert columns to numeric, handling special cases
  mutate(across(!c(Stock,Area), ~{
    x <- ifelse(tolower(.) %in% c("n/a", "n/a"), NA, .)
    # First remove asterisks
    x <- gsub("\\*", "", x)
    # Then remove commas
    x <- gsub(",", "", x)
    # Convert to numeric
    as.numeric(x)
  }))

# fill stock names
goa_specs <- goa_specs %>%
  fill(Stock, .direction = "down")

# melt
goa_specs_long <- goa_specs %>%
  pivot_longer(-c(Stock,Area), names_to = "Var_Year", values_to = "mt") %>%
  separate(col = Var_Year,
           into = c("Var", "Year"),
           sep = "_") %>%
  mutate(Year = as.numeric(Year))

# keep Total GOA only, and keep ABC and TAC only
goa_specs_tot <- goa_specs_long %>% 
  filter(Area %in% c("Total", "Total (GW)", "GW"), Var %in% c("ABC","TAC","OFL")) %>%
  drop_na() %>%
  mutate(Area = "GOA") %>%
  group_by(Stock,Area,Var,Year)%>%
  summarise(mt = sum(mt))

# plot for AMSS
tac_sp <- unique(goa_specs_tot$Stock)
length(unique(goa_specs_tot$Stock))
colors <- c(viridis(11)[2:10], inferno(11)[2:10], cividis(10)[2:9])
goa_specs_tot$Var <- factor(goa_specs_tot$Var, levels = c("OFL","ABC","TAC"))
p <- goa_specs_tot %>%
  filter(Year >1991) %>%
  filter(Var != "OFL") %>%
  ggplot(aes(x = Year, y = mt, fill = Stock))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = 800000, linetype = "dashed", color = "red")+
  theme_bw()+
  labs(x = "", y = "Catch (mt)", fill = "")+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'))+
  guides(fill = guide_legend(nrow = 7))+
  facet_grid(~Var)
p

# add catch for AMSS
catch_data <- read.csv("data/Groundfish Total Catch.csv", fileEncoding = 'UTF-8-BOM')

# there are a lot of non TAC species reported here, as well as by-catch species, species that are in the FMP but are not groundfish, etc.
# For the purpose of comparing to ABC/TAC plots, we will filter only the species that have a TAC in the harvest specifications
# Map species in the catch to species in the harvest specification data set
tac_key <- read.csv("data/tac_catch_key.csv", header = T)

# process the catch data so that it can be mapped to the harvest specification data
catch_data_short <- catch_data %>%
  select(Year, Species.Group.Name, Catch..mt.) %>%
  left_join(tac_key, by = c("Species.Group.Name" = "Catch_sp")) %>%
  group_by(Year, TAC_sp) %>%
  summarise(mt = sum(Catch..mt., na.rm = T)) %>%
  mutate(Area = "GOA", Var = "Catch") %>%
  select(TAC_sp, Area, Var, Year, mt) %>%
  rename(Stock = TAC_sp)

goa_specs_tot_catch <- goa_specs_tot %>% rbind(catch_data_short) %>% drop_na()

# order factors
goa_specs_tot_catch$Var <- factor(goa_specs_tot_catch$Var, levels = c("OFL", "ABC", "TAC", "Catch"))

p <- goa_specs_tot_catch %>%
  filter(Year >1991) %>%
  filter(Var %in% c("TAC","Catch")) %>%
  ggplot(aes(x = Year, y = mt, fill = Stock))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = 800000, linetype = "dashed", color = "red")+
  theme_bw()+
  labs(x = "", y = "Catch (mt)", fill = "")+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'))+
  guides(fill = guide_legend(nrow = 7))+
  facet_grid(~Var)
p

# pivot wider
goa_ratio <- goa_specs_tot %>%
  pivot_wider(id_cols = c(Stock,Area,Year), names_from = Var, values_from = mt) %>%
  mutate(ratio = TAC/ABC)

# view
goa_ratio %>%
  ggplot(aes(x = Year, y = ratio))+
  geom_line()+
  facet_wrap(~Stock)


# Plot specs for GOA ------------------------------------------------------
goa_specs_tot_catch %>%
  filter(Var != "OFL") %>%
  filter(!(Stock %in% c("Octopus","Squids","Sharks","Other Species"))) %>%
  filter(Year >= 2010) %>%
  ggplot(aes(x = Year, y = mt, color = Var))+
  geom_line(size = 1, alpha = 0.5)+
  #scale_color_manual(values = c("#2E86AB", "#F24236"))+
  facet_wrap(~Stock, scales = 'free')

# Determine the OY weights ----------------------------------------
# assume that attainment (catch/TAC) is a measure of how much a stock is valued, i.e. it reflects fisher's behavior OR it reflects bycatch patterns of choke species
# First need to map species in this data to the Atlantis groups
# Should this be catch/TAC? TAC/ABC? Catch/ABC?
# Catch/TAC is a measure of a species's attainment
# In the EBS, OY acts when going from ABC to TAC, but this is the GOA and it may not work so well

grp <- read.csv("data/GOA_Groups.csv")

key <- data.frame("Stock" = unique(goa_specs_tot_catch$Stock),
                  "Code" = c("ATF","DFS","SKB","FFD","RFP","FHS","SKL","RFS","OCT","FFS","RFD","SKO",NA,"POP","COD","RFP","POL","REX","RFS","SBF","SCU","FFS",NA,"RFS","RFS","SQD","THO"))

# clean stock names and map to Atlantis stocks
goa_specs_clean <- goa_specs_tot_catch %>%
  filter(!(Stock %in% c("Octopus","Squids","Sharks","Other Species"))) %>%
  left_join(key, by = "Stock") %>%
  group_by(Code,Area,Var,Year) %>%
  summarise(mt = sum(mt)) %>%
  ungroup() %>%
  left_join(grp %>% select(Code,Name), by = "Code") %>%
  filter(Year >= 2010, Year < 2020) %>%
  group_by(Code,Name,Area,Var) %>%
  summarize(mean_mt = mean(mt, na.rm = T)) %>%
  ungroup() 

# look at attainment first
catch_to_tac <- goa_specs_clean %>%
  filter(Var %in% c("TAC","Catch")) %>%
  pivot_wider(names_from = Var, values_from = mean_mt) %>%
  mutate(ratio = Catch / TAC) %>%
  arrange(-ratio) %>%
  mutate(w = n():1)

write.csv(catch_to_tac, "data/w_from_harvest_specs.csv", row.names = F)

# alternatively, do TAC/ABC
# in principle this is what is done in the BSAI under the OY
# tac_to_abc <- goa_specs_clean %>%
#   filter(Var %in% c("TAC","ABC")) %>%
#   pivot_wider(names_from = Var, values_from = mean_mt) %>%
#   mutate(ratio = TAC / ABC) %>%
#   arrange(-ratio) %>%
#   mutate(w = n():1)
# in reality, it ends up being pretty odd - deep water flatfish ends up being some of the highest "value" stocks 
# part of the reason is that TAC=ABC for many stocks, including non-target stocks
# this is because there is no constraining cap in the GOA at present so managers tend to not worry about rescaling ABC for much of the flatfish - attainment is so low anyway
# so use the attainment metric instead


# two designs: a ramp and a binary scheme for w
# Ramp: weights are a simple progression from least to most valuable, or a linear transformation thereof
ramp <- catch_to_tac %>%
  arrange(-w) 
ramp

# tt <- ramp$ratio
# max(tt - c(tt[2:length(tt)],0))
#   
# see it here:
# A tibble: 18 × 3
# Code  Name                        w
# <chr> <chr>                   <int>
# 1 COD   Cod                        18
# 2 SBF   Sablefish                  17
# 3 POP   Pacific_ocean_perch        16
# 4 POL   Pollock                    15
# 5 SKO   Skate_other                14
# 6 RFP   Rockfish_pelagic_shelf     13
# 7 RFS   Rockfish_slope             12
# 8 DFS   Shallow_demersal           11
# 9 SKB   Skate_big                  10
# 10 RFD   Rockfish_demersal_shelf     9
# 11 THO   Thornyhead                  8
# 12 SKL   Skate_longnose              7
# 13 ATF   Arrowtooth_flounder         6
# 14 REX   Rex_sole                    5
# 15 SCU   Sculpins                    4
# 16 FFS   Flatfish_shallow            3
# 17 FHS   Flathead_sole               2
# 18 FFD   Flatfish_deep               1


# print out for prm
wgts_print <- grp %>%
  select(Code) %>%
  left_join(ramp, by = "Code") %>%
  select(Code,w) %>%
  mutate(w = replace_na(w,0))

# Binary: start from ramp and land somewhere for the demarcation
# it was agreed upon that starting arrowtooth and below would be low value
# It was suggested to use a ternary scheme instead: commercial species, choke species, and non-target species
# for now I'd keep it simple and cut between skate longnose and atf
binary <- ramp %>%
  mutate(w = ifelse(w > 6, 10, 1)) # what values should we actually use here? use the shiny app

# cutting here makes sense because:
# just based on the ratio column, there is a gap between SKL and ATF
# this way all choke species (skates) have 10; also all those expensive rockfish in the ex-vessel price data have 10
# based on TAC the split is roughly 60/40 for species of value VS no value
binary %>%
  group_by(w) %>%
  summarise(tot_TAC = sum(TAC))
