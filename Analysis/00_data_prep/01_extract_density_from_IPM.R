################################################################################
# Extract Population Density Parameters from Integrated Population Model
################################################################################
# 
# Purpose: Extract female density, poults per brood (ppb), and hen-with-brood 
#          (hwb) estimates from the integrated population model (IPM) output
#          (Winter et al., in review). These values are used to scale BBS 
#          indices to absolute densities.
#
# Input:  O_24_abundance_summary.rds - IPM abundance estimates by WMU and year
#         O_24_ppb_summary.rds - Poults per brood estimates
#         O_24_hwb_summary.rds - Hen-with-brood ratio estimates
#
# Output: Female density values by WMU region for 2024
#         ppb and hwb summary statistics
#
# Author: Veronica A. Winter
# Date: December 2025
################################################################################
rm(list=ls())
gc()

library(dplyr)

################################################################################
# 1. Load IPM Abundance Data
################################################################################

# Load abundance summary from IPM (contains density estimates by sex, year, WMU)
x <- readRDS("Data/IPM/Winter_etal_IPM_abundance_summary.rds")
head(x)

################################################################################
# 2. Calculate Total Abundance by Year and WMU
################################################################################

# Sum across all sex/age classes to get total abundance
abundance <- x %>% 
  mutate(year = case_when(
    year == "1" ~ "2020",
    year == "2" ~ "2021",
    year == "3" ~ "2022",
    year == "4" ~ "2023",
    year == "5" ~ "2024"
  )) %>% 
  group_by(year, wmu) %>% 
  reframe(total = sum(density_value))

################################################################################
# 3. Extract Female Density for Fall 2024 (Key Parameter)
################################################################################

# Female density is used as the scaling factor to convert BBS indices to 
# absolute densities in the logistic growth model
abundance_fall <- x %>% 
  filter(sex == "Female") %>%  # Only adult females
  mutate(year = case_when(
    year == "1" ~ "2020",
    year == "2" ~ "2021",
    year == "3" ~ "2022",
    year == "4" ~ "2023",
    year == "5" ~ "2024"
  )) %>% 
  mutate(WMU_Group = paste("Region", wmu, sep = " ")) %>% 
  group_by(year, WMU_Group) %>% 
  filter(year == "2024") %>%  # Most recent year
  reframe(total = sum(density_value)) %>% 
  arrange(WMU_Group)

# Median female density across all regions
median(abundance_fall$total)
# Result: 1.68 females per km^2

################################################################################
# 4. Extract Poults Per Brood (ppb)
################################################################################

# ppb represents reproductive output and is used in the fecundity function
ppb <- readRDS("Data/IPM/Winter_etal_IPM_ppb_summary.rds")

ppb_overall <- ppb %>% 
  mutate(year = case_when(
    year == "1" ~ "2020",
    year == "2" ~ "2021",
    year == "3" ~ "2022",
    year == "4" ~ "2023",
    year == "5" ~ "2024"
  )) %>%  
  group_by(year, wmu) %>% 
  reframe(total = sum(median_value)) %>% 
  filter(year == "2024")

median(ppb_overall$total)

# Result: 3.75 poults per brood

################################################################################
# 5. Extract Hen-With-Brood Ratio (hwb)
################################################################################

# hwb is the proportion of hens successfully raising broods
# Combined with ppb, this gives poult density: P = F * hwb * ppb
hwb <- readRDS("Data/IPM/Winter_etal_IPM_hwb_summary.rds")

hwb_overall <- hwb %>% 
  mutate(year = case_when(
    year == "1" ~ "2020",
    year == "2" ~ "2021",
    year == "3" ~ "2022",
    year == "4" ~ "2023",
    year == "5" ~ "2024"
  )) %>% 
  group_by(year, wmu) %>% 
  reframe(total = sum(median_value)) %>% 
  filter(year == "2024")

median(hwb_overall$total)
# Result: 0.77 hen-with-brood ratio (by WMU)

################################################################################
# Notes for Next Steps
################################################################################

# These IPM-derived density values are used in:
#   02_BBS_logistic_growth.R - to scale BBS indices to absolute densities
#   03_calculate_inflection_points.R - to derive F̄ and P̄ parameters
#
# Key relationship:
#   Female density at inflection point (F̄) = K_density / 2
#   where K_density is scaled using abundance_fall values
################################################################################