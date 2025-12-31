################################################################################
# Create Population Trend Scenarios for Decision Model
################################################################################
#
# Purpose: 
#   Generate Pbar (poults per hen at inflection) values for three population
#   trend scenarios (Increase, Stable, Decrease) using:
#   1. Observed median Pbar from PA WMUs as "stable" baseline
#   2. ±10% adjustments for increase/decrease scenarios
#   3. Weather model refinements for April decision timing
#
# Rationale:
#   - Stable = median observed Pbar across PA regions
#   - Increase = +10% from stable (favorable conditions)
#   - Decrease = -10% from stable (unfavorable conditions)
#   - April values refined by temperature-recruitment relationship
#
# Input:  
#   region_parameters_for_mdp.csv (BBS-derived parameters, script 02)
#   april_scaled_weather2.rds (April weather data, script 04)
#   Fitted weather-recruitment model (script 05)
#
# Output:
#   Scenario-specific parameter files for MATLAB MDP model
#
# Author: Veronica A. Winter
# Date: December 2025
################################################################################

rm(list = ls())
gc()

library(glmmTMB)
library(dplyr)

################################################################################
# 1. CALCULATE BASELINE "STABLE" PBAR FROM OBSERVED DATA
################################################################################

cat("\n========== CALCULATING BASELINE PBAR ==========\n")

# Load observed Pbar values from BBS analysis (script 02)
# These are empirically-derived inflection point parameters
region_params <- read.csv('Data/region_parameters_for_mdp.csv')

# Calculate median Pbar across all PA WMU regions
# This represents "average" or "stable" conditions
Pbar_stable_base <- median(region_params$Pbar, na.rm = TRUE)
range(region_params$Pbar)

cat("Observed Pbar values across regions:\n")
print(region_params %>% select(Region, Pbar))

cat("\nMedian Pbar (STABLE baseline):", round(Pbar_stable_base, 3), "\n")

################################################################################
# 2. CREATE ±10% SCENARIOS FOR INCREASE/DECREASE
################################################################################

cat("\n========== CREATING POPULATION TREND SCENARIOS ==========\n")

# Stable = median observed (baseline)
# Because median density is high, reducing the baseline to a range that matches 
# more with the literature
Pbar_stable <- 1.8

# Increase = +10% from stable (favorable reproduction)
Pbar_increase <- Pbar_stable + (Pbar_stable * 0.10)

# Decrease = -10% from stable (poor reproduction)
Pbar_decrease <- Pbar_stable - (Pbar_stable * 0.10)

cat("Population trend scenarios (±10% from median):\n")
cat("  Stable (baseline):", round(Pbar_stable, 3), "\n")
cat("  Increase (+18%):", round(Pbar_increase, 3), "\n")
cat("  Decrease (-18%):", round(Pbar_decrease, 3), "\n")

# Similarly adjust Fbar (female density at inflection)
Fbar_stable_base <- median(region_params$Fbar, na.rm = TRUE)
# Since Fbar has a wide range, choose value that would coincide with PBar = 1.8
Fbar_stable <- 2.5
Fbar_increase <- Fbar_stable + (Fbar_stable * 0.10)
Fbar_decrease <- Fbar_stable - (Fbar_stable * 0.10)

cat("\nFemale density at inflection (Fbar):\n")
cat("  Stable:", round(Fbar_stable, 3), "\n")
cat("  Increase:", round(Fbar_increase, 3), "\n")
cat("  Decrease:", round(Fbar_decrease, 3), "\n")

################################################################################
# 3. FIT WEATHER-RECRUITMENT MODEL FOR APRIL REFINEMENT
################################################################################

cat("\n========== FITTING WEATHER-RECRUITMENT MODEL ==========\n")

# Load April weather data
month <- "april"
scaled_weather <- readRDS(paste0("Data/Rec_data/", month, "_scaled_weather2.rds"))

# Load observed PPB data
ph_df <- readRDS("Data/Rec_data/ph_df_aug31.rds")

# Join data
ph_train <- ph_df %>%
  left_join(scaled_weather, by = c("Year", "MU"))

# Fit model: PPB ~ temperature
ph_formula <- PHratio ~ 0 + scale_avg_temperature + (1|MU)
m.refined.ph <- glmmTMB(ph_formula, data = ph_train, 
                        family = gaussian(link = "identity"))

# Extract temperature coefficient
temp_coef <- fixef(m.refined.ph)$cond["scale_avg_temperature"]

cat("Temperature effect on PPB:", round(temp_coef, 4), "\n")
# Temperature effect on PPB: 0.0427
################################################################################
# 4. DEFINE TEMPERATURE SCENARIOS FOR APRIL
################################################################################

cat("\n========== DEFINING TEMPERATURE SCENARIOS ==========\n")

# Calculate temperature quantiles from historical data
temp_quantiles <- quantile(ph_train$scale_avg_temperature, 
                           probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

temp_cold <- temp_quantiles[1]   # 10th percentile (cold spring)
temp_avg <- temp_quantiles[2]    # 50th percentile (average spring)
temp_warm <- temp_quantiles[3]   # 90th percentile (warm spring)

cat("Temperature scenarios (standardized):\n")
cat("  Cold spring (10th percentile):", round(temp_cold, 3), "\n")
cat("  Average spring (median):", round(temp_avg, 3), "\n")
cat("  Warm spring (90th percentile):", round(temp_warm, 3), "\n")

################################################################################
# 5. CALCULATE WEATHER ADJUSTMENTS TO PBAR
################################################################################

cat("\n========== CALCULATING WEATHER ADJUSTMENTS ==========\n")

# Calculate PPB effect of cold/warm spring relative to average
# Effect = coefficient × (scenario_temp - average_temp)
effect_cold <- temp_coef * (temp_cold - temp_avg)
effect_warm <- temp_coef * (temp_warm - temp_avg)

cat("Weather model effects on PPB:\n")
cat("  Cold spring effect:", round(effect_cold, 4), "\n")
cat("  Warm spring effect:", round(effect_warm, 4), "\n")

# Apply adjustments to baseline Pbar values
# For April, we're refining the estimate based on observed spring weather

# Cold April (unfavorable weather)
Pbar_stable_cold <- Pbar_stable + effect_cold
Pbar_increase_cold <- Pbar_increase + effect_cold
Pbar_decrease_cold <- Pbar_decrease + effect_cold

# Warm April (favorable weather)
Pbar_stable_warm <- Pbar_stable + effect_warm
Pbar_increase_warm <- Pbar_increase + effect_warm
Pbar_decrease_warm <- Pbar_decrease + effect_warm

################################################################################
# 6. SUMMARY OF ALL SCENARIO VALUES
################################################################################

cat("\n========== SUMMARY OF SCENARIO PARAMETERS ==========\n")

# Create summary table
scenario_table <- data.frame(
  Decision_Timing = rep(c("Jan/Sept", "April_Cold", "April_Warm"), each = 3),
  Population_Trend = rep(c("Decrease", "Stable", "Increase"), 3),
  Fbar = c(
    # Jan/Sept (no weather info)
    Fbar_decrease, Fbar_stable, Fbar_increase,
    # April Cold (weather info: cold)
    Fbar_decrease, Fbar_stable, Fbar_increase,
    # April Warm (weather info: warm)
    Fbar_decrease, Fbar_stable, Fbar_increase
  ),
  Pbar = c(
    # Jan/Sept (no weather info)
    Pbar_decrease, Pbar_stable, Pbar_increase,
    # April Cold (cold spring observed)
    Pbar_decrease_cold, Pbar_stable_cold, Pbar_increase_cold,
    # April Warm (warm spring observed)
    Pbar_decrease_warm, Pbar_stable_warm, Pbar_increase_warm
  )
) %>%
  mutate(
    Fbar = round(Fbar, 3),
    Pbar = round(Pbar, 3)
  )

print(scenario_table)

# Save summary table
write.csv(scenario_table, "Data/scenario_parameter_summary.csv", row.names = FALSE)

################################################################################
# 7. INTERPRET ADJUSTMENTS
################################################################################

cat("\n========== INTERPRETATION ==========\n")

cat("Baseline (Jan/Sept, no weather info):\n")
cat("  Stable: Pbar =", round(Pbar_stable, 3), "(median observed)\n")
cat("  Increase: Pbar =", round(Pbar_increase, 3), "(+10%)\n")
cat("  Decrease: Pbar =", round(Pbar_decrease, 3), "(-10%)\n\n")

cat("April Cold Spring adjustments:\n")
cat("  Adjustment:", round(effect_cold, 4), "\n")
cat("  Stable: Pbar =", round(Pbar_stable_cold, 3), "\n")
cat("  Increase: Pbar =", round(Pbar_increase_cold, 3), "\n")
cat("  Decrease: Pbar =", round(Pbar_decrease_cold, 3), "\n\n")

cat("April Warm Spring adjustments:\n")
cat("  Adjustment:", round(effect_warm, 4), "\n")
cat("  Stable: Pbar =", round(Pbar_stable_warm, 3), "\n")
cat("  Increase: Pbar =", round(Pbar_increase_warm, 3), "\n")
cat("  Decrease: Pbar =", round(Pbar_decrease_warm, 3), "\n\n")

################################################################################
# 8. CREATE PARAMETER FILES FOR MATLAB
################################################################################

cat("\n========== CREATING PARAMETER FILES FOR MDP ==========\n")

# Load other required parameters
params <- readRDS("Data/LogisticGrowthModParams_2023.rds")

# Set slope to be the same in all scenarios
range(params$r)
# Slope will be 0.2

################################################################################
# END OF SCRIPT
################################################################################

cat("\n========== COMPLETE ==========\n")
cat("All scenario parameter files created!\n")
cat("Ready for MATLAB MDP analysis\n\n")
