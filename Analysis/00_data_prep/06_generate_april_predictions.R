################################################################################
# Generate April Recruitment Predictions for Decision Model
################################################################################
#
# Purpose: 
#   Use weather-recruitment model to predict Pbar (poults per hen at 
#   inflection point) in April, before recruitment is directly observed
#   in September. This enables earlier decision-making with updated
#   biological information.
#
# Method:
#   1. Load fitted weather-recruitment model (from script 05)
#   2. Predict PPB for April 2023 based on spring weather
#   3. Calculate updated Pbar = PPB × HWB
#   4. Create region-specific parameter sets for April MDP model
#
# Input:  
#   april23_scaled_weather.rds (April 2023 weather data)
#   m.refined.ppb (fitted model from 05_fit_weather_recruitment_model.R)
#   LogisticGrowthModParams_2023.rds (BBS parameters)
#   0a_findingDensity.R (IPM density estimates)
#
# Output:
#   rec_weather_prediction.csv (April-predicted parameters)
#   params_ppb_summary_april2023.csv (complete parameter table)
#   Comparison with September parameters
#
# Author: Veronica A. Winter
# Date: 2025
################################################################################

rm(list = ls())
gc()

library(glmmTMB)
library(dplyr)

################################################################################
# 1. LOAD APRIL WEATHER DATA
################################################################################

cat("\n========== LOADING APRIL WEATHER DATA ==========\n")

# Load April 2023 weather data
# This is scaled/weighted weather from script 04
april23_weather_df <- readRDS("Data/Rec_data/april23_scaled_weather.rds")

cat("April 2023 weather data loaded\n")
cat("  Groups:", paste(unique(april23_weather_df$Group), collapse = ", "), "\n")
cat("  Date range:", min(april23_weather_df$date), "to", max(april23_weather_df$date), "\n")

################################################################################
# 2. LOAD FITTED WEATHER-RECRUITMENT MODEL
################################################################################

cat("\n========== FITTING PPB MODEL ==========\n")

# Specify weather data to use
month <- "april"
scaled_weather <- readRDS(paste0("Data/Rec_data/", month, "_scaled_weather2.rds"))

# Load observed PPB data
ppb_df <- readRDS("Data/Rec_data/ph_df_aug31.rds")

# Join with weather
ppb <- ppb_df %>%
  left_join(scaled_weather, by = c("Year", "MU"))

# Define and fit model (same as script 05)
ppb_formula <- PHratio ~ 0 + scale_avg_temperature + (1|MU)

m.refined.ppb <- glmmTMB(ppb_formula, data = ppb, 
                         family = gaussian(link = "identity"))

cat("Model fitted:\n")
print(summary(m.refined.ppb))

################################################################################
# 3. PREDICT PPB FOR APRIL 2023
################################################################################

cat("\n========== PREDICTING PPB FOR APRIL 2023 ==========\n")

# Use fitted model to predict PPB based on April weather
april_preds <- predict(
  m.refined.ppb,
  newdata = april23_weather_df,
  type = "response",
  se.fit = TRUE,
  re.form = NULL  # Include WMU random effects
)

# Attach predictions to data
april23_weather_df$PPB_pred <- as.numeric(april_preds$fit)
april23_weather_df$PPB_se <- april_preds$se.fit
april23_weather_df$PPB_lower <- april_preds$fit - 1.96 * april_preds$se.fit
april23_weather_df$PPB_upper <- april_preds$fit + 1.96 * april_preds$se.fit

# Calculate mean PPB per region (average across all days in April)
mean_april_ppb <- april23_weather_df %>% 
  filter(!Group == "Group 10") %>%  # Exclude Group 10 (insufficient data)
  group_by(Group) %>% 
  summarise(PPB = mean(PPB_pred, na.rm = TRUE), .groups = 'drop')

cat("April PPB predictions:\n")
print(mean_april_ppb)

################################################################################
# 4. LOAD DENSITY AND HWB DATA
################################################################################

cat("\n========== LOADING ADDITIONAL PARAMETERS ==========\n")

# Load density estimates from IPM
if (file.exists("Analysis/01_extract_density_from_IPM.R")) {
  message("Loading IPM density estimates...")
  source("Analysis/01_extract_density_from_IPM.R")
  
  k_density <- abundance_fall %>% 
    select(WMU_Group, total)
  
  hwb_ratio <- hwb_overall$total  # Hen-with-brood ratio
  ppb_value <- mean_april_ppb$PPB  # April-predicted PPB
  
} else {
  stop("ERROR: Cannot find density estimation script.")
}

cat("Parameters loaded:\n")
cat("  HWB ratio:", round(mean(hwb_ratio), 3), "\n")
cat("  Mean predicted PPB:", round(mean(ppb_value), 3), "\n")

################################################################################
# 5. CALCULATE APRIL-BASED INFLECTION POINT PARAMETERS
################################################################################

cat("\n========== CALCULATING APRIL PARAMETERS ==========\n")

# Load BBS logistic growth parameters
params <- readRDS("Data/LogisticGrowthModParams_2023.rds")

# Join with density data
params_df <- as.data.frame(params) %>% 
  mutate(WMU_Group = rownames(params)) %>% 
  left_join(k_density, by = "WMU_Group") %>% 
  rename(K_density = total)

# Fill missing values with mean
params_df$K_density <- ifelse(is.na(params_df$K_density), 
                              mean(params_df$K_density, na.rm = TRUE), 
                              params_df$K_density)

# Calculate inflection point parameters
# Key difference from September: using April-predicted PPB
params_df <- params_df %>%
  mutate(
    # Scaling factor: BBS index to density
    scaling_factor = K_density / K,
    
    # Female density at inflection (F̄)
    N_star = (K / 2) * scaling_factor,
    
    # Region identifier
    Region = WMU_Group,
    
    # April-predicted poults per brood
    ppb = ppb_value,
    
    # Poult density at inflection (P̄)
    # P̄ = PPB × HWB (April prediction × observed HWB)
    Pbar = ppb * hwb_ratio
  )

# Save April parameter table
write.csv(params_df, "Data/params_ppb_summary_april2023.csv", row.names = FALSE)

cat("✓ April parameters calculated\n")

################################################################################
# 6. GENERATE DENSITY-DEPENDENT CURVES (OPTIONAL VISUALIZATION)
################################################################################

# Load helper function
source("Analysis/00_data_prep/xx_pph_density_fun.R")

# Create sequence for plotting
density_seq <- seq(0, max(params_df$K_density, na.rm = TRUE) * 1.2, 
                   length.out = 100)

# Generate PPH curves for each region
pph_curves <- data.frame()

for (i in 1:nrow(params_df)) {
  region_curve <- data.frame(
    WMU_Group = params_df$WMU_Group[i],
    female_density = density_seq,
    pph = calculate_pph(
      female_density = density_seq,
      N_star = params_df$N_star[i],
      Pbar = params_df$Pbar[i],
      slope = params_df$r[i]
    )
  )
  pph_curves <- rbind(pph_curves, region_curve)
}

################################################################################
# 7. CREATE REGION-SPECIFIC PARAMETER SETS FOR MDP MODEL
################################################################################

cat("\n========== CREATING MDP PARAMETER SETS ==========\n")

regionParams <- list()

for (i in 1:nrow(params_df)) {
  region_name <- params_df$WMU_Group[i]
  regionParams[[region_name]] <- list(
    Fbar = params_df$N_star[i],        # F̄: Female density at inflection
    Pbar = params_df$Pbar[i],          # P̄: Poult density at inflection (April)
    slope = params_df$r[i],            # Growth rate
    K_density = params_df$K_density[i] # Carrying capacity
  )
}

# Save as RDS for R
saveRDS(regionParams, "Data/rec_weather_prediction.rds")

# Create CSV version for MATLAB
region_params_df <- do.call(rbind, lapply(names(regionParams), function(region) {
  params <- regionParams[[region]]
  data.frame(
    WMU_Group = region,
    Fbar = params$Fbar,
    Pbar = params$Pbar,
    slope = params$slope,
    K_density = params$K_density
  )
})) %>%
  filter(!WMU_Group == "Group 10")  # Exclude Group 10

write.csv(region_params_df, "Data/rec_weather_prediction.csv", row.names = FALSE)

cat("✓ Parameter sets saved\n")

################################################################################
# 8. SUMMARY AND NEXT STEPS
################################################################################

cat("\n========== SUMMARY ==========\n")
cat("✓ April PPB predictions generated using weather model\n")
cat("✓ April Pbar calculated (PPB × HWB)\n")
cat("✓ Region-specific parameters created for MDP\n")
cat("✓ Comparison with September parameters complete\n")

cat("\nOutputs created:\n")
cat("  1. rec_weather_prediction.csv - April parameters for MATLAB\n")
cat("  2. params_ppb_summary_april2023.csv - Complete April parameter table\n")
cat("  3. april_vs_september_comparison.csv - Parameter comparison\n")

cat("\nNext steps:\n")
cat("  - Use April parameters in MATLAB decision model\n")
cat("  - Compare optimal decisions: January vs April vs September\n")
cat("  - Calculate value of information (VOI) for each timing\n")
cat("  - Evaluate utility gain from waiting for better information\n")

cat("\n========================================\n\n")

################################################################################
# END OF SCRIPT
################################################################################