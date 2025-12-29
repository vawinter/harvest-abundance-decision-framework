################################################################################
# BBS Logistic Growth Model and Density Scaling
################################################################################
# 
# Purpose: Fit logistic growth models to BBS indices and scale to absolute
#          densities to derive inflection point parameters (F̄, P̄, slope)
#          for each WMU region
#
# Method: 
#   1. Use bbsBayes2 with custom PA WMU stratification
#   2. Fit first-difference model to BBS data (1966-2023)
#   3. Fit logistic growth model with random effects by WMU group
#   4. Scale BBS indices to density using IPM estimates
#   5. Calculate inflection point parameters
#
# Input:  PGC_BNDWildlifeManagementUnits2024.shp (PA WMU boundaries)
#         O_24_abundance_summary.rds (from 01_extract_density_from_IPM.R)
#         BBS data (downloaded via bbsBayes2)
#
# Output: region_parameters_for_mdp.csv (F̄, P̄, slope, K by region)
#         Diagnostic plots and model summaries
#
# Author: Veronica A. Winter
# Date: May 2025
# Updated: August 2025 (density scaling), December 2025 (cleaned for MS)
################################################################################

rm(list = ls())
gc()

library(bbsBayes2)
library(dplyr)
library(ggplot2)
library(sf)
library(cmdstanr)
library(nlme)

################################################################################
# 1. LOAD HELPER FUNCTIONS
################################################################################

# Load density-dependent fecundity function
# This calculates poults per hen (PPH) as a function of female density
source("xx_pph_density_fun.R")

################################################################################
# 2. CREATE CUSTOM WMU STRATIFICATION
################################################################################

# Load Pennsylvania WMU boundaries
pa_wmu <- st_read("Data/PA_GIS/PGC_BNDWildlifeManagementUnits2024.shp")

# Group WMUs into biologically meaningful regions
# Grouping based on geographic proximity and habitat characteristics
groups <- pa_wmu %>%
  mutate(strata_name = case_when(
    substr(WMU_ID, 1, 1) == "1" ~ "Group 1",        # Northern tier
    WMU_ID %in% c("2A", "2C", "2D", "2E") ~ "Group 2",  # North-central
    WMU_ID == "2B" ~ "Group 3",                      # Central mountains
    WMU_ID %in% c("2F", "2G", "2H", "3A", "3B") ~ "Group 4",  # South-central
    WMU_ID == "3C" ~ "Group 5",                      # Southeast
    WMU_ID == "3D" ~ "Group 6",                      # Southwest
    WMU_ID %in% c("4A", "4B", "4D") ~ "Group 7",     # West-central
    WMU_ID %in% c("4C", "4E") ~ "Group 8",           # Northwest
    WMU_ID %in% c("5A", "5B") ~ "Group 9",           # Allegheny Plateau north
    WMU_ID %in% c("5C", "5D") ~ "Group 10"           # Allegheny Plateau south
  ))

# Dissolve individual WMUs into groups
pa_groups <- groups %>% 
  group_by(strata_name) %>% 
  summarize(geometry = st_union(geometry)) %>%
  # Transform to WGS84 for compatibility with BBS
  st_transform(crs = 4326) %>% 
  ungroup()

# Ensure valid geometries
pa_groups <- st_make_valid(pa_groups)

if (any(!st_is_valid(pa_groups))) {
  warning("Some geometries are invalid after processing")
}

# Create map of WMU groups
groups_map <- ggplot() +
  geom_sf(data = pa_groups, aes(fill = strata_name), 
          color = "black", alpha = 0.6) + 
  geom_sf_text(data = pa_groups, aes(label = strata_name), size = 4) +
  labs(x = "", y = "", fill = "WMU Group",
       title = "Pennsylvania WMU Groups for BBS Analysis") +
  theme_classic() +
  theme(legend.position = "right")

ggsave("Dataviz/group_map.png", plot = groups_map,
       width = 8, height = 5, dpi = 700, bg = "white")

################################################################################
# 3. DOWNLOAD AND PREPARE BBS DATA
################################################################################

# Download BBS data (only needed once)
# Uncomment if first time running:
# fetch_bbs_data()

# Create stratification using custom WMU boundaries for Wild Turkey
strat <- stratify(by = "custom", 
                  strata_custom = pa_groups, 
                  species = "Wild Turkey")

# Prepare data for analysis
wild_turkey <- prepare_data(strat, 
                            min_year = 1966,
                            max_year = 2023)

################################################################################
# 4. FIT FIRST-DIFFERENCE MODEL TO BBS DATA
################################################################################

# The first-difference model treats abundance as a random walk
# This allows flexible trends without assuming a specific curve shape
# Model: θ_t = θ_{t-1} + ε_t, where ε_t ~ Normal(0, σ²)

wt_mod <- prepare_model(wild_turkey,
                        model = "first_diff")

# Run the model
# This can take 10-30 minutes depending on your system
message("Running BBS model... (this may take a while)")
model_run <- run_model(wt_mod, 
                       iter_sampling = 2000, 
                       iter_warmup = 2000, 
                       parallel_chains = 2,
                       output_basename = "WildTurkey_Regions",
                       output_dir = "Results/")

# Check convergence (Rhat should be < 1.01)
converge_n_smooth <- get_convergence(model_run) %>%
  arrange(-rhat)

message("Convergence diagnostics:")
print(converge_n_smooth)

if (any(converge_n_smooth$rhat > 1.01)) {
  warning("Some parameters have Rhat > 1.01. Consider increasing iterations.")
}

################################################################################
# 5. GENERATE INDICES AND TRENDS
################################################################################

# Generate annual population indices for each region
indices <- generate_indices(model_run)

# Plot indices with observed means
p <- plot_indices(indices = indices,
                  add_observed_means = TRUE)

# Display plots for each region
print(p)

# Generate trends across all years
trends_all <- generate_trends(indices)
final_trends <- bind_rows(trends_all)

# Extract index values for modeling
final_indices <- indices$indices

# Save indices for future use
saveRDS(final_indices, "Data/BBS-WMU_group.rds")

################################################################################
# 6. FIT LOGISTIC GROWTH MODEL
################################################################################

# Prepare data for logistic growth modeling
pa_data <- final_indices %>% 
  filter(!region %in% "continent") %>% 
  rename(WMU_Group = region) %>% 
  # Remove groups with insufficient data
  group_by(WMU_Group) %>%
  filter(n() >= 10) %>%  # Need at least 10 years
  ungroup() %>%
  # Create time variables
  mutate(
    Years = year - min(year),
    Years_Centered = Years - mean(Years)
  ) %>% 
  select(index, WMU_Group, Years, Years_Centered) %>%
  rename(Index = index)

message(paste("Modeling", length(unique(pa_data$WMU_Group)), "WMU groups"))

# Logistic growth model: N(t) = K / (1 + exp(-r * t))
# K = carrying capacity, r = intrinsic growth rate

# Step 1: Fit fixed-effects model for starting values
logistic_nls <- nls(
  Index ~ K / (1 + exp(-r * Years_Centered)),
  data = pa_data,
  start = list(K = max(pa_data$Index), r = 0.1)
)

message("Fixed-effects starting values:")
print(coef(logistic_nls))

# Step 2: Fit random-effects model (K and r vary by WMU group)
logistic_mixed_model <- nlme(
  Index ~ K / (1 + exp(-r * Years_Centered)),
  data = pa_data,
  fixed = K + r ~ 1,
  random = K + r ~ 1 | WMU_Group,
  start = list(fixed = coef(logistic_nls)),
  control = nlmeControl(msMaxEval = 2000, pnlsMaxIter = 100)
)

# Extract group-specific parameters
params <- round(coef(logistic_mixed_model), 4)
carrying_capacity <- params$K
growth_rate <- params$r

# Add predictions to data for plotting
pa_data$predicted <- predict(logistic_mixed_model, newdata = pa_data)

# Calculate R-squared
residuals <- pa_data$Index - pa_data$predicted
SS_tot <- sum((pa_data$Index - mean(pa_data$Index, na.rm = TRUE))^2, na.rm = TRUE)
SS_res <- sum(residuals^2, na.rm = TRUE)
R_squared <- 1 - (SS_res / SS_tot)

message(paste("Logistic model R-squared:", round(R_squared, 3)))

# Save model parameters
saveRDS(params, "Data/LogisticGrowthModParams_2023.rds")

################################################################################
# 7. SCALE BBS INDICES TO ABSOLUTE DENSITIES
################################################################################

# Load density data from IPM (from 01_extract_density_from_IPM.R)
# This provides the scaling factor to convert BBS indices to real densities

if (file.exists("Analysis/01_extract_density_from_IPM.R")) {
  message("Loading IPM density estimates...")
  source("Analysis/01_extract_density_from_IPM.R")
  
  k_density <- abundance_fall %>% 
    select(WMU_Group, total)
  
  hwb_ratio <- hwb_overall$total
  ppb_value <- ppb_overall$total
  
} else {
  stop("ERROR: Cannot find density estimation script. Run 01_extract_density_from_IPM.R first.")
}

message("Scaling BBS indices to absolute densities...")

# Join logistic parameters with density data
params_df <- as.data.frame(params) %>% 
  mutate(WMU_Group = rownames(params)) %>% 
  left_join(k_density, by = "WMU_Group") %>% 
  rename(K_density = total)

# Fill missing values with mean (conservative approach)
params_df$K_density <- ifelse(is.na(params_df$K_density), 
                              mean(params_df$K_density, na.rm = TRUE), 
                              params_df$K_density)

# Calculate scaled parameters
params_df <- params_df %>%
  mutate(
    # Scaling factor: converts BBS index units to density units
    scaling_factor = K_density / K,
    
    # Female density at inflection point (F̄)
    # Inflection occurs at K/2 in logistic growth
    N_star = (K / 2) * scaling_factor,
    
    # Add poults per brood
    ppb = ppb_value,
    
    # Poult density at inflection (P̄)
    # P̄ = ppb × hen-with-brood ratio
    Pbar = ppb * hwb_ratio,
    
    # Region name for clarity
    Region = WMU_Group
  )

# Save complete parameters
write.csv(params_df, "Data/params_ppb_summary_2023.csv", row.names = FALSE)

message("Parameter calculation complete!")

################################################################################
# 8. GENERATE DENSITY-DEPENDENT REPRODUCTION CURVES
################################################################################

message("Generating density-dependent reproduction curves...")

# Create sequence of female densities for plotting
density_seq <- seq(0, max(params_df$K_density, na.rm = TRUE) * 1.2, 
                   length.out = 100)

# Calculate PPH curves for each region
# Using the reverse sigmoidal fecundity function
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
# 9. CREATE REGION-SPECIFIC PARAMETER SETS FOR MDP MODEL
################################################################################

message("Creating region-specific parameter sets for decision model...")

regionParams <- list()

for (i in 1:nrow(params_df)) {
  region_name <- params_df$WMU_Group[i]
  regionParams[[region_name]] <- list(
    Fbar = params_df$N_star[i],      # F̄: Female density at inflection
    Pbar = params_df$Pbar[i],        # P̄: Poults per hen at inflection
    slope = params_df$r[i],          # Intrinsic growth rate
    K_density = params_df$K_density[i]  # Carrying capacity
  )
}

# Exclude Group 10 if insufficient data
regionParams <- regionParams[!names(regionParams) %in% "Group 10"]

# Save as RDS for R
saveRDS(regionParams, "Data/region_parameters_for_mdp.rds")

# Create CSV version for MATLAB
region_params_df <- do.call(rbind, lapply(names(regionParams), function(region) {
  params <- regionParams[[region]]
  data.frame(
    WMU_Group = region,
    Fbar = params$Fbar,
    Pbar = params$Pbar,
    slope = round(params$slope, 2),
    K_density = params$K_density
  )
})) %>%
  filter(!WMU_Group == "Region 10")

write.csv(region_params_df, "Data/region_parameters_for_mdp.csv", row.names = FALSE)

# Add metadata
attr(regionParams, "generated_on") <- Sys.time()
attr(regionParams, "data_source") <- "BBS analysis with IPM density scaling"

################################################################################
# 10. CREATE VISUALIZATIONS
################################################################################

message("Creating visualizations...")

# Scale BBS data to density units for plotting
pa_data_scaled <- pa_data %>%
  left_join(params_df %>% select(WMU_Group, scaling_factor), 
            by = "WMU_Group") %>%
  mutate(Density = Index * scaling_factor)

# Scale predictions
predictions <- predict(logistic_mixed_model, newdata = pa_data)
pa_data_scaled$predictions_scaled <- predictions * pa_data_scaled$scaling_factor

# Add calendar years
min_year <- min(final_indices$year)
pa_data_scaled$Year <- pa_data_scaled$Years + min_year

#------------------------------------------------------------------------------
# Plot 1: Population density over time
#------------------------------------------------------------------------------
DD_scaled <- ggplot(pa_data_scaled, aes(x = Year, y = Density, color = WMU_Group)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predictions_scaled), linewidth = 1) +
  geom_hline(data = params_df, aes(yintercept = K_density, color = WMU_Group),
             linetype = "dashed", alpha = 0.5) +
  labs(title = "Wild Turkey Population Density Over Time",
       y = expression("Female Density (birds/km"^2*")"), 
       x = "Year",
       color = "WMU Group") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

#------------------------------------------------------------------------------
# Plot 2: Density-dependent reproduction curves
#------------------------------------------------------------------------------
pph_plot <- ggplot(pph_curves, aes(x = female_density, y = pph, color = WMU_Group)) +
  geom_line(linewidth = 1) +
  # Mark inflection points
  geom_point(data = params_df, aes(x = N_star, y = Pbar), size = 3) +
  labs(title = "Density-Dependent Reproduction",
       y = "Poults Per Hen", 
       x = expression("Female Density (birds/km"^2*")"),
       color = "WMU Group") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

#------------------------------------------------------------------------------
# Plot 3: Standardized curves for comparison
#------------------------------------------------------------------------------
standardized_curves <- pph_curves %>%
  group_by(WMU_Group) %>%
  mutate(
    # Standardize density relative to carrying capacity
    std_density = female_density / 
      params_df$K_density[params_df$WMU_Group == WMU_Group[1]],
    # Standardize PPH relative to maximum
    std_pph = pph / max(pph, na.rm = TRUE)
  ) %>%
  ungroup()

combined_plot <- ggplot() +
  geom_line(data = standardized_curves, 
            aes(x = std_density, y = std_pph, color = WMU_Group),
            linewidth = 1) +
  # Reference lines at inflection (0.5, 0.5)
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray50") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50") +
  labs(title = "Standardized Density-Dependent Relationships",
       x = "Relative Female Density (proportion of K)",
       y = "Relative Reproductive Output",
       color = "WMU Group") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

# Save plots
ggsave("Dataviz/turkey_density_over_time.png", plot = DD_scaled,
       width = 10, height = 7, dpi = 300, bg = "white")
ggsave("Dataviz/turkey_reproduction_curves.png", plot = pph_plot,
       width = 10, height = 7, dpi = 300, bg = "white")
ggsave("Dataviz/turkey_standardized_curves.png", plot = combined_plot,
       width = 10, height = 7, dpi = 300, bg = "white")

################################################################################
# 11. PRINT SUMMARY STATISTICS
################################################################################

params_summary <- params_df %>%
  select(WMU_Group, K, r, K_density, N_star, Pbar) %>%
  arrange(WMU_Group)

message("\n========== KEY PARAMETERS FOR DECISION MODEL ==========")
print(params_summary)
message("=======================================================\n")

# Create text summary
param_text <- paste0(
  "Parameters for Density-Dependent Wild Turkey Model\n\n",
  "Overall Mean Values:\n",
  "- Carrying Capacity (K): ", round(mean(params_df$K_density, na.rm = TRUE), 2), " birds/km²\n",
  "- Inflection Point (F̄): ", round(mean(params_df$N_star, na.rm = TRUE), 2), " birds/km²\n",
  "- Poults per Hen at Inflection (P̄): ", round(mean(params_df$Pbar, na.rm = TRUE), 2), "\n",
  "- Growth Rate (r): ", round(mean(params_df$r, na.rm = TRUE), 2), "\n\n",
  "Density-Dependent Fecundity Function:\n",
  "P(F) = [2η × P̄] / [η + 1 + (η-1) × (F/F̄)^η]\n",
  "where η = 1 - 2(F̄/P̄) × slope\n\n",
  "Data Sources:\n",
  "- BBS data: 1966-2023\n",
  "- IPM density estimates: Winter et al. (in review)\n",
  "- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n"
)

writeLines(param_text, "Dataviz/model_parameters_summary.txt")

message("\n✓ Analysis complete!")
message("✓ Parameters saved to: Data/region_parameters_for_mdp.csv")
message("✓ Plots saved to: Dataviz/")
message("✓ Summary text saved to: Data/model_parameters_summary.txt")

################################################################################
# END OF SCRIPT
################################################################################