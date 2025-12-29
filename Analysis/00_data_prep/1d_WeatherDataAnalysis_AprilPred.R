# Suppose you have a dataset with April temperatures and MU/year info
# e.g., April_weather_df has columns: MU, Year, scale_avg_temperature
rm(list=ls())
gc()
library(glmmTMB)
library(dplyr)

# Get April weather
## run 1a for 2024

# Format df
# run 1b for 2024
april23_weather_df <- readRDS("Data/Rec_data/april23_scaled_weather.rds")

#-----------------------------------------------------------------------------
# Load data
#-----------------------------------------------------------------------------
month <- "april"
scaled_weather <- readRDS(paste0("Data/Rec_data/", month, "_scaled_weather2.rds"))

#------------------------------------------------------------------------
# 1. Poults per brood (PPB) Analysis
#------------------------------------------------------------------------
ppb_df <- readRDS("Data/Rec_data/ph_df_aug31.rds")

# Join the datasets
ppb <- ppb_df %>%
  left_join(scaled_weather) 

# Create the model formula
ppb_formula <- PHratio ~ 0 + scale_avg_temperature + (1|MU)

# Fit the refined model
m.refined.ppb <- glmmTMB(ppb_formula, data = ppb, family = gaussian(link = "identity"))

# Summary of the model
summary(m.refined.ppb)

#------------------------------------------------------------------------
# 2. # Use the fitted model (m.refined.pph)
#------------------------------------------------------------------------
april_preds <- predict(
  m.refined.ppb,
  newdata = april23_weather_df,
  type = "response",
  se.fit = TRUE,
  re.form = NULL  # includes random effect for MU
)

# Attach results
april23_weather_df$PPB_pred <- as.numeric(april_preds$fit)
april23_weather_df$PPB_se <- april_preds$se.fit
april23_weather_df$PPB_lower <- april_preds$fit - 1.96 * april_preds$se.fit
april23_weather_df$PPB_upper <- april_preds$fit + 1.96 * april_preds$se.fit

# Get mean PPB per group in that year
mean_april_ppb <- april23_weather_df %>% 
  group_by(Group) %>% 
  filter(!Group == "Group 10") %>% 
  summarise(PPB = mean(PPB_pred, na.rm = TRUE))

# Format file to match 0b csv file formatting
# Load density data for scaling BBS indices to actual densities
# If this file doesn't exist, we'll create placeholder values
if (file.exists("Analysis/0a_findingDensity.R")) {
  message("Loading actual density estimates...")
  source("Analysis/0a_findingDensity.R")
  k_density <- abundance_fall[, 2:3] #grab only 2023
  hwb_ratio <- hwb_overall$total
  ppb_value <- mean_april_ppb$PPB # prediction from above 
} else {
  message("WARNING: Density estimate file not found, using placeholder values.")
  message("Replace these with your actual density estimates before using results!")
  
  # Placeholder values for HWB and PPB
  hwb_ratio <- hwb_ratio  # Around 0.79 based on your comment
  ppb_value <- ppb_value   # Around 4.1 based on your comment
}

message("Scaling BBS indices to density values...")

# Join parameters with density data
params <- readRDS("Data/LogisticGrowthModParams_2023.rds")
params_df <- as.data.frame(params) %>% 
  mutate(WMU_Group = rownames(params)) %>% 
  left_join(k_density, by = "WMU_Group") %>% 
  rename(K_density = total)

# Fill in missing values with the mean
params_df$K_density <- ifelse(is.na(params_df$K_density), 
                              mean(params_df$K_density, na.rm = TRUE), 
                              params_df$K_density)

# Calculate scaling factors and inflection points
params_df <- params_df %>%
  mutate(
    # Calculate scaling factor to convert index to density
    scaling_factor = K_density / K,
    
    # Calculate female density at inflection point (half carrying capacity)
    N_star = K/2 * scaling_factor,
    
    # Add region name for plotting
    Region = WMU_Group,
    
    # Add poults per brood from your data
    ppb = ppb_value,
    
    # Calculate poults per hen (Pbar) at inflection point
    Pbar = ppb * hwb_ratio
  )

# Save the complete parameters dataframe
write.csv(params_df, "Data/params_ppb_summary_april2023.csv", row.names = FALSE)

message("Generating density-dependent reproduction curves...")

# Create a sequence of density values for plotting the DD curve
density_seq <- seq(0, max(params_df$K_density, na.rm = TRUE) * 1.2, length.out = 100)

# Generate PPH curves for each region
pph_curves <- data.frame()

for(i in 1:nrow(params_df)) {
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

# Create a list of region-specific parameter sets for the decision model
message("Creating region-specific parameter sets for the decision model...")
regionParams <- list()
for(i in 1:nrow(params_df)) {
  region_name <- params_df$WMU_Group[i]
  regionParams[[region_name]] <- list(
    Fbar = params_df$N_star[i],        # Female density at inflection point
    Pbar = params_df$Pbar[i],          # Poults per hen at inflection point
    slope = params_df$r[i] / 4,        # Scaled slope parameter
    K_density = params_df$K_density[i] # Carrying capacity
  )
}

# Save the region-specific parameters for use in the MATLAB model
saveRDS(regionParams, "Data/rec_weather_prediction.rds")

# Create dataframe version for CSV export
region_params_df <- do.call(rbind, lapply(names(regionParams), function(region) {
  params <- regionParams[[region]]
  data.frame(
    WMU_Group = region,
    Fbar = params$Fbar,
    Pbar = params$Pbar,
    slope = params$slope,
    K_density = params$K_density
  )
}))

write.csv(region_params_df, "Data/rec_weather_prediction.csv", row.names = FALSE)

# Read in data
x <- read.csv('Data/rec_weather_prediction.csv') %>% 
  filter(!WMU_Group == "Group 10")
print(x)

y <- read.csv('Data/region_parameters_for_mdp.csv')
print(y)

# Calculate the difference in Pbar
diff_Pbar <- x$Pbar - y$Pbar
median(diff_Pbar)

1.8 - median(diff_Pbar)
max(x$Fbar)
