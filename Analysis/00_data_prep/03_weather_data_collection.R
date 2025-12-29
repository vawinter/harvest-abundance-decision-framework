################################################################################
# Download Weather Covariates for Recruitment Prediction
################################################################################
#
# Purpose: Download spring weather data (March-June) from Daymet for random
#          points within each WMU region. These covariates are used to predict
#          recruitment in the April decision timing model.
#
# Method: For each WMU region, sample 30 random points and download daily
#         weather data (temperature, precipitation) for spring months
#
# Input:  WMUgroup_pts.rds (random points by WMU group)
#         Generated from '00_PA_region_random_pts.R'
#
# Output: weather_covar_groups_combined.rds (spring weather by region and year)
#         Intermediate files saved as backup during processing
#
# Data Source: Daymet (https://daymet.ornl.gov/)
#              1 km resolution gridded climate data for North America
#
# Author: Veronica A. Winter
# Date: January 2025
################################################################################

rm(list = ls())
gc()

# Set seed for reproducibility
set.seed(33456)

library(daymetr)   # Download Daymet weather data
library(dplyr)
library(lubridate)
library(stringr)

################################################################################
# 1. LOAD SAMPLING POINTS
################################################################################

# Load random points for each WMU group
# Each group has 30 points sampled across its geographic extent
grp <- readRDS("Data/PA_GIS/WMUgroup_pts.rds")

# Years to download weather data for
start_year <- 2017
end_year <- 2023

################################################################################
# 2. DOWNLOAD WEATHER DATA FOR EACH REGION
################################################################################

# Initialize empty data frame for all weather data
group_weather <- data.frame()

# Get list of unique regions
regions <- unique(grp$Group)
total_regions <- length(regions)

cat("\n========== DOWNLOADING WEATHER DATA ==========\n")
cat("Total regions to process:", total_regions, "\n")
cat("Years:", start_year, "to", end_year, "\n")
cat("Months: March - June (spring recruitment period)\n")
cat("==============================================\n\n")

# Loop through each WMU region
for (i in 1:length(regions)) {
  
  cat("Processing", regions[i], "(", i, "of", total_regions, ")\n")
  
  # Filter points for current region
  dat <- grp %>% 
    filter(Group == regions[i])
  
  # Initialize data frame for this region
  group_temp_weather <- data.frame()
  
  # Loop through each sampling point within the region
  for (g in seq_len(nrow(dat))) {
    
    dat_row <- dat[g, ]
    
    # Wrap in error handling to prevent one failed point from stopping entire process
    tryCatch({
      
      # Extract coordinates
      lon <- dat_row$X
      lat <- dat_row$Y
      
      # Download Daymet data for this location
      # Returns daily values for: tmin, tmax, prcp, srad, vp, swe, dayl
      temp_daymet <- download_daymet(
        lat = lat,
        lon = lon,
        start = start_year,
        end = end_year,
        silent = TRUE  # Suppress progress messages
      )
      
      # Extract and transform weather data
      temp_daymet_df <- temp_daymet$data %>% 
        mutate(
          # Convert day-of-year to calendar date
          date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"),
          # Add separate date components
          day = day(date),
          month = month(date),
          year = year(date),
          # Add region identifier
          group = dat_row$Group
        ) %>% 
        # Filter to spring months only (March-June)
        # This is the critical period for poult survival
        filter(month %in% c(3, 4, 5, 6))
      
      # Append to region's weather data
      group_temp_weather <- bind_rows(group_temp_weather, temp_daymet_df)
      
    }, error = function(e) {
      # Log errors but continue processing
      cat("  ERROR: Point", g, "failed -", conditionMessage(e), "\n")
    })
  }
  
  # Add region data to main dataset
  group_weather <- bind_rows(group_weather, group_temp_weather)
  
  # Save intermediate results after each region (backup in case of crash)
  saveRDS(group_weather, 
          paste0("Data/Rec_data/weather_covar_groups_temp_", i, ".rds"))
  
  # Print progress
  cat("  Completed", regions[i], "-", nrow(group_temp_weather), "rows\n")
  cat("  Cumulative total:", nrow(group_weather), "rows\n\n")
}

# Save initial combined file
saveRDS(group_weather, "Data/Rec_data/2023_weather_covar_groups.rds")

cat("\n✓ Initial download complete!\n")

################################################################################
# 3. COMBINE AND VERIFY DATA
################################################################################

cat("\n========== COMBINING DATA FILES ==========\n")

# Find all intermediate files
weather_files <- list.files(
  path = "Data/Rec_data/", 
  pattern = "weather_covar_groups_temp_\\d+\\.rds$",
  full.names = TRUE
)

cat("Found", length(weather_files), "intermediate files\n")

# Initialize empty data frame
combined_weather <- data.frame()

# Read and combine all files
for (file in weather_files) {
  cat("Reading:", basename(file), "\n")
  temp_data <- readRDS(file)
  combined_weather <- bind_rows(combined_weather, temp_data)
}

# Remove any duplicate rows
combined_weather <- distinct(combined_weather)

################################################################################
# 4. DATA QUALITY CHECKS
################################################################################

cat("\n========== DATA QUALITY CHECKS ==========\n")

# Check number of days per group per year
# Should have ~120-122 days (March-June = 31+30+31+30 = 122 days per point)
day_counts <- combined_weather %>%
  group_by(group, year) %>%
  summarise(
    n_days = n_distinct(date),
    n_points = n(),
    .groups = 'drop'
  ) %>%
  arrange(group, year)

print(day_counts)

# Check for missing groups
expected_groups <- paste("Group", 1:10)
actual_groups <- unique(combined_weather$group)
missing_groups <- setdiff(expected_groups, actual_groups)

if (length(missing_groups) > 0) {
  warning("Missing data for groups: ", paste(missing_groups, collapse = ", "))
} else {
  cat("✓ All groups present\n")
}

# Summary statistics
cat("\n========== SUMMARY ==========\n")
cat("Total rows:", nrow(combined_weather), "\n")
cat("Groups:", paste(unique(combined_weather$group), collapse = ", "), "\n")
cat("Date range:", min(combined_weather$date), "to", max(combined_weather$date), "\n")
cat("Variables:", paste(names(combined_weather), collapse = ", "), "\n")

# Check for missing values in key variables
missing_summary <- combined_weather %>%
  summarise(
    missing_tmin = sum(is.na(tmin..deg.c.)),
    missing_tmax = sum(is.na(tmax..deg.c.)),
    missing_prcp = sum(is.na(prcp..mm.day.))
  )

if (sum(missing_summary) > 0) {
  warning("Missing values detected:")
  print(missing_summary)
} else {
  cat("✓ No missing values in key variables\n")
}

################################################################################
# 5. SAVE FINAL COMBINED FILE
################################################################################

# Save the combined file
output_file <- "Data/Rec_data/weather_covar_groups_combined.rds"
saveRDS(combined_weather, output_file)

cat("\n✓ Combined file saved as:", output_file, "\n")

################################################################################
# 6. CLEANUP INTERMEDIATE FILES (OPTIONAL)
################################################################################

cat("\n========== CLEANUP ==========\n")

response <- readline("Delete intermediate files? (yes/no): ")

if (tolower(trimws(response)) == "yes") {
  file.remove(weather_files)
  cat("✓ Intermediate files deleted\n")
} else {
  cat("Intermediate files preserved in Data/Rec_data/\n")
}

cat("\n========== PROCESS COMPLETE ==========\n")

################################################################################
# NOTES FOR NEXT STEPS
################################################################################

# This weather data is used in:
#   - April decision timing model (recruitment prediction)
#   - Weather-based updates to Pbar parameter
#
# Key variables for recruitment prediction:
#   - tmin..deg.c.: Minimum temperature (cold stress on poults)
#   - tmax..deg.c.: Maximum temperature (heat stress)
#   - prcp..mm.day.: Precipitation (wet conditions reduce survival)
#
# Typical workflow:
#   1. Run this script to download weather (do annually)
#   2. Calculate summary statistics (e.g., mean spring temp, total precip)
#   3. Fit recruitment model: poults ~ weather + density
#   4. Use predictions to update Pbar in April decision model

################################################################################
# END OF SCRIPT
################################################################################