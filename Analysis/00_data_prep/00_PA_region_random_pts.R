################################################################################
# Generate Random Sampling Points for Weather Data Collection
################################################################################
#
# Purpose: Create a grid of 300 random points within each WMU region for
#          downloading weather data. These points provide spatially 
#          representative weather sampling across each region.
#
# Method: For each WMU group, sample points using a regular grid pattern
#         to ensure even spatial coverage
#
# Input:  PGC_BNDWildlifeManagementUnits2021.shp (PA WMU boundaries)
#
# Output: WMUgroup_pts.rds (coordinates of sampling points by WMU group)
#         Map visualization showing point distribution
#
# Author: Veronica A. Winter
# Date: January 2025
################################################################################

rm(list = ls())
gc()

# Set seed for reproducibility
set.seed(33456)

library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

################################################################################
# 1. LOAD AND PREPARE WMU BOUNDARIES
################################################################################

# Load Pennsylvania WMU shapefile
pa_wmu <- st_read("Data/PA_GIS/PGC_BNDWildlifeManagementUnits2021.shp")

# Group individual WMUs into biologically meaningful regions
# Same grouping as used in BBS analysis
groups <- pa_wmu %>%
  mutate(WMU_Group = case_when(
    substr(WMU_ID, 1, 1) == "1" ~ "Group 1",                # Northern tier
    WMU_ID %in% c("2A", "2C", "2D", "2E") ~ "Group 2",     # North-central
    WMU_ID == "2B" ~ "Group 3",                             # Central mountains
    WMU_ID %in% c("2F", "2G", "2H", "3A", "3B") ~ "Group 4", # South-central
    WMU_ID == "3C" ~ "Group 5",                             # Southeast
    WMU_ID == "3D" ~ "Group 6",                             # Southwest
    WMU_ID %in% c("4A", "4B", "4D") ~ "Group 7",           # West-central
    WMU_ID %in% c("4C", "4E") ~ "Group 8",                 # Northwest
    WMU_ID %in% c("5A", "5B") ~ "Group 9",                 # Allegheny north
    WMU_ID %in% c("5C", "5D") ~ "Group 10"                 # Allegheny south
  ))

# Dissolve individual WMUs into unified group polygons
pa_groups <- groups %>% 
  group_by(WMU_Group) %>% 
  summarize(geometry = st_union(geometry)) %>%
  ungroup()

################################################################################
# 2. GENERATE SAMPLING POINTS FOR EACH WMU GROUP
################################################################################

cat("\n========== GENERATING SAMPLING POINTS ==========\n")

# Initialize empty data frame for all points
group_pts <- data.frame()

# Loop through each WMU group
for (i in 1:10) {
  
  cat("Sampling Group", i, "\n")
  
  # Filter to current group
  grp <- pa_groups %>% 
    filter(WMU_Group == WMU_Group[i]) %>% 
    # Transform to planar coordinate system for distance-based sampling
    # PA State Plane (EPSG:2272) uses feet as units
    st_transform(crs = 2272)
  
  # Sample points using regular grid pattern
  # Regular sampling ensures even spatial coverage
  sample_points <- st_sample(
    grp, 
    size = 300,         # Number of points per region
    type = "regular"    # Regular grid (not random)
  ) %>%
    # Convert to sf object
    st_sf(geometry = ., crs = 2272)
  
  # Transform back to geographic coordinates (lat/lon)
  # NAD83 (EPSG:4269) for compatibility with Daymet
  sample_points <- st_transform(sample_points, crs = 4269)
  
  # Extract coordinates as data frame
  points_as_coordinates <- st_coordinates(sample_points) %>%
    as.data.frame()
  
  # Add group identifier
  points_as_coordinates$Group <- pa_groups$WMU_Group[i]
  
  # Append to master data frame
  group_pts <- bind_rows(group_pts, points_as_coordinates)
  
  cat("  Generated", nrow(points_as_coordinates), "points\n")
}

cat("\nTotal points generated:", nrow(group_pts), "\n")

################################################################################
# 3. SAVE SAMPLING POINTS
################################################################################

# Save as RDS for use in weather download script
saveRDS(group_pts, "../turkey_SDP/Data/WMUgroup_pts.rds")

cat("✓ Points saved to: ../turkey_SDP/Data/WMUgroup_pts.rds\n")

################################################################################
# 4. VISUALIZE SAMPLING POINT DISTRIBUTION
################################################################################

cat("\n========== CREATING VISUALIZATION ==========\n")

# Convert points back to sf object for plotting
group_pts_sf <- st_as_sf(group_pts, 
                         coords = c("X", "Y"), 
                         crs = 4269)

# Create map showing point distribution
sampling_map <- ggplot() +
  # Add WMU group boundaries
  geom_sf(data = pa_groups, 
          aes(fill = WMU_Group), 
          alpha = 0.3, 
          color = "black",
          linewidth = 0.5) +
  # Add sampling points
  geom_sf(data = group_pts_sf, 
          aes(color = Group), 
          size = 1.5,
          alpha = 0.6) +
  # Labels and formatting
  labs(
    title = "Weather Sampling Points by WMU Region",
    subtitle = paste0(nrow(group_pts), " total points (300 per region)"),
    fill = "WMU Group",
    color = "Sampling\nGroup"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

print(sampling_map)

# Save the map
ggsave("Dataviz/weather_sampling_points_map.png", 
       plot = sampling_map,
       width = 10, 
       height = 8, 
       dpi = 300,
       bg = "white")

cat("✓ Map saved to: Dataviz/weather_sampling_points_map.png\n")

################################################################################
# 5. PRINT SUMMARY STATISTICS
################################################################################

cat("\n========== SUMMARY STATISTICS ==========\n")

# Points per group
points_per_group <- group_pts %>%
  group_by(Group) %>%
  summarise(n_points = n()) %>%
  arrange(Group)

print(points_per_group)

# Coordinate ranges (verify reasonable lat/lon values)
cat("\nCoordinate ranges:\n")
cat("Longitude: ", round(min(group_pts$X), 2), "to", round(max(group_pts$X), 2), "\n")
cat("Latitude: ", round(min(group_pts$Y), 2), "to", round(max(group_pts$Y), 2), "\n")

# Expected ranges for Pennsylvania:
# Longitude: approximately -80.5 to -74.7
# Latitude: approximately 39.7 to 42.3

################################################################################
# 6. USAGE NOTES
################################################################################

cat("\n========== NEXT STEPS ==========\n")
cat("These points are used in:\n")
cat("  03_download_weather_data.R\n")
cat("\nEach point will download daily weather data from Daymet\n")
cat("Total API calls per year: ", nrow(group_pts), " (one per point)\n")
cat("Expected download time: ~", round(nrow(group_pts) * 0.5 / 60, 1), "minutes\n")
cat("===============================\n\n")

################################################################################
# NOTES ON SAMPLING DESIGN
################################################################################

# Why 300 points per region?
#   - Provides good spatial coverage without excessive API calls
#   - Regular grid ensures no clustering or gaps
#   - Sufficient to capture spatial weather variation within regions
#
# Why regular grid vs. random?
#   - Regular grid has more even spatial coverage
#   - Reduces risk of sampling bias
#   - Still provides representative weather sampling
#
# Coordinate system choices:
#   - PA State Plane (2272) for distance-based sampling
#   - NAD83 (4269) for geographic coordinates (matches Daymet)
#
# Alternative approach (if 300 points is too many):
#   - Reduce to 30-50 points per region
#   - Use stratified random sampling by elevation or habitat
#   - Weight points by land area if regions vary in size

################################################################################
# END OF SCRIPT
################################################################################