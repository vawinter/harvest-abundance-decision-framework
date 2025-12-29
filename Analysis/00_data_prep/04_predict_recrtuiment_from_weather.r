################################################################################
# Predict Recruitment from Weather Covariates
################################################################################
#
# Purpose: 
#   1. Fit models for poults-per-brood (PPB) and hen-with-brood (HWB) from survey data
#   2. Process spring weather data with temporal weighting around nest initiation
#   3. Create weather-weighted predictors for April recruitment forecasting
#
# Method:
#   - GLMMs for PPB and HWB with year, day-of-year, and WMU random effects
#   - Weather data weighted by proximity to median nest initiation date (DOY 100)
#   - Normalized weighted predictors for temperature, precipitation, and snow
#
# Input:  
#   PA.rds (poult-hen survey data)
#   2023_weather_covar_groups_combined.rds (from 03_download_weather_data.R)
#   20250131_NestAttempts_allbirds.csv (nest timing data)
#
# Output: 
#   ph_df_aug31.rds (predicted poults per brood)
#   hwb_df_aug31.rds (predicted hen-with-brood ratio)
#   pph_df_aug31.rds (predicted poults per hen = PPB × HWB)
#   april23_scaled_weather.rds (weather covariates for April model)
#
# Author: Veronica A. Winter
# Date: February 2025
################################################################################

rm(list = ls())
gc()

set.seed(33456)

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(glmmTMB)
library(car)

################################################################################
# PART 1: MODEL POULTS PER BROOD (PPB)
################################################################################

cat("\n========== MODELING POULTS PER BROOD ==========\n")

# Load survey data
df.pa <- readRDS("Data/IPM/PA_sss.rds")

# Filter low-quality observations
df1 <- df.pa %>%
  # Remove sightings with >25% unknown sex/age
  filter(Unknown / Total < 0.25) %>%
  # Remove unrealistic PPB values (>16 poults/hen)
  filter(PHratio <= 16) %>%
  # Remove likely misidentified groups (8+ hens with 0 poults)
  filter(!(Hens >= 8 & Poults == 0))

# Create scaled day-of-year and quadratic term
df1$doy.scale <- scale(df1$doy)
df1$doy.2 <- df1$doy.scale^2

# Filter to analysis years and positive PPB
df.all <- df1 %>%
  filter(Year >= 2017, PHratio > 0) %>%
  mutate(
    Year = as.factor(Year),
    MU = as.factor(MU)
  )

cat("Data summary:\n")
cat("  Years:", paste(unique(df.all$Year), collapse = ", "), "\n")
cat("  WMU groups:", paste(unique(df.all$MU), collapse = ", "), "\n")
cat("  Total observations:", nrow(df.all), "\n\n")

# Fit GLMM for poults per brood
# Gamma distribution appropriate for continuous positive data
# Random intercept by WMU accounts for regional variation
m.ppb.all <- glmmTMB(
  PHratio ~ 0 + Year + doy.scale + doy.2 + (1|MU),  
  data = df.all, 
  family = Gamma(link = log)
)

cat("PPB model summary:\n")
print(summary(m.ppb.all))

################################################################################
# PART 2: PREDICT PPB AT AUGUST 31 (END OF BROOD-REARING SEASON)
################################################################################

# Convert August 31 (DOY 244) to scaled units
aug31 <- (244 - attr(df1$doy.scale, "scaled:center")) / 
  attr(df1$doy.scale, "scaled:scale")
aug31.2 <- aug31^2

# Create prediction grid for all year × WMU combinations
pred.pop <- expand.grid(
  Year = unique(df.all$Year), 
  doy.scale = aug31, 
  doy.2 = aug31.2, 
  MU = unique(df.all$MU)
)

# Predict PPB with confidence intervals
final4.all <- cbind(
  pred.pop,
  predict(m.ppb.all, newdata = pred.pop, se.fit = TRUE, type = "response")
) %>%
  mutate(
    lcl95 = fit - 1.96 * se.fit,
    ucl95 = fit + 1.96 * se.fit,
    class = factor(MU, levels = 1:10)
  )

# Visualize PPB predictions
ppb_plot <- ggplot(final4.all) +
  geom_errorbar(aes(x = Year, ymin = lcl95, ymax = ucl95, 
                    group = MU, color = class), 
                width = 0.2, position = position_dodge(0.2)) +
  geom_line(aes(x = Year, y = fit, group = MU, color = class), 
            position = position_dodge(0.2)) +
  geom_point(aes(x = Year, y = fit, group = MU, color = class), 
             position = position_dodge(0.2)) +
  labs(y = "Poults per Brood", x = "Year", color = "WMU Group",
       title = "Predicted Poults per Brood (August 31)") +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

print(ppb_plot)

# Format for saving
ph_df <- final4.all %>%
  mutate(
    Year = as.factor(Year),
    MU = as.factor(MU)
  ) %>%
  rename(PHratio = fit) %>%
  select(-c(class, se.fit, lcl95, ucl95))

saveRDS(ph_df, "Data/Rec_data/ph_df_aug31.rds")
cat("\n✓ PPB predictions saved\n")

################################################################################
# PART 3: MODEL HEN-WITH-BROOD RATIO (HWB)
################################################################################

cat("\n========== MODELING HEN-WITH-BROOD RATIO ==========\n")

# Prepare data for HWB analysis
df.all <- df1 %>%
  filter(Year > 2016) %>%
  mutate(
    Year = as.factor(Year),
    MU = as.factor(MU)
  )

# Expand data: each hen becomes a separate observation (0/1 for HWB)
df1.all <- uncount(df.all, Hens)

# Fit binomial GLMM for probability of hen having brood
m.hwb.all <- glmmTMB(
  HWB ~ 1 + Year + doy.scale + doy.2 + (1|MU), 
  dispformula = ~0,  # No overdispersion
  data = df1.all, 
  family = binomial()
)

cat("HWB model summary:\n")
print(summary(m.hwb.all))

################################################################################
# PART 4: PREDICT HWB AT AUGUST 31
################################################################################

# Same prediction grid as PPB
pred.pop <- expand.grid(
  Year = unique(df1.all$Year), 
  doy.scale = aug31, 
  doy.2 = aug31.2, 
  MU = unique(df1.all$MU)
)

# Predict HWB with confidence intervals
hens <- cbind(
  pred.pop,
  predict(m.hwb.all, newdata = pred.pop, se.fit = TRUE, type = "response")
) %>%
  mutate(
    lcl95 = fit - 1.96 * se.fit,
    ucl95 = fit + 1.96 * se.fit,
    class = MU
  )

# Visualize HWB predictions
hwb_plot <- ggplot(hens) +
  geom_errorbar(aes(x = Year, ymin = lcl95, ymax = ucl95, 
                    group = MU, color = class), 
                width = 0.2, position = position_dodge(0.2)) +
  geom_line(aes(x = Year, y = fit, group = MU, color = class), 
            position = position_dodge(0.2)) +
  geom_point(aes(x = Year, y = fit, group = MU, color = class), 
             position = position_dodge(0.2)) +
  labs(y = "Proportion with Poults", x = "Year", color = "WMU Group",
       title = "Predicted Hen-With-Brood Ratio (August 31)") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

print(hwb_plot)

# Format for saving
ph_df2 <- hens %>%
  mutate(
    Year = as.factor(Year),
    MU = as.factor(MU)
  ) %>%
  rename(HWB = fit) %>%
  select(-c(class, se.fit, lcl95, ucl95))

saveRDS(ph_df2, "Data/Rec_data/hwb_df_aug31.rds")
cat("\n✓ HWB predictions saved\n")

################################################################################
# PART 5: CALCULATE POULTS PER HEN (PPH = PPB × HWB)
################################################################################

cat("\n========== CALCULATING POULTS PER HEN ==========\n")

# Combine PPB and HWB to get overall recruitment metric
ph <- ph_df2 %>%
  left_join(ph_df, by = c("Year", "MU", "doy.scale", "doy.2")) %>%
  mutate(pph = PHratio * HWB)

# Summary statistics
pph_summary <- ph %>%
  group_by(Year) %>%
  summarise(
    mean_pph = mean(pph, na.rm = TRUE),
    min_pph = min(pph, na.rm = TRUE),
    max_pph = max(pph, na.rm = TRUE)
  )

cat("PPH summary by year:\n")
print(pph_summary)

saveRDS(ph, "Data/Rec_data/pph_df_aug31.rds")
cat("\n✓ PPH predictions saved\n")

################################################################################
# PART 6: DETERMINE MEDIAN NEST INITIATION DATE
################################################################################

cat("\n========== ANALYZING NEST TIMING ==========\n")

# Load nest attempt data
nests <- read.csv("Data/20250131_NestAttempts_allbirds.csv", header = TRUE)

# Format dates and extract day-of-year
nests <- nests %>%
  mutate(
    startI = as.Date(startI),
    year = year(startI),
    moday = yday(startI)
  )

# Calculate median nest initiation by year
nest_med_date <- nests %>% 
  group_by(year) %>% 
  summarise(median_doy = median(moday, na.rm = TRUE))

cat("Median nest initiation dates:\n")
print(nest_med_date)
cat("\nOverall median: ~DOY 138 (May 19)\n")

# For weather weighting, we'll center on DOY 100 (April 10)
# This is earlier than nest initiation to capture pre-laying conditions
median_doy <- 100

################################################################################
# PART 7: PROCESS AND WEIGHT WEATHER DATA
################################################################################

cat("\n========== PROCESSING WEATHER DATA ==========\n")

# Load weather data from Daymet download
combined_weather <- readRDS("Data/Rec_data/2023_weather_covar_groups_combined.rds")

# Calculate daily weather averages across sampling points within each region
daily_weather <- combined_weather %>%
  rename(Group = group) %>%
  group_by(Group, date) %>%
  summarise(
    # Average min and max temperature
    avg_temperature = mean((tmax..deg.c. + tmin..deg.c.) / 2, na.rm = TRUE),
    # Precipitation in mm
    avg_precipitation = mean(prcp..mm.day., na.rm = TRUE),
    # Snow water equivalent
    avg_snow = mean(swe..kg.m.2., na.rm = TRUE),
    .groups = "drop"
  )

cat("Weather data summary:\n")
cat("  Date range:", min(daily_weather$date), "to", max(daily_weather$date), "\n")
cat("  Groups:", paste(unique(daily_weather$Group), collapse = ", "), "\n\n")

# Define temporal window for weighting (±14 days around median DOY)
adj_doy <- 14

# Process weather data with temporal weighting
scaled_weather <- daily_weather %>%
  mutate(
    month = as.factor(month(date)),
    Year = as.factor(year(date)),
    doy = yday(date),
    MU = as.factor(str_extract(Group, "\\d+")),
    
    # Standardize temperature (mean = 0, sd = 1)
    scale_avg_temperature = as.numeric(scale(avg_temperature, 
                                             center = TRUE, scale = TRUE)),
    
    # Log-transform and scale precipitation and snow
    # Log transformation handles right skew in precip/snow data
    log_avg_snow = log(avg_snow + 1),
    scale_log_avg_snow = as.numeric(scale(log(avg_snow + 1), 
                                          center = FALSE, scale = TRUE)),
    log_avg_precip = log(avg_precipitation + 1),
    scale_log_avg_precip = as.numeric(scale(log(avg_precipitation + 1), 
                                            center = FALSE, scale = TRUE))
  ) %>%
  # Filter to critical period: ±14 days around DOY 100
  filter(doy >= (median_doy - adj_doy), doy <= (median_doy + adj_doy)) %>%
  # Apply temporal weighting
  group_by(Group, Year) %>%
  mutate(
    # Weight function: peak at median_doy, declining to 0.2 at edges
    weight = case_when(
      doy < (median_doy - adj_doy) ~ 0.2,
      doy >= (median_doy - adj_doy) & doy < median_doy ~ 
        0.2 + ((doy - (median_doy - adj_doy)) / adj_doy) * (1 - 0.2),
      doy == median_doy ~ 1.0,
      doy > median_doy & doy <= (median_doy + adj_doy) ~ 
        1.0 - ((doy - median_doy) / adj_doy) * (1 - 0.2),
      doy > (median_doy + adj_doy) ~ 0.2
    ),
    
    # Normalize weights to sum to 1 within each group-year
    normalized_weight = weight / sum(weight, na.rm = TRUE),
    
    # Apply weights to scaled weather variables
    weighted_temp = scale_avg_temperature * normalized_weight,
    weighted_snow = scale_log_avg_snow * normalized_weight,
    weighted_precip = scale_log_avg_precip * normalized_weight
  ) %>%
  # Normalize weighted predictors to sum to 1
  mutate(
    norm_weighted_temp = weighted_temp / sum(weighted_temp, na.rm = TRUE),
    norm_weighted_snow = weighted_snow / sum(weighted_snow, na.rm = TRUE),
    norm_weighted_precip = weighted_precip / sum(weighted_precip, na.rm = TRUE)
  ) %>%
  ungroup()

################################################################################
# PART 8: VERIFY WEATHER WEIGHTING
################################################################################

cat("\n========== VERIFYING WEATHER WEIGHTING ==========\n")

# Check that normalized weights sum to 1
weight_check <- scaled_weather %>%
  group_by(Group, Year) %>%
  summarise(
    sum_norm_temp = sum(norm_weighted_temp, na.rm = TRUE),
    sum_norm_snow = sum(norm_weighted_snow, na.rm = TRUE),
    sum_norm_precip = sum(norm_weighted_precip, na.rm = TRUE),
    .groups = "drop"
  )

cat("Normalized weight sums (should all be ~1.0):\n")
print(summary(weight_check))

# Visualize weighting scheme for one example
weight_plot <- scaled_weather %>%
  filter(Year == "2023", Group == "Group 1") %>%
  ggplot(aes(x = doy, y = weight)) + 
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = median_doy, linetype = "dashed", 
             color = "red", linewidth = 1) +
  annotate("text", x = median_doy + 2, y = 0.9, 
           label = paste("Median DOY =", median_doy),
           hjust = 0, color = "red") +
  labs(
    title = "Temporal Weighting Scheme for Weather Covariates",
    subtitle = "Example: Group 1, 2023",
    x = "Day of Year",
    y = "Weight"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(weight_plot)

# Verify temperature scaling
cat("\nTemperature scaling check:\n")
cat("  Mean:", round(mean(scaled_weather$scale_avg_temperature, na.rm = TRUE), 4), 
    "(should be ~0)\n")
cat("  SD:", round(sd(scaled_weather$scale_avg_temperature, na.rm = TRUE), 4), 
    "(should be ~1)\n")

################################################################################
# PART 9: SAVE PROCESSED WEATHER DATA
################################################################################

saveRDS(scaled_weather, "Data/Rec_data/april23_scaled_weather.rds")

cat("\n✓ Scaled weather data saved\n")

################################################################################
# SUMMARY AND NEXT STEPS
################################################################################

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Outputs created:\n")
cat("  1. ph_df_aug31.rds - Poults per brood predictions\n")
cat("  2. hwb_df_aug31.rds - Hen-with-brood predictions\n")
cat("  3. pph_df_aug31.rds - Poults per hen (PPB × HWB)\n")
cat("  4. april23_scaled_weather.rds - Weather covariates\n")
cat("\nNext steps:\n")
cat("  - Fit recruitment model: pph ~ weather + density\n")
cat("  - Use predictions to update Pbar in April decision model\n")
cat("  - Compare April vs September decision timing\n")
cat("=======================================\n\n")

################################################################################
# END OF SCRIPT
################################################################################