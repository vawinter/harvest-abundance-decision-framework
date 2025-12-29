################################################################################
# Calculate Utility Weights from Choice Model
################################################################################
# 
# Purpose: Convert raw choice model preferences into normalized utility weights
#          for use in the Markov Decision Process model. Creates three scenarios
#          reflecting uncertainty in hunter preferences under declining populations.
#
# Input:  2024 Fall Turkey Choice Model Dataset - Veronica Model - arc 5.22.2025 CSV.csv
#         (Raw choice experiment data)
#
# Output: norm_weights_popsize_scenario_A.csv (Even split: 0 vs 1 week)
#         norm_weights_popsize_scenario_B.csv (3:1 prefer ≥1 week)
#         norm_weights_popsize_scenario_C.csv (Unanimous closure)
#
# Author: Veronica A. Winter
# Date: December 2025
################################################################################

rm(list = ls())
gc()

library(brglm2)
library(dplyr)
library(ggplot2)
library(cowplot)

################################################################################
# 1. SET SCENARIO AND VISUALIZATION PARAMETERS
################################################################################

# Which scenario to calculate
# Options: "A", "B", or "C"
scenario <- "B"

# Color scheme for consistent visualization across season lengths
length_colors <- c("0" = "#e41a1c",  # Red (closed season)
                   "1" = "#4daf4a",  # Green (1 week)
                   "2" = "#377eb8",  # Blue (2 weeks)
                   "3" = "#984ea3")  # Purple (3 weeks)

################################################################################
# 2. LOAD AND PREPARE CHOICE MODEL DATA
################################################################################

# Load raw choice experiment data
df <- read.csv("Data/2024 Fall Turkey Choice Model Dataset - Veronica Model - arc 5.22.2025 CSV.csv")

# Convert variables to factors
df$Season.Length <- as.factor(df$Season.Length)

# Recode population size from numeric to interpretable labels
df$Population.Size <- ifelse(df$Population.Size == 1, "Decrease",
                             ifelse(df$Population.Size == 2, "Stable", "Increase"))
df$Population.Size <- as.factor(df$Population.Size)

# Order factor levels logically
df$Population.Size <- factor(df$Population.Size, 
                             levels = c('Decrease', 'Stable', 'Increase'))

################################################################################
# 3. CALCULATE RAW PREFERENCE PROPORTIONS
################################################################################

# Calculate proportion of times each season length was chosen
# under each population trend
out <- df %>% 
  group_by(Population.Size, Season.Length) %>%
  summarize(mn = mean(Response.Indicator), .groups = 'drop')

print("Raw preference proportions:")
print(out)

# Plot raw (unmodeled) data
plot_raw <- ggplot(data = out) +
  geom_line(aes(x = Population.Size, y = mn, 
                group = Season.Length, color = Season.Length),
            linewidth = 1) +
  geom_point(aes(x = Population.Size, y = mn, 
                 group = Season.Length, color = Season.Length), 
             size = 3) +
  scale_color_manual(values = length_colors) +
  labs(x = "Population Trend", 
       y = "Proportion of Times Preferred", 
       color = "Season Length\n(weeks)",
       title = "Raw Survey Data") +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

################################################################################
# 4. FILL IN MISSING COMBINATIONS (SCENARIO-SPECIFIC)
################################################################################

# The choice experiment did not include all combinations
# Here we fill in missing values based on three alternative scenarios
# that reflect uncertainty about hunter preferences under declining populations

# Create full grid of all possible combinations
full_data <- expand.grid(
  Population.Size = factor(c("Decrease", "Stable", "Increase"),
                           levels = c("Decrease", "Stable", "Increase")),
  Season.Length = factor(0:3)
)

# Join with observed data and fill missing values
out_full <- full_data %>%
  left_join(out, by = c("Population.Size", "Season.Length")) %>%
  mutate(mn2 = case_when(
    
    # ==== DECLINING POPULATION SCENARIOS (differs by A/B/C) ====
    
    # Scenario A: Even split between 0 and 1 week
    # Reflects moderate uncertainty about conservation concern
    # Population.Size == "Decrease" & Season.Length == "0" ~ 0.4,
    # Population.Size == "Decrease" & Season.Length == "1" ~ 0.4,
    
    # Scenario B: 3:1 ratio favoring at least 1 week
    # Reflects strong hunting tradition even during decline
    Population.Size == "Decrease" & Season.Length == "0" ~ 0.3,
    Population.Size == "Decrease" & Season.Length == "1" ~ 0.6,
    
    # Scenario C: Strong preference for closure
    # Reflects high conservation priority
    # Population.Size == "Decrease" & Season.Length == "0" ~ 0.8,
    # Population.Size == "Decrease" & Season.Length == "1" ~ 0.3,
    
    # ==== STABLE AND INCREASING POPULATIONS (same across scenarios) ====
    
    # Under stable populations: no support for closed season
    Population.Size == "Stable" & Season.Length == "0" ~ 0,
    
    # Under increasing populations: no support for closed season
    Population.Size == "Increase" & Season.Length == "0" ~ 0,
    
    # Under increasing populations: strong preference for 3 weeks
    Population.Size == "Increase" & Season.Length == "3" ~ 0.87,
    
    # All other combinations: use observed data
    TRUE ~ mn
  ))

# Ensure factor levels are preserved
out_full$Population.Size <- factor(out_full$Population.Size, 
                                   levels = c('Decrease', 'Stable', 'Increase'))

print("Filled preference proportions:")
print(out_full)

# Plot expanded dataset
plot_filled <- ggplot(data = out_full) +
  geom_point(aes(x = Population.Size, y = mn2, 
                 group = Season.Length, color = Season.Length), 
             size = 3) +
  geom_line(aes(x = Population.Size, y = mn2, 
                group = Season.Length, color = Season.Length),
            linewidth = 1) +
  scale_color_manual(values = length_colors) +
  labs(x = "Population Trend", 
       y = "Proportion of Times Preferred", 
       color = "Season Length\n(weeks)",
       title = paste("Scenario", scenario, "- Filled Data")) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

# Focus plot: declining population preferences
plot_decline <- out_full %>%
  filter(Population.Size == "Decrease") %>%
  ggplot(aes(x = Season.Length, y = mn2)) +
  geom_point(size = 4) +
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  labs(x = "Season Length (weeks)",
       y = "Preference",
       title = "Preferences When Population is DECREASING") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16))

# Show comparison
plot_grid(plot_raw, plot_filled, ncol = 2)

################################################################################
# 5. NORMALIZE WEIGHTS WITHIN EACH POPULATION TREND
################################################################################

# Normalize so weights sum to 1.0 within each population trend
# This ensures proper probability interpretation in the MDP model
out_norm <- out_full %>% 
  group_by(Population.Size) %>% 
  mutate(
    # Sum of raw weights within group
    sum_weights = sum(mn2), 
    # Normalized weights (sum to 1.0)
    normalized_weights = round(mn2 / sum_weights, 2),
    # Verify sum (should be 1.0)
    sum_check = sum(round(normalized_weights, 2))
  ) %>% 
  ungroup() %>% 
  # Keep only necessary columns
  select(-c(mn, sum_weights))

print("Normalized weights:")
print(out_norm)

# Verify normalization
print("Sum checks (should all be 1.0):")
print(out_norm %>% group_by(Population.Size) %>% summarize(sum = sum(normalized_weights)))

# Remove check column
out_norm$sum_check <- NULL

# Plot normalized weights
plot_normalized <- ggplot(data = out_norm) +
  geom_line(aes(x = Population.Size, y = normalized_weights, 
                group = Season.Length, color = Season.Length),
            linewidth = 1) +
  geom_point(aes(x = Population.Size, y = normalized_weights, 
                 group = Season.Length, color = Season.Length), 
             size = 3) +
  scale_color_manual(values = length_colors) +
  labs(x = "Population Trend", 
       y = "Normalized Weight", 
       color = "Season Length\n(weeks)",
       title = paste("Scenario", scenario, "- Normalized Weights")) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

# Show comparison of filled vs normalized
plot_grid(plot_filled, plot_normalized, ncol = 2)

################################################################################
# 6. SAVE NORMALIZED WEIGHTS FOR MDP MODEL
################################################################################

# Save to CSV for use in MATLAB decision model
output_file <- paste0("Data/norm_weights_popsize_scenario_", scenario, ".csv")
write.csv(out_norm, output_file, row.names = FALSE)

cat("\nNormalized weights saved to:", output_file, "\n")

################################################################################
# 7. COMPARE ALL THREE SCENARIOS (OPTIONAL)
################################################################################

# Read all three scenario files
out_normA <- read.csv("Data/norm_weights_popsize_scenario_A.csv") %>% 
  arrange(Population.Size)

out_normB <- read.csv("Data/norm_weights_popsize_scenario_B.csv") %>% 
  arrange(Population.Size)

out_normC <- read.csv("Data/norm_weights_popsize_scenario_C.csv") %>% 
  arrange(Population.Size)

# Combine all scenarios for comparison
combined <- bind_rows(
  out_normA %>% mutate(Scenario = "A"),
  out_normB %>% mutate(Scenario = "B"),
  out_normC %>% mutate(Scenario = "C")
) %>%
  select(Population.Size, Scenario, Season.Length, normalized_weights) %>%
  arrange(Population.Size, Scenario, Season.Length)

cat("\nCombined scenario comparison:\n")
print(combined)

# Create wider format for easier reading
library(tidyr)
table_wide <- combined %>%
  pivot_wider(
    names_from = Season.Length,
    values_from = normalized_weights,
    names_prefix = "Week_"
  ) %>%
  arrange(Population.Size, Scenario)

cat("\nWide format table:\n")
print(table_wide)

# Create publication-ready table
library(knitr)
table_publication <- kable(
  table_wide, 
  col.names = c("Population Trend", "Scenario", 
                "0 weeks", "1 week", "2 weeks", "3 weeks"),
  digits = 2,
  caption = "Normalized utility weights by population trend, scenario, and season length"
)

print(table_publication)

# Save combined table
write.csv(combined, "Data/manuscript_weights_table.csv", row.names = FALSE)

# Copy to clipboard for Word (optional - requires clipr package)
# library(clipr)
# table_for_word <- table_wide %>%
#   rename(
#     "Population Trend" = Population.Size,
#     "0 weeks" = Week_0,
#     "1 week" = Week_1,
#     "2 weeks" = Week_2,
#     "3 weeks" = Week_3
#   ) %>%
#   mutate(across(contains("week"), ~round(., 2)))
# write_clip(table_for_word)

################################################################################
# 8. KEY FINDINGS SUMMARY
################################################################################

cat("\n========== SCENARIO DEFINITIONS ==========\n")
cat("Scenario A: Even split (50% prefer 0 weeks, 50% prefer 1 week)\n")
cat("Scenario B: 3:1 ratio (75% prefer ≥1 week even if declining)\n")
cat("Scenario C: Unanimous closure (80% support 0 weeks)\n")
cat("\nNote: Scenarios only differ under DECREASING populations.\n")
cat("Under Stable and Increasing populations, all scenarios are identical.\n")
cat("==========================================\n\n")

################################################################################
# END OF SCRIPT
################################################################################