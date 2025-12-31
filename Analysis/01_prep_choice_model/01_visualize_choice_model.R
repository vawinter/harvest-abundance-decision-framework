################################################################################
# Hunter Preference Choice Model - Data Visualization
################################################################################
# 
# Purpose: Visualize hunter preferences for season length under different 
#          population trends (Increasing, Stable, Decreasing) from discrete 
#          choice experiment survey
#
# Input:  2024 Fall Turkey Choice Model Dataset - Veronica Model - arc 5.22.2025 CSV.csv
#         (Survey data with Response.Indicator, Population.Size, Season.Length)
#
# Output: Line plot showing proportion of times each season length was 
#         preferred under each population trend
#
# Author: Duane R. Diefenbach (original)
#         Veronica A. Winter (modifications)
# Date: May 2025
################################################################################
rm(list=ls())
gc()

library(brglm2)
library(dplyr)
library(ggplot2)

################################################################################
# 1. Load and Prepare Data
################################################################################

# Load choice model dataset
# Response.Indicator: 1 if option was chosen, 0 otherwise
# Population.Size: 1=Decrease, 2=Stable, 3=Increase
# Season.Length: Number of weeks (0, 1, 2, or 3)
df <- read.csv("Data/2024 Fall Turkey Choice Model Dataset.csv")

################################################################################
# 2. Recode Variables
################################################################################

# Convert Season.Length to factor for categorical analysis
df$Season.Length <- as.factor(df$Season.Length)

# Recode Population.Size from numeric to interpretable labels
df$Population.Size <- ifelse(df$Population.Size == 1, "Decrease",
                             ifelse(df$Population.Size == 2, "Stable", "Increase"))
df$Population.Size <- as.factor(df$Population.Size)

# Order population size factor levels logically
df$Population.Size <- factor(df$Population.Size, 
                             levels = c('Decrease', 'Stable', 'Increase'))

################################################################################
# 3. Calculate Mean Preference by Population Trend and Season Length
################################################################################

# Calculate proportion of times each season length was preferred
# under each population trend
out <- df %>% 
  group_by(Population.Size, Season.Length) %>%
  summarize(mn = mean(Response.Indicator), .groups = 'drop')

# Display results
print(out)

################################################################################
# 4. Create Visualization
################################################################################

# Plot preference patterns
# x-axis: Population trend (Decrease, Stable, Increase)
# y-axis: Proportion of times preferred (0-1)
# Lines: Different season lengths (0, 1, 2, 3 weeks)
ggplot(data = out) +
  geom_line(aes(x = Population.Size, 
                y = mn, 
                group = Season.Length, 
                color = Season.Length),
            linewidth = 1.2) +
  labs(x = "Population Trend", 
       y = "Proportion of Times Preferred",
       color = "Season Length\n(weeks)",
       title = "Hunter Preferences for Season Length by Population Trend") +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

################################################################################
# 5. Key Findings
################################################################################

# Interpretation:
# - Under decreasing populations: Hunters strongly prefer shorter seasons
# - Under stable populations: Moderate preference for 2-3 week seasons
# - Under increasing populations: Strong preference for 3-week seasons
#
# These preferences are expanded in:
#   02_calculate_preferences.R
################################################################################