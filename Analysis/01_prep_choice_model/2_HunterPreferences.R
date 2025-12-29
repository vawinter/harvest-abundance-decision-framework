# Turkey Hunting Season Length Utility Extrapolation
# Based on choice model results & DRD modeling

rm(list=ls())
gc()
#---------------X
library(brglm2)
library(dplyr)
library(ggplot2)
library(cowplot)

# Set scenario
scenario <- "C"

# set colors for easy consistent, easy visualization
length <- c("0" = "#e41a1c", "1" = "#4daf4a", "2" = "#377eb8", "3" = "#984ea3")

# Load and prepare data
df <- read.csv("Data/2024 Fall Turkey Choice Model Dataset - Veronica Model - arc 5.22.2025 CSV.csv")
df$Season.Length <- as.factor(df$Season.Length)
df$Population.Size <- ifelse(df$Population.Size==1,"Decrease",
                              ifelse(df$Population.Size==2,"Stable","Increase"))
df$Population.Size <- as.factor(df$Population.Size)
df$Population.Size <- factor(df$Population.Size, levels=(c('Decrease','Stable','Increase')))

# Calculate summary statistics by population and season length
out <- df %>% group_by(Population.Size,Season.Length) %>%
              summarize(mn = mean(Response.Indicator))
out
# Plot un-modeled data
start <- ggplot(data=out) +
  geom_line(aes(x=Population.Size, y=mn, group=Season.Length, color=Season.Length)) +
  geom_point(aes(x=Population.Size, y=mn, group=Season.Length, color=Season.Length), size = 2) +
  scale_color_manual(values = length) +
  labs(x="Population Trend", y="Proportion of times the preferred option", color = "") +
  ylim(0, 1) +
  theme_classic()+
  theme(legend.position = "top")

# Next: Need to find reasonable values for missing weeks within population sizes
# First, need to expand the data set to include all possible weeks within pop. sizes
# Expand 'out'----

# Define all possible combinations
full_data <- expand.grid(
  Population.Size = factor(c("Decrease", "Stable", "Increase")),
  Season.Length = factor(0:3)
)

# Join and fill in missing values
out_full <- full_data %>%
  left_join(out, by = c("Population.Size", "Season.Length")) %>%
  # Create a new column that fills in the missing data w. values
  mutate(mn2 = case_when(## Scenario A: Even split in hunters wanting a closed season (0) and 1 week
                          # Population.Size == "Decrease" & Season.Length == "0" ~ 0.4,
                          # Population.Size == "Decrease" & Season.Length == "1" ~ 0.4,
                         # Scenario B: 3:1 hunters want at least 1 week even if population is declining
                         Population.Size == "Decrease" & Season.Length == "0" ~ 0.3,
                         Population.Size == "Decrease" & Season.Length == "1" ~ 0.6,
                         # ## Scenario C: Closed season, and forego a year of hunting
                         # Population.Size == "Decrease" & Season.Length == "0" ~ 0.8,
                         # Population.Size == "Decrease" & Season.Length == "1" ~ 0.3,
    
                         # Others that need to stay regardless of A/B/C
                         Population.Size == "Stable" & Season.Length == "0" ~ 0,
                         Population.Size == "Increase" & Season.Length == "0" ~ 0,
                         Population.Size == "Increase" & Season.Length == "3" ~ 0.87,
                         TRUE ~ mn))
out_full$Population.Size <- as.factor(out_full$Population.Size)
out_full$Population.Size <- factor(out_full$Population.Size, levels=(c('Decrease','Stable','Increase')))
out_full

# Plot expanded data set
mid <- ggplot(data=out_full) +
  geom_point(aes(x=Population.Size, y=mn2, group=Season.Length, color=Season.Length), size = 2) +
  geom_line(aes(x=Population.Size, y=mn2, group=Season.Length, color=Season.Length)) +
  scale_color_manual(values = length) +
  labs(x="Population Trend", y="Proportion of times the preferred option", color = "") +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "top")


out_full %>%
  filter(Population.Size == "Decrease") %>%
  ggplot(aes(x = Season.Length, y = mn2)) +
  geom_point(size = 4) +
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  labs(x = "Weeks",
       y = "Peference",
       title = "Preference of Weeks When Population is DECREASING") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),       # tick labels
    axis.title = element_text(size = 14),      # axis titles
    legend.text = element_text(size = 12),     # legend labels (if you add one later)
    legend.title = element_text(size = 13),     # legend title
    plot.title = element_text(hjust = 0.5, size = 16)  # centers the title
  )


# plot together
plot_grid(start, mid)

# Now, normalize all values within a population size group
out_norm <- out_full %>% 
  # Group by sizes
  group_by(Population.Size) %>% 
  mutate(
    # Sum the raw weights
    sum_weights = sum(mn2), 
    # Normalized weights within a group
    normalized_weights = round(mn2/sum_weights, 2),
    # check everything sums to 1
   sum = sum(round(normalized_weights, 2))
     ) %>% 
  ungroup() %>% 
  # remove extra columns
  dplyr::select(-c(mn, sum_weights))
out_norm

# Remove check column post checking
out_norm$sum <- NULL

# Plot
end <- ggplot(data=out_norm) +
  geom_line(aes(x=Population.Size, y=normalized_weights, group=Season.Length, color=Season.Length)) +
  geom_point(aes(x=Population.Size, y=normalized_weights, group=Season.Length, color=Season.Length), size = 2) +
  scale_color_manual(values = length) +
  labs(x="Population Trend", y="Proportion of times the preferred option", color = "") +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "top")

# plot together
plot_grid(mid, end)

# Save weights
write.csv(out_norm, paste0("Data/norm_weights_popsize_scenario_", scenario, ".csv"),
          row.names = F)

# Done!
out_normA <- read.csv("Data/norm_weights_popsize_scenario_A.csv")
xA <- out_normA %>% 
  arrange(Population.Size)
xA

out_normB <- read.csv("Data/norm_weights_popsize_scenario_B.csv")
xB <- out_normB %>% 
  arrange(Population.Size)
xB

out_normC <- read.csv("Data/norm_weights_popsize_scenario_C.csv")
xC <- out_norm %>% 
  arrange(Population.Size) 
xC


#-------_X
library(dplyr)
library(tidyr)

# Read and prepare all three scenarios
out_normA <- read.csv("Data/norm_weights_popsize_scenario_A.csv") %>% 
  arrange(Population.Size)  
out_normB <- read.csv("Data/norm_weights_popsize_scenario_B.csv")%>% 
  arrange(Population.Size) 
out_normC <- read.csv("Data/norm_weights_popsize_scenario_C.csv")%>% 
  arrange(Population.Size) 

# Combine all scenarios
combined <- bind_rows(
  out_normA %>% mutate(Scenario = "A"),
  out_normB %>% mutate(Scenario = "B"),
  out_normC %>% mutate(Scenario = "C")
) %>%
  select(Population.Size, Scenario, Season.Length, normalized_weights) %>%
  arrange(Population.Size, Scenario, Season.Length)

# View the table
print(combined)

# Optional: Create a wider format for easier reading
table_wide <- combined %>%
  pivot_wider(
    names_from = Season.Length,
    values_from = normalized_weights,
    names_prefix = "Week_"
  ) %>%
  arrange(Population.Size, Scenario)

print(table_wide)


# # Survival
# x <- readRDS("Data/IPM/2024/O_24_comb-survival_summary.rds") %>% 
#   group_by(Region, sex, age_class) %>% 
#   # mean survival
#   summarise(avg_surv = mean(median_value, na.rm = T)) %>% 
#   # grab ROI
#   filter(Region %in% c("Region 5", "Region 6"))
# x

# Save as CSV
write.csv(combined, "manuscript_weights_table.csv", row.names = FALSE)

# Or create a publication-ready table using knitr
library(knitr)
kable(table_wide, 
      col.names = c("Population Trend", "Scenario", "0 weeks", "1 week", "2 weeks", "3 weeks"),
      digits = 2,
      caption = "Normalized utility weights by population trend, scenario, and season length")
library(clipr)

# Create your table
table_for_word <- table_wide %>%
  rename(
    "Population Trend" = Population.Size,
    "0 weeks" = Week_0,
    "1 week" = Week_1,
    "2 weeks" = Week_2,
    "3 weeks" = Week_3
  ) %>%
  mutate(across(contains("week"), ~round(., 2)))

# Copy to clipboard
write_clip(table_for_word)
