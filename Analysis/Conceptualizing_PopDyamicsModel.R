################################################################################
# Conceptual R Implementation of PennTurkeyModel (MATLAB)
################################################################################
#
# Purpose: 
#   Test and visualize the turkey population dynamics model before 
#   implementing in MATLAB MDP framework. This R version mirrors the
#   biological model structure in PennTurkeyModel.m
#
# Model Structure:
#   - State variables: Males (M), Jakes (J), Females (F) at December
#   - Transitions: December → September → December (annual cycle)
#   - Decision: Season length (0-3 weeks)
#   - Stochasticity: Survival (v), Fecundity (w), Mast (O)
#
# Use Cases:
#   1. Verify biological model logic before MATLAB implementation
#   2. Explore parameter sensitivity
#   3. Visualize population trajectories
#   4. Compare season length effects
#
# Author: Veronica A. Winter
# Date: November 2025
################################################################################

rm(list = ls())
gc()

library(tidyverse)

################################################################################
# 1. CORE POPULATION MODEL FUNCTION
################################################################################

turkey_pop_model <- function(
    # Initial populations (December t)
  Md_init,          # Male density
  Jd_init,          # Jake density  
  Fd_init,          # Female density
  
  # Management decision
  season_length = 2,  # Weeks of fall harvest (0, 1, 2, or 3)
  
  # Environmental conditions
  mast_level = 12.998,  # Oak mast abundance
  
  # Fecundity parameters (density-dependent reproduction)
  Fbar = 2.5,       # Female density at inflection point
  Pbar = 1.8,       # Poults per hen at inflection point
  slope = 0.2,      # Slope of fecundity curve at inflection
  
  # Survival parameters (spring/summer)
  gammam = 0.409,   # Male survival rate (Jan-Aug)
  gammaj = 0.653,   # Jake survival rate (Jan-Aug)
  gammaf = 0.542,   # Female survival rate (Jan-Aug)
  
  # Stochastic multipliers (set to 1 for deterministic)
  v = 1,            # Summer survival noise
  w = 1             # Poult production noise
) {
  
  #---------------------------------------------------------------------------
  # Calculate fecundity shape parameter
  #---------------------------------------------------------------------------
  # eta controls the shape of the density-dependent fecundity curve
  eta <- 1 + 2 * Fbar * slope / Pbar
  
  #---------------------------------------------------------------------------
  # Calculate harvest rate
  #---------------------------------------------------------------------------
  # Harvest depends on season length (L) and mast abundance (O)
  # More mast → lower harvest (turkeys dispersed)
  # Longer season → higher harvest
  harvest_rate <- function(L, O) {
    base_rate <- 0.07 - 0.00686 * O + 0.0175 * L
    max(0.01, base_rate)  # Minimum 1% harvest
  }
  
  H <- harvest_rate(season_length, mast_level)
  
  #---------------------------------------------------------------------------
  # Summer transitions (December t → September t)
  #---------------------------------------------------------------------------
  # Adults survive spring/summer with category-specific rates
  Ms <- gammam * Md_init * v  # Males (spring breeding, predation)
  Js <- gammaj * Jd_init * v  # Jakes (less vulnerable than adults)
  Fs <- gammaf * Fd_init * v  # Females (nesting mortality)
  
  #---------------------------------------------------------------------------
  # Density-dependent poult production
  #---------------------------------------------------------------------------
  # Reverse sigmoidal function: high production at low density
  # P(F) = [2η × P̄] / [η + 1 + (η-1) × (F/F̄)^η]
  Ps <- (2 * eta * Pbar / (eta + 1 + (eta - 1) * (Fd_init / Fbar)^eta)) * w
  
  #---------------------------------------------------------------------------
  # Winter transitions (September t → December t)
  #---------------------------------------------------------------------------
  # Post-harvest survival (hunting season is Sept-Nov)
  Md_next <- Ms + Js                        # Males: adults + recruited jakes
  Jd_next <- (1 - H) * Fs * Ps / 2         # Jakes: male poults that survive harvest
  Fd_next <- (1 - H) * Fs * (1 + Ps / 2)   # Females: adults + female poults
  
  #---------------------------------------------------------------------------
  # Return results
  #---------------------------------------------------------------------------
  list(
    # September populations
    Ms = Ms,
    Js = Js,
    Fs = Fs,
    Ps = Ps,
    
    # Harvest
    harvest_rate = H,
    
    # December populations (year t+1)
    Md_next = Md_next,
    Jd_next = Jd_next,
    Fd_next = Fd_next,
    
    # Summary metrics
    total_males_sept = Ms + Js,
    total_males_dec = Md_next,
    total_pop_dec = Md_next + Jd_next + Fd_next
  )
}

################################################################################
# 2. MULTI-YEAR PROJECTION FUNCTION
################################################################################

project_turkey_pop <- function(
    years,                # Number of years to project
    Md_init,              # Initial male density
    Jd_init,              # Initial jake density
    Fd_init,              # Initial female density
    season_length = 2,    # Season length (constant across years)
    mast_level = 12.998,  # Mast level (constant or vector)
    Fbar = 2.5,           # Fecundity parameters
    Pbar = 1.8,
    slope = 0.2
) {
  
  # Initialize storage
  trajectory <- tibble(
    year = 0:years,
    Md = numeric(years + 1),
    Jd = numeric(years + 1),
    Fd = numeric(years + 1),
    harvest_rate = numeric(years + 1),
    Ps = numeric(years + 1)
  )
  
  # Set initial conditions
  trajectory$Md[1] <- Md_init
  trajectory$Jd[1] <- Jd_init
  trajectory$Fd[1] <- Fd_init
  
  # Project forward year by year
  for (i in 1:years) {
    result <- turkey_pop_model(
      Md_init = trajectory$Md[i],
      Jd_init = trajectory$Jd[i],
      Fd_init = trajectory$Fd[i],
      season_length = season_length,
      mast_level = mast_level,
      Fbar = Fbar, 
      Pbar = Pbar, 
      slope = slope
    )
    
    # Store next year's populations
    trajectory$Md[i + 1] <- result$Md_next
    trajectory$Jd[i + 1] <- result$Jd_next
    trajectory$Fd[i + 1] <- result$Fd_next
    trajectory$harvest_rate[i] <- result$harvest_rate
    trajectory$Ps[i] <- result$Ps
  }
  
  return(trajectory)
}

################################################################################
# 3. EXAMPLE: SINGLE TIME STEP
################################################################################

cat("\n========== SINGLE TIME STEP EXAMPLE ==========\n")

result <- turkey_pop_model(
  Md_init = 1.5,
  Jd_init = 1.5,
  Fd_init = 3.0,
  season_length = 2
)

cat("Initial populations (December):\n")
cat("  Males:", 1.5, "  Jakes:", 1.5, "  Females:", 3.0, "\n\n")

cat("September populations:\n")
cat("  Males:", round(result$Ms, 3), 
    "  Jakes:", round(result$Js, 3), 
    "  Females:", round(result$Fs, 3), "\n")
cat("  Poults per hen:", round(result$Ps, 3), "\n\n")

cat("Harvest rate:", round(result$harvest_rate, 3), "\n\n")

cat("Next December populations:\n")
cat("  Males:", round(result$Md_next, 3), 
    "  Jakes:", round(result$Jd_next, 3), 
    "  Females:", round(result$Fd_next, 3), "\n")

################################################################################
# 4. EXAMPLE: 20-YEAR PROJECTION
################################################################################

cat("\n========== 20-YEAR PROJECTION ==========\n")

pop_trajectory <- project_turkey_pop(
  years = 20, 
  Md_init = 1.5, 
  Jd_init = 1.5, 
  Fd_init = 3.0,
  season_length = 2, 
  mast_level = 12.998,
  Fbar = 2.5, 
  Pbar = 1.8, 
  slope = 0.2
)

# Plot trajectory
trajectory_plot <- pop_trajectory %>%
  pivot_longer(cols = c(Md, Jd, Fd), 
               names_to = "class", 
               values_to = "density") %>%
  mutate(class = factor(class, 
                        levels = c("Md", "Jd", "Fd"),
                        labels = c("Males", "Jakes", "Females"))) %>%
  ggplot(aes(x = year, y = density, color = class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "Males" = "#043A50",
      "Jakes" = "#5C7391",
      "Females" = "#71130F"
    )
  ) +
  labs(
    title = "Turkey Population Projection: 20-Year Trajectory",
    subtitle = "Season length = 2 weeks, Baseline parameters",
    x = "Year", 
    y = "Density (birds per unit area)",
    color = "Age/Sex Class"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

print(trajectory_plot)
ggsave("Dataviz/turkey_trajectory_baseline.png", 
       trajectory_plot, width = 8, height = 5, dpi = 300)

################################################################################
# 5. COMPARE SEASON LENGTHS
################################################################################

cat("\n========== COMPARING SEASON LENGTHS ==========\n")

season_comparison <- map_df(0:3, function(L) {
  project_turkey_pop(
    years = 20, 
    Md_init = 1.5, 
    Jd_init = 1.5, 
    Fd_init = 3.0,
    season_length = L,
    mast_level = 12.998,
    Fbar = 2.5, 
    Pbar = 1.8, 
    slope = 0.2
  ) %>%
    mutate(season_length = L)
})

# Plot female density by season length
season_plot <- season_comparison %>%
  ggplot(aes(x = year, y = Fd, color = factor(season_length))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "0" = "#A05C4E",   # No hunting
      "1" = "#F1D2CA",   # 1 week
      "2" = "#8B9E90",   # 2 weeks
      "3" = "#2B4648"    # 3 weeks
    ),
    name = "Season Length\n(weeks)"
  ) +
  labs(
    title = "Female Density by Season Length",
    subtitle = "Longer seasons reduce female abundance",
    x = "Year", 
    y = "Female Density"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

print(season_plot)
ggsave("Dataviz/season_length_comparison.png", 
       season_plot, width = 8, height = 5, dpi = 300)

# Summary statistics
season_summary <- season_comparison %>%
  filter(year == 20) %>%
  group_by(season_length) %>%
  summarise(
    Final_Males = round(Md, 3),
    Final_Jakes = round(Jd, 3),
    Final_Females = round(Fd, 3),
    Final_Total = round(Md + Jd + Fd, 3)
  )

cat("\nFinal populations (year 20) by season length:\n")
print(season_summary)

################################################################################
# 6. POPULATION TREND SCENARIOS
################################################################################

cat("\n========== POPULATION TREND SCENARIOS ==========\n")

# Calculate ±10% adjustments from baseline
baseline_Pbar <- 1.8
baseline_Fbar <- 2.5

scenarios <- tribble(
  ~scenario,    ~Pbar,                      ~Fbar,
  "Decrease",   baseline_Pbar * 0.90,       baseline_Fbar * 0.90,
  "Stable",     baseline_Pbar,              baseline_Fbar,
  "Increase",   baseline_Pbar * 1.10,       baseline_Fbar * 1.10
)

cat("Population trend scenarios (±10% from baseline):\n")
print(scenarios)

# Project each scenario
scenario_trajectories <- map_df(1:nrow(scenarios), function(i) {
  project_turkey_pop(
    years = 20,
    Md_init = 1.5,
    Jd_init = 1.5,
    Fd_init = 3.0,
    season_length = 2,
    mast_level = 12.998,
    Fbar = scenarios$Fbar[i],
    Pbar = scenarios$Pbar[i],
    slope = 0.2
  ) %>%
    mutate(scenario = scenarios$scenario[i])
})

# Plot scenarios
scenario_plot <- scenario_trajectories %>%
  ggplot(aes(x = year, y = Fd, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "Decrease" = "#E74C3C",
      "Stable" = "#95A5A6",
      "Increase" = "#27AE60"
    )
  ) +
  labs(
    title = "Female Density by Population Trend",
    subtitle = "Scenarios represent ±10% variation in reproductive capacity",
    x = "Year",
    y = "Female Density",
    color = "Scenario"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

print(scenario_plot)
ggsave("Dataviz/population_trend_scenarios.png",
       scenario_plot, width = 8, height = 5, dpi = 300)

################################################################################
# 7. VERIFY MODEL MATCHES MATLAB
################################################################################

cat("\n========== MODEL VERIFICATION ==========\n")

cat("This R model should produce identical results to MATLAB when:\n")
cat("  1. Same initial conditions\n")
cat("  2. Same parameters (Fbar, Pbar, slope, gammas)\n")
cat("  3. Deterministic mode (v=1, w=1)\n")
cat("  4. Same mast level\n\n")

cat("Use this script to:\n")
cat("  ✓ Verify biological model logic\n")
cat("  ✓ Test parameter sensitivity\n")
cat("  ✓ Visualize population dynamics\n")
cat("  ✓ Debug before MATLAB implementation\n\n")

cat("Next step: Implement in MATLAB with MDPSolve for optimization\n")

################################################################################
# END OF SCRIPT
################################################################################