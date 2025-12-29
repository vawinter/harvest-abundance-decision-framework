##############################################################################X
#                 X--------# Chapter 2: MDP for PA Turkey #-------X
#
#     X--------# Goal: Testing Population Model from PennTurkeyModel #-------X
#                     Function in MatLab from P. Fackler 
#
# Created: 11/11/2025
#
# V. Winter
##############################################################################X
# clean env
rm(list=ls())
gc()

# Libraries
library(tidyverse)

##----------------------------------X
# Construct the function ----
##----------------------------------X
# Population model function
turkey_pop_model <- function(
  # Initial populations (January/December)
  Md_init, Jd_init, Fd_init,
  # Season and environment
  season_length = 2,  # 0-3 weeks
  mast_level = 12.998,  # oak mast
  # Fecundity parameters
  Fbar = 2.5,
  Pbar = 1.8,
  slope = 0.2,
  # Survival parameters
  gammam = 0.409,
  gammaj = 0.653,
  gammaf = 0.542,
  # Stochasticity (set to 1 for deterministic)
  v = 1,  # summer survival multiplier
  w = 1   # poult production multiplier
) {
  
  # Calculate fecundity shape parameter
  eta <- 1 + 2 * Fbar * slope / Pbar
  
  # Calculate harvest rate based on season length and mast
  harvest_rate <- function(L, O) {
    base_rate <- 0.07 - 0.00686 * O + 0.0175 * L
    max(0.01, base_rate)
  }
  
  H <- harvest_rate(season_length, mast_level)
  
  # Summer transitions (January -> September)
  Ms <- gammam * Md_init * v
  Js <- gammaj * Jd_init * v
  Fs <- gammaf * Fd_init * v
  
  # Density-dependent poult production
  Ps <- (2 * eta * Pbar / (eta + 1 + (eta - 1) * (Fd_init / Fbar)^eta)) * w
  
  # Winter transitions (September -> December next year)
  Md_next <- Ms + Js
  Jd_next <- (1 - H) * Fs * Ps / 2
  Fd_next <- (1 - H) * Fs * (1 + Ps / 2)
  
  # Return all values
  list(
    # September populations
    Ms = Ms, Js = Js, Fs = Fs, Ps = Ps,
    # Harvest
    harvest_rate = H,
    # December populations (next year)
    Md_next = Md_next,
    Jd_next = Jd_next,
    Fd_next = Fd_next,
    # Summary
    total_males_sept = Ms + Js,
    total_males_dec = Md_next,
    total_pop_dec = Md_next + Jd_next + Fd_next
  )
}

##----------------------------------X
# Multi-year projection function ---
##----------------------------------X
project_turkey_pop_model <- function(years, 
                                     Md_init, 
                                     Jd_init, 
                                     Fd_init,
                                     season_length = 2,
                                     mast_level,
                                     Fbar, Pbar, slope) {
  
  # Storage
  trajectory <- tibble(
    year = 0:years,
    Md = numeric(years + 1),
    Jd = numeric(years + 1),
    Fd = numeric(years + 1),
    harvest_rate = numeric(years + 1),
    Ps = numeric(years + 1)
  )
  
  # Initial conditions
  trajectory$Md[1] <- Md_init
  trajectory$Jd[1] <- Jd_init
  trajectory$Fd[1] <- Fd_init
  
  # Project forward
  for (i in 1:years) {
    result <- turkey_pop_model(
      Md_init = trajectory$Md[i],
      Jd_init = trajectory$Jd[i],
      Fd_init = trajectory$Fd[i],
      season_length = season_length,
      mast_level = mast_level,
      Fbar = Fbar, Pbar = Pbar, slope = slope
    )
    
    trajectory$Md[i + 1] <- result$Md_next
    trajectory$Jd[i + 1] <- result$Jd_next
    trajectory$Fd[i + 1] <- result$Fd_next
    trajectory$harvest_rate[i] <- result$harvest_rate
    trajectory$Ps[i] <- result$Ps
  }
  
  trajectory
}

##----------------------------------X
# Example: Run single time step ---
##----------------------------------X
result <- turkey_pop_model(
  Md_init = 1.5,
  Jd_init = 1.5,
  Fd_init = 3.0,
  season_length = 2)

print(result)
##----------------------------------X
# Example: 20-year projection ---
##----------------------------------X
## Set up inputs ----
##----------------------------------X
# Project out
pop_trajectory <- project_turkey_pop_model(years = 20, 
                                           Md_init = 1.5, 
                                            Jd_init = 1.5, 
                                            Fd_init = 3,
                                            season_length = 2, slope = 0.2,
                                            Fbar = 3, Pbar = 1.2, mast_level= 23)

# Plot results
pop_trajectory %>%
  pivot_longer(cols = c(Md, Jd, Fd), names_to = "class", values_to = "density") %>%
  ggplot(aes(x = year, y = density, color = class)) +
  scale_color_manual(values = c(
    "Fd"   = "#71130F",  # Color for adult female
    "Md"   = "#043A50",
    "Jd"   = "#5C7391"
  )) +
  geom_line(size = 0.5) +
  geom_point(size = 3) +
  labs(title = paste("Turkey Population Projection: ", length(pop_trajectory$year), "year projection", sep = " "),
       x = "Year", y = "Density",
       color = "Age/Sex Class") +
  theme_classic()

##----------------------------------X
# Compare different season lengths
map_df(0:3, function(L) {
  project_turkey_pop_model(years = 20, 
                           Md_init = 1.5, 
                           Jd_init = 1.5, 
                           Fd_init = 3.0,
                           season_length = L,
                           mast_level = 23,
                           Fbar = 2.75, Pbar = 1.98, slope = 0.2) %>%
    mutate(season_length = L)
}) %>%
  # plot output
  ggplot(aes(x = year, y = Fd, color = factor(season_length))) +
  scale_color_manual(values = c(
                    "0" = "#A05C4E", 
                    "1" = "#F1D2CA",  
                    "2" ="#8B9E90",  
                    "3" ="#2B4648")
                )  +
  geom_line(size = 0.5) +
  geom_point(size = 3) +
  labs(title = "Female Density by Season Length",
       x = "Year", y = "Female Density",
       color = "Season Length (weeks)") +
  theme_classic()

###--------------------------X
## 10% decrease
1.8-(1.8*0.10)
# Pbar = 1.62
2.5-(2.5*0.10)
# Fbar = 2.25
###------------X
# Stable
# Pbar = 1.8
# Fbar = 2.5
###------------X
## 10% increase
1.8+(1.8*0.10)
# Pbar = 1.98
2.5+(2.5*0.10)
# Fbar = 2.75
###--------------------------X

