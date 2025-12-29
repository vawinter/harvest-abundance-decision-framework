################################################################################
# Density-Dependent Fecundity Function for Wild Turkeys
################################################################################
#
# Purpose: Calculate poults per hen (PPH) as a function of female density
#          using a reverse sigmoidal (logistic-like) relationship
#
# Mathematical Form:
#   P(F) = [2η × P̄] / [η + 1 + (η - 1) × (F/F̄)^η]
#   where η = 1 - 2(F̄/P̄) × slope
#
# Biological Interpretation:
#   - At low density (F → 0): PPH approaches maximum (2η/(η+1)) × P̄
#   - At inflection point (F = F̄): PPH = P̄
#   - At high density (F → ∞): PPH → 0 (density-dependent suppression)
#
# References:
#   van Genuchten & Gupta (1993) - Reverse sigmoid function
#   Fackler (2014) - Application to wildlife population dynamics
#
# Author: Veronica A. Winter
# Date: May 2025
################################################################################

calculate_pph <- function(female_density, N_star, Pbar, slope) {
  #' Calculate Density-Dependent Poults Per Hen
  #'
  #' Computes reproductive output (poults per hen) as a function of
  #' female population density using a parameterized reverse sigmoidal curve
  #'
  #' @param female_density Numeric vector. Female density values (birds/km²)
  #' @param N_star Numeric. Female density at inflection point (F̄)
  #' @param Pbar Numeric. PPH value at the inflection point (P̄)
  #' @param slope Numeric. Slope of the curve at the inflection point (r)
  #'
  #' @return Numeric vector of PPH values corresponding to each density
  #'
  #' @details
  #' The function implements a reverse sigmoidal relationship where:
  #'   - High PPH at low density (abundant resources per hen)
  #'   - Declining PPH as density increases (competition, reduced resources)
  #'   - Inflection point occurs at (N_star, Pbar)
  #'
  #' The shape parameter η controls the steepness of the decline.
  #' Larger |slope| values → steeper decline around the inflection point.
  #'
  #' @examples
  #' # Calculate PPH for a range of densities
  #' densities <- seq(0, 10, length.out = 100)
  #' pph <- calculate_pph(densities, N_star = 3, Pbar = 2.5, slope = 0.2)
  #' plot(densities, pph, type = "l")
  
  # Calculate shape parameter η
  # η determines the steepness and asymmetry of the curve
  # Derived from requiring the curve to pass through (N_star, Pbar)
  # with slope 'slope' at that point
  eta <- 1 - 2 * (N_star / Pbar) * slope
  
  # Calculate PPH using reverse sigmoidal function
  # Numerator: Maximum reproductive output scaled by shape
  # Denominator: Controls the rate of decline with density
  pph <- (2 * eta * Pbar) / 
    (eta + 1 + (eta - 1) * (female_density / N_star)^eta)
  
  return(pph)
}

################################################################################
# VALIDATION AND VISUALIZATION
################################################################################

# Only run example if script is sourced interactively
if (interactive()) {
  
  message("\n========== PPH Function Example ==========\n")
  
  # Create sequence of density values for visualization
  density_seq <- seq(0, 10, length.out = 100)
  
  # Example parameters (typical values for PA wild turkeys)
  N_star <- 2.5   # Inflection point density (F̄)
  Pbar <- 2.0     # PPH at inflection (P̄)
  slope <- 0.2    # Growth rate at inflection
  
  # Calculate PPH values
  pph_values <- calculate_pph(
    female_density = density_seq,
    N_star = N_star,
    Pbar = Pbar,
    slope = slope
  )
  
  # Create data frame for plotting
  plot_data <- data.frame(
    female_density = density_seq,
    pph = pph_values
  )
  
  # Mark inflection point
  inflection_point <- data.frame(
    x = N_star, 
    y = Pbar
  )
  
  # Visualize the relationship
  library(ggplot2)
  
  p <- ggplot(plot_data, aes(x = female_density, y = pph)) +
    # Main curve
    geom_line(linewidth = 1.2, color = "steelblue") +
    # Inflection point
    geom_point(data = inflection_point, aes(x = x, y = y), 
               color = "red", size = 4) +
    geom_text(data = inflection_point, aes(x = x, y = y, 
                                           label = paste0("(", round(x, 1), ", ", round(y, 1), ")")),
              vjust = -1, color = "red") +
    # Reference lines
    geom_vline(xintercept = N_star, linetype = "dashed", 
               color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = Pbar, linetype = "dashed", 
               color = "gray50", alpha = 0.5) +
    # Labels
    labs(
      title = "Density-Dependent Reproductive Output",
      subtitle = paste0("F̄ = ", N_star, " birds/km², ",
                        "P̄ = ", Pbar, " poults/hen, ",
                        "slope = ", slope),
      x = expression("Female Density (birds/km"^2*")"),
      y = "Poults Per Hen (PPH)"
    ) +
    # Formatting
    theme_classic() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 10.5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(pph_values) * 1.1))
  
  print(p)
  
  # Print key values
  message("Key values:")
  message(paste("  Max PPH (at F=0):", round(max(pph_values), 2)))
  message(paste("  PPH at inflection:", round(Pbar, 2)))
  message(paste("  PPH at K (F=", round(N_star * 2, 1), "):", 
                round(calculate_pph(N_star * 2, N_star, Pbar, slope), 2)))
  message("\n==========================================\n")
}

################################################################################
# PARAMETER SENSITIVITY TESTING
################################################################################

test_pph_sensitivity <- function() {
  #' Test how PPH function responds to parameter changes
  #'
  #' Creates comparison plots showing effects of varying N_star, Pbar, and slope
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  density_seq <- seq(0, 10, length.out = 100)
  
  # Test varying N_star (inflection point)
  n_star_values <- c(2, 3, 4)
  data_n_star <- expand.grid(
    female_density = density_seq,
    N_star = n_star_values
  ) %>%
    mutate(
      pph = calculate_pph(female_density, N_star, Pbar = 2.5, slope = 0.2),
      parameter = paste("F̄ =", N_star)
    )
  
  # Test varying Pbar (PPH at inflection)
  pbar_values <- c(1.5, 2.0, 2.5)
  data_pbar <- expand.grid(
    female_density = density_seq,
    Pbar = pbar_values
  ) %>%
    mutate(
      pph = calculate_pph(female_density, N_star = 3, Pbar, slope = 0.2),
      parameter = paste("P̄ =", Pbar)
    )
  
  # Test varying slope
  slope_values <- c(0.1, 0.2, 0.3)
  data_slope <- expand.grid(
    female_density = density_seq,
    slope = slope_values
  ) %>%
    mutate(
      pph = calculate_pph(female_density, N_star = 3, Pbar = 2.5, slope),
      parameter = paste("slope =", slope)
    )
  
  # Create comparison plots
  p1 <- ggplot(data_n_star, aes(x = female_density, y = pph, 
                                color = parameter, group = parameter)) +
    geom_line(linewidth = 1) +
    labs(title = "Effect of Varying F̄", x = "Female Density", y = "PPH",
         color = "Parameter") +
    theme_classic()
  
  p2 <- ggplot(data_pbar, aes(x = female_density, y = pph, 
                              color = parameter, group = parameter)) +
    geom_line(linewidth = 1) +
    labs(title = "Effect of Varying P̄", x = "Female Density", y = "PPH",
         color = "Parameter") +
    theme_classic()
  
  p3 <- ggplot(data_slope, aes(x = female_density, y = pph, 
                               color = parameter, group = parameter)) +
    geom_line(linewidth = 1) +
    labs(title = "Effect of Varying Slope", x = "Female Density", y = "PPH",
         color = "Parameter") +
    theme_classic()
  
  # Combine plots if cowplot is available
  if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(p1, p2, p3, ncol = 3)
  } else {
    print(p1)
    print(p2)
    print(p3)
  }
}

# Run sensitivity test if called interactively
if (interactive() && exists("run_sensitivity") && run_sensitivity) {
  test_pph_sensitivity()
}

################################################################################
# END OF SCRIPT
################################################################################