################################################################################
# Fecundity Variation Analysis
# Coefficient of Variation (CV) to Standard Deviation (SD) Conversion
################################################################################
# Purpose: Explore how changes in CV affect fecundity variation for wild turkey
#          population models across different decision timing scenarios
# Author: Veronica A. Winter
# Date: Janaury 2026
################################################################################

library(ggplot2)

# Parameters ------------------------------------------------------------------
mean_fecundity <- 1  # Burr3 distribution has mean = 1
cvw_january <- 1.116 / 4.96  # Current January CV value

# Calculate derived values
sd_january <- cvw_january * mean_fecundity  # SD = CV * mean
se_january <- sd_january / sqrt(1)  # SE = SD / sqrt(n), n=1

# Display current parameters --------------------------------------------------
cat("Current January Parameters\n")
cat("==========================\n")
cat(sprintf("  Mean fecundity: %.4f\n", mean_fecundity))
cat(sprintf("  CV:             %.4f\n", cvw_january))
cat(sprintf("  SD:             %.4f\n", sd_january))
cat(sprintf("  SE:             %.4f\n\n", se_january))

# Explore range of CV values -------------------------------------------------
cv_values <- seq(0, 0.5, by = 0.025)

cv_results <- data.frame(
  cv = cv_values,
  sd = cv_values * mean_fecundity,
  se = (cv_values * mean_fecundity) / sqrt(1)
)

cat("\nCV to SD Conversion Table\n")
cat("=========================\n")
print(cv_results, row.names = FALSE)

# Scenario comparison ---------------------------------------------------------
comparison <- data.frame(
  scenario = c("September (known)", 
               "April (low var)",
               "April (high var)", 
               "January (baseline)",
               "Doubled variation"),
  cv = c(0, 
         cvw_january * 0.5, 
         cvw_january * 2.0,
         cvw_january, 
         cvw_january * 2),
  multiplier = c(0, 0.5, 2.0, 1.0, 2.0)
)

comparison$sd <- comparison$cv * mean_fecundity
comparison$se <- comparison$sd / sqrt(1)

cat("\n\nScenario Comparison\n")
cat("===================\n")
print(comparison, row.names = FALSE)

# Visualization ---------------------------------------------------------------
p <- ggplot(cv_results, aes(x = cv, y = sd)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_vline(xintercept = cvw_january, 
             linetype = "dashed", 
             color = "red", 
             linewidth = 0.8) +
  annotate("text", 
           x = cvw_january + 0.03, 
           y = max(cv_results$sd) * 0.85, 
           label = "January\n(baseline)", 
           color = "red",
           size = 4) +
  # Add scenario markers
  geom_point(data = comparison, 
             aes(x = cv, y = sd), 
             color = "darkgreen", 
             size = 3, 
             shape = 17) +
  labs(x = "Coefficient of Variation (CV)",
       y = "Standard Deviation (SD)",
       title = "Fecundity Variation: CV to SD Conversion",
       subtitle = "Mean = 1 (Burr3 distribution)",
       caption = "Triangles indicate decision timing scenarios") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 12),
        plot.caption = element_text(hjust = 0, face = "italic"))

print(p)

# Save plot (optional) --------------------------------------------------------
# ggsave("fecundity_cv_analysis.png", p, width = 8, height = 6, dpi = 300)

# Summary statistics ----------------------------------------------------------
cat("\n\nSummary Statistics\n")
cat("==================\n")
cat(sprintf("CV range explored: %.4f to %.4f\n", 
            min(cv_values), max(cv_values)))
cat(sprintf("SD range: %.4f to %.4f\n", 
            min(cv_results$sd), max(cv_results$sd)))
cat(sprintf("\nJanuary baseline:\n"))
cat(sprintf("  CV = %.4f corresponds to SD = %.4f\n", cvw_january, sd_january))
cat(sprintf("\nApril scenarios:\n"))
cat(sprintf("  Low variation (50%%):  CV = %.4f, SD = %.4f\n", 
            cvw_january * 0.5, cvw_january * 0.5))
cat(sprintf("  High variation (200%%): CV = %.4f, SD = %.4f\n", 
            cvw_january * 2.0, cvw_january * 2.0))

################################################################################
# Notes:
# - CV (Coefficient of Variation) = SD / mean
# - For mean = 1, CV = SD numerically, but they represent different concepts
# - September: cvw = 0 (known recruitment, no variation)
# - January: cvw = baseline (1.116 / 4.96 ≈ 0.225)
# - April Low: cvw = 0.5 × baseline (reduced uncertainty with weather model)
# - April High: cvw = 2.0 × baseline (increased uncertainty)
################################################################################