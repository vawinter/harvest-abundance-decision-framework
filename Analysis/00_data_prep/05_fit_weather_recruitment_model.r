################################################################################
# Fit Weather-Recruitment Models with Cross-Validation
################################################################################
#
# Purpose: 
#   1. Fit GLMMs relating spring temperature to recruitment metrics (HWB, PPB, PPH)
#   2. Perform rigorous cross-validation (leave-one-year-out, leave-one-region-out)
#   3. Calculate effect sizes and biological significance
#   4. Generate predictions for April decision model
#
# Method:
#   - GLMMs with temperature as predictor, WMU random effects
#   - Cross-validation to assess predictive performance
#   - Effect size analysis (standardized coefficients, Cohen's f²)
#
# Input:  
#   april_scaled_weather2.rds (from 04_predict_recruitment_from_weather.R)
#   hwb_df_aug31.rds (hen-with-brood predictions)
#   ph_df_aug31.rds (poults per brood predictions)
#   pph_df_aug31.rds (poults per hen predictions)
#
# Output:
#   Model summaries and predictions
#   Cross-validation results and plots
#   Effect size analyses
#
# Author: Veronica A. Winter
# Date: May 2025
################################################################################

rm(list = ls())
gc()

set.seed(33456)

library(ggplot2)
library(dplyr)
library(glmmTMB)
library(car)
library(caret)
library(patchwork)
library(viridis)

################################################################################
# HELPER FUNCTIONS
################################################################################

#------------------------------------------------------------------------------
calc_performance <- function(actual, predicted) {
  #' Calculate Model Performance Metrics
  #'
  #' @param actual Numeric vector of observed values
  #' @param predicted Numeric vector of predicted values
  #' @return List with RMSE, MAE, and R-squared
  
  # Remove NA values
  valid_indices <- !is.na(actual) & !is.na(predicted)
  actual <- actual[valid_indices]
  predicted <- predicted[valid_indices]
  
  # Calculate residuals
  residuals <- actual - predicted
  
  # Root Mean Square Error
  rmse <- sqrt(mean(residuals^2))
  
  # Mean Absolute Error
  mae <- mean(abs(residuals))
  
  # R-squared (proportion of variance explained)
  ss_total <- sum((actual - mean(actual))^2)
  ss_residual <- sum(residuals^2)
  r_squared <- if (ss_total > 0) 1 - (ss_residual / ss_total) else 0
  
  return(list(RMSE = rmse, MAE = mae, R_squared = r_squared))
}

#------------------------------------------------------------------------------
loyo_cv <- function(data, formula, year_var = "Year") {
  #' Leave-One-Year-Out Cross-Validation
  #'
  #' Fits model on all years except one, predicts the held-out year
  #' Tests temporal generalizability
  #'
  #' @param data Data frame with response, predictors, and year variable
  #' @param formula Model formula
  #' @param year_var Name of year column (default "Year")
  #' @return List with performance metrics and predictions
  
  years <- unique(data[[year_var]])
  n_years <- length(years)
  
  # Initialize storage
  loyo_rmse <- numeric(n_years)
  loyo_mae <- numeric(n_years)
  loyo_r2 <- numeric(n_years)
  loyo_year <- years
  
  all_actuals <- c()
  all_predicted <- c()
  all_years <- c()
  
  cat("\nPerforming leave-one-year-out cross-validation...\n")
  
  # Loop through each year
  for (i in 1:n_years) {
    test_year <- years[i]
    cat("  Testing year", test_year, "\n")
    
    # Split data
    test_indices <- which(data[[year_var]] == test_year)
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]
    
    # Fit model on training data
    model <- glmmTMB(formula, data = train_data, 
                     family = gaussian(link = "identity"))
    
    # Predict on test year
    predictions <- predict(model, newdata = test_data, type = "response", 
                           allow.new.levels = TRUE)
    
    # Extract actual values
    response_var <- as.character(formula[[2]])
    actuals <- test_data[[response_var]]
    
    # Store predictions
    all_actuals <- c(all_actuals, actuals)
    all_predicted <- c(all_predicted, predictions)
    all_years <- c(all_years, rep(test_year, length(actuals)))
    
    # Calculate metrics for this fold
    performance <- calc_performance(actuals, predictions)
    loyo_rmse[i] <- performance$RMSE
    loyo_mae[i] <- performance$MAE
    loyo_r2[i] <- performance$R_squared
  }
  
  # Results by year
  results_by_year <- data.frame(
    Year = loyo_year,
    RMSE = loyo_rmse,
    MAE = loyo_mae,
    R_squared = loyo_r2
  )
  
  # All predictions
  predictions_df <- data.frame(
    Actual = all_actuals,
    Predicted = all_predicted,
    Year = all_years
  )
  
  return(list(
    avg_RMSE = mean(loyo_rmse, na.rm = TRUE),
    avg_MAE = mean(loyo_mae, na.rm = TRUE),
    avg_R2 = mean(loyo_r2, na.rm = TRUE),
    by_year = results_by_year,
    predictions = predictions_df
  ))
}

#------------------------------------------------------------------------------
loro_cv <- function(data, formula, region_var = "MU") {
  #' Leave-One-Region-Out Cross-Validation
  #'
  #' Fits model on all regions except one, predicts the held-out region
  #' Tests spatial generalizability
  #'
  #' @param data Data frame with response, predictors, and region variable
  #' @param formula Model formula
  #' @param region_var Name of region column (default "MU")
  #' @return List with performance metrics and predictions
  
  regions <- unique(data[[region_var]])
  n_regions <- length(regions)
  
  # Initialize storage
  loro_rmse <- numeric(n_regions)
  loro_mae <- numeric(n_regions)
  loro_r2 <- numeric(n_regions)
  loro_region <- regions
  
  all_actuals <- c()
  all_predicted <- c()
  all_regions <- c()
  
  cat("\nPerforming leave-one-region-out cross-validation...\n")
  
  # Loop through each region
  for (i in 1:n_regions) {
    test_region <- regions[i]
    cat("  Testing region", test_region, "\n")
    
    # Split data
    test_indices <- which(data[[region_var]] == test_region)
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]
    
    # Fit model on training data
    model <- glmmTMB(formula, data = train_data, 
                     family = gaussian(link = "identity"))
    
    # Predict on test region
    predictions <- predict(model, newdata = test_data, type = "response", 
                           allow.new.levels = TRUE)
    
    # Extract actual values
    response_var <- as.character(formula[[2]])
    actuals <- test_data[[response_var]]
    
    # Store predictions
    all_actuals <- c(all_actuals, actuals)
    all_predicted <- c(all_predicted, predictions)
    all_regions <- c(all_regions, rep(test_region, length(actuals)))
    
    # Calculate metrics for this fold
    performance <- calc_performance(actuals, predictions)
    loro_rmse[i] <- performance$RMSE
    loro_mae[i] <- performance$MAE
    loro_r2[i] <- performance$R_squared
  }
  
  # Results by region
  results_by_region <- data.frame(
    Region = loro_region,
    RMSE = loro_rmse,
    MAE = loro_mae,
    R_squared = loro_r2
  )
  
  # All predictions
  predictions_df <- data.frame(
    Actual = all_actuals,
    Predicted = all_predicted,
    Region = all_regions
  )
  
  return(list(
    avg_RMSE = mean(loro_rmse, na.rm = TRUE),
    avg_MAE = mean(loro_mae, na.rm = TRUE),
    avg_R2 = mean(loro_r2, na.rm = TRUE),
    by_region = results_by_region,
    predictions = predictions_df
  ))
}

#------------------------------------------------------------------------------
plot_cv_results <- function(cv_results, title, color_var) {
  #' Plot Cross-Validation Results
  #'
  #' Creates scatter plot of actual vs predicted values
  #'
  #' @param cv_results Output from loyo_cv or loro_cv
  #' @param title Plot title
  #' @param color_var Variable to color points by ("Year" or "Region")
  
  p <- ggplot(cv_results$predictions, 
              aes(x = Actual, y = Predicted, color = as.factor(get(color_var)))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = title,
      x = "Actual Values",
      y = "Predicted Values",
      color = color_var
    ) +
    theme_classic() +
    scale_color_viridis_d() +
    coord_equal(ratio = 1) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    )
  
  return(p)
}

#------------------------------------------------------------------------------
extract_effect_sizes <- function(model, model_name) {
  #' Extract Effect Sizes from glmmTMB Model
  #'
  #' Gets standardized coefficient, SE, CI, and p-value
  #'
  #' @param model Fitted glmmTMB model
  #' @param model_name Name for output table
  #' @return Data frame with effect size information
  
  # Get fixed effects summary
  fixed_effects <- summary(model)$coefficients$cond
  
  # Extract coefficient and statistics
  coef_val <- fixed_effects[1, "Estimate"]
  se_val <- fixed_effects[1, "Std. Error"]
  ci_lower <- coef_val - 1.96 * se_val
  ci_upper <- coef_val + 1.96 * se_val
  p_val <- fixed_effects[1, "Pr(>|z|)"]
  
  # Create output
  effect_df <- data.frame(
    Model = model_name,
    Coefficient = coef_val,
    SE = se_val,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    P_value = p_val,
    Significant = p_val < 0.05
  )
  
  return(effect_df)
}

################################################################################
# 1. LOAD DATA
################################################################################

cat("\n========== LOADING DATA ==========\n")

# Specify which weather month to use
month <- "april"

# Load weather covariates
scaled_weather <- readRDS(paste0("Data/Rec_data/", month, "_scaled_weather2.rds"))

# Load recruitment metrics
hwb_df <- readRDS("Data/Rec_data/hwb_df_aug31.rds")
ppb_df <- readRDS("Data/Rec_data/ph_df_aug31.rds")
pph_df <- readRDS("Data/Rec_data/pph_df_aug31.rds")

cat("✓ Data loaded successfully\n")

################################################################################
# 2. HEN-WITH-BROOD (HWB) MODEL
################################################################################

cat("\n========== HWB MODEL ==========\n")

# Join recruitment and weather data
hwb <- hwb_df %>%
  left_join(scaled_weather, by = c("Year", "MU"))

# Define model formula
# Predicting HWB from standardized temperature with WMU random effects
hwb_formula <- HWB ~ 0 + scale_avg_temperature + (1|MU)

# Fit model
m.refined.hwb <- glmmTMB(hwb_formula, data = hwb, 
                         family = gaussian(link = "identity"))

cat("Model summary:\n")
print(summary(m.refined.hwb))

# Cross-validation - leave one year out
hwb_loyo <- loyo_cv(hwb, hwb_formula)
cat("\nLeave-One-Year-Out Cross-Validation:\n")
cat("  Average RMSE:", round(hwb_loyo$avg_RMSE, 3), "\n")
cat("  Average MAE:", round(hwb_loyo$avg_MAE, 3), "\n")
cat("  Average R²:", round(hwb_loyo$avg_R2, 3), "\n")

# Create validation plot
hwb_loyo_plot <- plot_cv_results(hwb_loyo, 
                                 "HWB - Leave-One-Year-Out Validation", 
                                 "Year")

################################################################################
# 3. POULTS-PER-BROOD (PPB) MODEL
################################################################################

cat("\n========== PPB MODEL ==========\n")

# Join data
ppb <- ppb_df %>%
  left_join(scaled_weather, by = c("Year", "MU"))

# Define model formula
ppb_formula <- PHratio ~ 0 + scale_avg_temperature + (1|MU)

# Fit model
m.refined.ppb <- glmmTMB(ppb_formula, data = ppb, 
                         family = gaussian(link = "identity"))

cat("Model summary:\n")
print(summary(m.refined.ppb))

# Cross-validation
ppb_loyo <- loyo_cv(ppb, ppb_formula)
cat("\nLeave-One-Year-Out Cross-Validation:\n")
cat("  Average RMSE:", round(ppb_loyo$avg_RMSE, 3), "\n")
cat("  Average MAE:", round(ppb_loyo$avg_MAE, 3), "\n")
cat("  Average R²:", round(ppb_loyo$avg_R2, 3), "\n")

# Create validation plot
ppb_loyo_plot <- plot_cv_results(ppb_loyo, 
                                 "PPB - Leave-One-Year-Out Validation", 
                                 "Year")

################################################################################
# 4. POULTS-PER-HEN (PPH) MODEL
################################################################################

cat("\n========== PPH MODEL ==========\n")

# Join data
pph <- pph_df %>%
  left_join(scaled_weather, by = c("Year", "MU"))

# Define model formula
pph_formula <- pph ~ 0 + scale_avg_temperature + (1|MU)

# Fit model
m.refined.pph <- glmmTMB(pph_formula, data = pph, 
                         family = gaussian(link = "identity"))

cat("Model summary:\n")
print(summary(m.refined.pph))

# Cross-validation
pph_loyo <- loyo_cv(pph, pph_formula)
cat("\nLeave-One-Year-Out Cross-Validation:\n")
cat("  Average RMSE:", round(pph_loyo$avg_RMSE, 3), "\n")
cat("  Average MAE:", round(pph_loyo$avg_MAE, 3), "\n")
cat("  Average R²:", round(pph_loyo$avg_R2, 3), "\n")

# Create validation plot
pph_loyo_plot <- plot_cv_results(pph_loyo, 
                                 "PPH - Leave-One-Year-Out Validation", 
                                 "Year")

################################################################################
# 5. COMBINED VALIDATION PLOTS
################################################################################

cat("\n========== CREATING VALIDATION VISUALIZATIONS ==========\n")

# Combine year validation plots
year_validation_plots <- (hwb_loyo_plot + ppb_loyo_plot + pph_loyo_plot) + 
  plot_layout(ncol = 3) +
  plot_annotation(title = "Leave-One-Year-Out Cross-Validation Results")

# Save
ggsave("Dataviz/year_validation_plots.png", year_validation_plots, 
       width = 12, height = 5, dpi = 300, bg = "white")

cat("✓ Validation plots saved\n")

################################################################################
# 6. GENERATE MODEL PREDICTIONS
################################################################################

cat("\n========== GENERATING PREDICTIONS ==========\n")

# HWB predictions with confidence intervals
hwb_preds <- predict(m.refined.hwb, newdata = hwb, type = "response", 
                     se.fit = TRUE, re.form = NULL)
hwb$HWB_pred <- hwb_preds$fit
hwb$HWB_se <- hwb_preds$se.fit
hwb$HWB_lower <- hwb$HWB_pred - 1.96 * hwb$HWB_se
hwb$HWB_upper <- hwb$HWB_pred + 1.96 * hwb$HWB_se

# PPB predictions with confidence intervals
ppb_preds <- predict(m.refined.ppb, newdata = ppb, type = "response", 
                     se.fit = TRUE, re.form = NULL)
ppb$PPB_pred <- ppb_preds$fit
ppb$PPB_se <- ppb_preds$se.fit
ppb$PPB_lower <- ppb$PPB_pred - 1.96 * ppb$PPB_se
ppb$PPB_upper <- ppb$PPB_pred + 1.96 * ppb$PPB_se

# PPH predictions with confidence intervals
pph_preds <- predict(m.refined.pph, newdata = pph, type = "response", 
                     se.fit = TRUE, re.form = NULL)
pph$PPH_pred <- pph_preds$fit
pph$PPH_se <- pph_preds$se.fit
pph$PPH_lower <- pph$PPH_pred - 1.96 * pph$PPH_se
pph$PPH_upper <- pph$PPH_pred + 1.96 * pph$PPH_se

cat("✓ Predictions generated\n")

################################################################################
# 7. VISUALIZE TEMPERATURE-RECRUITMENT RELATIONSHIPS
################################################################################

cat("\n========== CREATING RELATIONSHIP PLOTS ==========\n")

# HWB vs temperature
hwb_plot <- ggplot(hwb, aes(x = scale_avg_temperature, y = HWB_pred, 
                            color = as.factor(MU))) +
  geom_point(aes(y = HWB), size = 2, alpha = 0.3) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = HWB_lower, ymax = HWB_upper, fill = as.factor(MU)), 
              alpha = 0.1, color = NA) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = "none") +
  labs(
    x = " ", 
    y = "Predicted HWB",
    color = "Region",
    title = paste0("HWB Model (LOYO R² = ", round(hwb_loyo$avg_R2, 2), ")")
  ) +
  theme_classic() +
  theme(legend.position = "none")

# PPB vs temperature
ppb_plot <- ggplot(ppb, aes(x = scale_avg_temperature, y = PPB_pred, 
                            color = as.factor(MU))) +
  geom_point(aes(y = PHratio), size = 2, alpha = 0.3) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PPB_lower, ymax = PPB_upper, fill = as.factor(MU)), 
              alpha = 0.1, color = NA) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = "none") +
  labs(
    x = " ", 
    y = "Predicted PPB",
    color = "Region",
    title = paste0("PPB Model (LOYO R² = ", round(ppb_loyo$avg_R2, 2), ")")
  ) +
  theme_classic() +  
  theme(legend.position = "none")

# PPH vs temperature
pph_plot <- ggplot(pph, aes(x = scale_avg_temperature, y = PPH_pred, 
                            color = as.factor(MU))) +
  geom_point(aes(y = pph), size = 2, alpha = 0.3) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PPH_lower, ymax = PPH_upper, fill = as.factor(MU)), 
              alpha = 0.1, color = NA) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = "none") +
  labs(
    x = "Standardized Average Temperature", 
    y = "Predicted PPH",
    color = "Region",
    title = paste0("PPH Model (LOYO R² = ", round(pph_loyo$avg_R2, 2), ")")
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- hwb_plot / ppb_plot / pph_plot + 
  plot_layout(heights = c(1, 1, 1.3)) +
  plot_annotation(title = "Temperature-Recruitment Relationships")

# Save
ggsave("Dataviz/temperature_recruitment_relationships.png", combined_plot, 
       width = 8, height = 10, dpi = 300, bg = "white")

cat("✓ Relationship plots saved\n")

################################################################################
# 8. EFFECT SIZE ANALYSIS
################################################################################

cat("\n========== EFFECT SIZE ANALYSIS ==========\n")

# Extract effect sizes for all models
hwb_effects <- extract_effect_sizes(m.refined.hwb, "HWB")
ppb_effects <- extract_effect_sizes(m.refined.ppb, "PPB")
pph_effects <- extract_effect_sizes(m.refined.pph, "PPH")

# Combine
all_effects <- rbind(hwb_effects, ppb_effects, pph_effects)

# Add effect size interpretation
effect_interpretation <- all_effects %>%
  mutate(
    Abs_Effect = abs(Coefficient),
    Effect_Size_Category = case_when(
      Abs_Effect < 0.1 ~ "Negligible",
      Abs_Effect < 0.3 ~ "Small", 
      Abs_Effect < 0.5 ~ "Medium",
      TRUE ~ "Large"
    ),
    Direction = ifelse(Coefficient > 0, "Positive", "Negative"),
    Interpretation = paste0(Direction, " ", Effect_Size_Category, " effect"),
    CI_text = paste0("(", round(CI_lower, 3), ", ", round(CI_upper, 3), ")")
  )

cat("\nEffect size summary:\n")
print(effect_interpretation)

# Create effect size plot
effect_size_plot <- ggplot(effect_interpretation, 
                           aes(x = Model, y = Coefficient)) +
  geom_hline(yintercept = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5), 
             linetype = "dotted", alpha = 0.3, color = "gray70") +
  geom_point(aes(color = Effect_Size_Category, shape = Significant), size = 5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.15, size = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.8) +
  scale_color_viridis_d(name = "Effect Size") +
  scale_shape_manual(
    values = c("FALSE" = 1, "TRUE" = 16), 
    name = "Significant\n(p < 0.05)",
    labels = c("No", "Yes")
  ) +
  labs(
    title = "Temperature Effects on Turkey Recruitment",
    subtitle = "Standardized coefficients showing biological significance",
    x = "Recruitment Metric",
    y = "Standardized Effect Size\n(Change per 1 SD increase in temperature)",
    caption = "Error bars show 95% CI\nReference lines at ±0.1, ±0.3, ±0.5"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 9, color = "gray50"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  coord_flip()

# Save
ggsave("Dataviz/temperature_effect_sizes.png", effect_size_plot, 
       width = 10, height = 6, dpi = 300, bg = "white")

cat("✓ Effect size plots saved\n")

################################################################################
# 9. SAVE RESULTS
################################################################################

cat("\n========== SAVING RESULTS ==========\n")

# Validation summary table
validation_summary <- data.frame(
  Metric = c("HWB", "PPB", "PPH"),
  LOYO_RMSE = c(hwb_loyo$avg_RMSE, ppb_loyo$avg_RMSE, pph_loyo$avg_RMSE),
  LOYO_MAE = c(hwb_loyo$avg_MAE, ppb_loyo$avg_MAE, pph_loyo$avg_MAE),
  LOYO_R2 = c(hwb_loyo$avg_R2, ppb_loyo$avg_R2, pph_loyo$avg_R2)
)

write.csv(validation_summary, "Results/validation_summary.csv", row.names = FALSE)

# Effect size table
write.csv(effect_interpretation, "Results/effect_size_summary.csv", row.names = FALSE)

# Model summaries
sink("Results/model_summaries.txt")
cat("========== HWB MODEL ==========\n")
print(summary(m.refined.hwb))
cat("\n\n========== PPB MODEL ==========\n")
print(summary(m.refined.ppb))
cat("\n\n========== PPH MODEL ==========\n")
print(summary(m.refined.pph))
sink()

cat("✓ Results saved\n")

################################################################################
# 10. SUMMARY
################################################################################

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Models fitted:\n")
cat("  1. HWB ~ temperature (LOYO R² =", round(hwb_loyo$avg_R2, 3), ")\n")
cat("  2. PPB ~ temperature (LOYO R² =", round(ppb_loyo$avg_R2, 3), ")\n")
cat("  3. PPH ~ temperature (LOYO R² =", round(pph_loyo$avg_R2, 3), ")\n\n")

cat("Key findings:\n")
for (i in 1:nrow(effect_interpretation)) {
  cat("  ", effect_interpretation$Model[i], ": ", 
      effect_interpretation$Interpretation[i], "\n")
}

cat("\nOutputs saved:\n")
cat("  - Dataviz/temperature_recruitment_relationships.png\n")
cat("  - Dataviz/temperature_effect_sizes.png\n")
cat("  - Dataviz/year_validation_plots.png\n")
cat("  - Results/validation_summary.csv\n")
cat("  - Results/effect_size_summary.csv\n")
cat("  - Results/model_summaries.txt\n")

cat("\nNext steps:\n")
cat("  - Use predictions to update Pbar in April decision model\n")
cat("  - Compare April vs September decision timing in MDP\n")
cat("=======================================\n\n")

################################################################################
# END OF SCRIPT
################################################################################