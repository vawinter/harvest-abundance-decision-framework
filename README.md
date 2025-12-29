# A Decision Framework for Balancing Hunting Opportunity and Population Abundance in Wild Turkey Management

## Overview

This repository contains code for evaluating optimal harvest season lengths for Pennsylvania wild turkeys using a Markov Decision Process (MDP) framework. The analysis balances two competing management objectives:

1. **Maximizing male abundance** for spring hunting  
2. **Maximizing fall hunting opportunity** (season length)

Key contributions include:
- Integration of **stakeholder-derived utility weights** from hunter preference surveys
- Evaluation of **decision timing** (January, April, September) to quantify the value of information
- Weather-informed **recruitment predictions** for April decision updates

## Analysis Components

### 01_prep_choice_model

Analyzes discrete choice experiment survey data to quantify hunter preferences for season length under different population trends.

**Scripts:**
- `01_visualize_choice_model.R` - Visualize raw preference patterns  
- `02_calculate_utility_weights.R` - Calculate normalized utility weights

**Key Outputs:**
- `norm_weights_popsize_scenario_[A/B/C].csv` - Utility weights for three scenarios reflecting uncertainty in support for closed seasons during population decline

**Methods:** Bias-reduced logistic regression (brglm2 package)

---

### 02_prep_population_parameters

Derives density-dependent population parameters from multiple data sources.

#### Data Preparation Scripts

**00_generate_sampling_points.R**
- Creates 300 random points per WMU region for spatially representative weather sampling

**01_extract_density_from_IPM.R**  
- Extracts female density, poults per brood (PPB), hen-with-brood (HWB) from integrated population model (Winter et al., in review)
- Provides scaling factors to convert BBS indices to absolute densities

#### Parameter Estimation Scripts

**02_BBS_logistic_growth.R**
- Fits logistic growth models to Breeding Bird Survey (BBS) data (1966-2023)
- Estimates carrying capacity (K), intrinsic growth rate (r), inflection points  
- Scales BBS indices to absolute densities using IPM estimates

**Key Parameters:** F̄ (female density at inflection), P̄ (poult density at inflection), slope

#### Weather-Based Recruitment Prediction Scripts

**03_download_weather_data.R**
- Downloads spring weather data (March-June) from Daymet for 3,000 locations

**04_predict_recruitment_from_weather.R**
- Fits GLMMs for hen-with-brood and poults-per-brood  
- Applies temporal weighting centered on nest initiation (day-of-year 100)
- Creates weather covariates for April predictions

**05_fit_weather_recruitment_model.R**
- Fits temperature-recruitment models with cross-validation (leave-one-year-out, leave-one-region-out)
- Effect size analysis with standardized coefficients

**06_generate_april_predictions.R**
- Predicts April recruitment based on spring temperature
- Compares April vs September parameter estimates
- Calculates value of information from waiting

#### Scenario Creation Scripts

**07_create_population_scenarios.R**
- Creates baseline scenarios from observed median parameters
- Generates ±10% adjustments for increase/decrease population trends  
- Applies weather model refinements for April timing
- Produces 9 scenario parameter files (3 decision timings × 3 population trends)

**Decision Timing Scenarios:**

| Timing | Information | Pbar Source | Mast Uncertainty |
|--------|------------|-------------|------------------|
| **January** | Lowest | Baseline (±10%) | High (osig = 2.3) |
| **April** | Medium | Weather-predicted | High (osig = 2.3) |
| **September** | Highest | Observed recruitment | Low (osig = 1.0) |

#### Helper Functions

**99_pph_density_fun.R**
- Density-dependent fecundity function (reverse sigmoidal)
- P(F) = [2η × P̄] / [η + 1 + (η-1) × (F/F̄)^η]

---

### 03_MDP_execution

Implements stochastic dynamic programming to determine optimal season lengths.

**Scripts:**
1. `01_main_decision_model.m` - Core MDP implementation  
2. `02_helper_functions/` - PennTurkeyModel, plotting utilities
3. `03_run_all_scenarios.m` - Batch execution across scenarios
4. `04_calculate_VOI.m` - Value of information analysis

**Model Components:**
- **State variables:** Male (M), Jake (J), Female (F) densities
- **Decision variable:** Season length (0, 1, 2, or 3 weeks)
- **Utility function:** U(L, M+J) = (L + α₁)^uw × (M+J + α₂)^(1-uw)
  - uw = 0.2 (utility weight balancing season length vs. male abundance)
  - α₁, α₂ = hunter preference-derived scaling factors

**Key Outputs:**
- Optimal season length by population state
- Long-run population distributions  
- Expected utility for each decision timing
- Value of information (VOI) calculations

**Methods:** MDPSolve toolbox (Fackler, 2014), MATLAB

---

## Software Requirements

### R (version 4.0+)
- **Packages:** dplyr, tidyr, ggplot2, glmmTMB, brglm2, sf, bbsBayes2, daymetr, patchwork, viridis, caret

### MATLAB (R2020a+)
- **Toolboxes:** MDPSolve, Statistics and Machine Learning Toolbox

---

## Data Requirements

### Input Data (not included in repository)

1. **Hunter preference survey**
   - `2024 Fall Turkey Choice Model Dataset - Veronica Model - arc 5.22.2025 CSV.csv`

2. **Integrated Population Model outputs** (Winter et al., in review)
   - `O_24_abundance_summary.rds` - Density estimates
   - `O_24_ppb_summary.rds` - Poults per brood
   - `O_24_hwb_summary.rds` - Hen-with-brood ratios

3. **BBS data**
   - Downloaded via `bbsBayes2::fetch_bbs_data()`

4. **Pennsylvania WMU boundaries**
   - `PGC_BNDWildlifeManagementUnits2024.shp`

5. **Nest timing data**
   - `20250131_NestAttempts_allbirds.csv`

### Generated Data Files

Intermediate and final outputs saved to `Data/` subdirectories:
- `norm_weights_popsize_scenario_*.csv` - Hunter preference weights
- `region_parameters_for_mdp.csv` - BBS-derived parameters
- `mdp_params_*.csv` - Scenario-specific parameters (9 files)

---

## Running the Analysis

### Complete Workflow

```R
# 0. Prepare data (run in order 00-06)
Rscript 00_data_prep/00_generate_sampling_points.R
Rscript 00_data_prep/01_extract_density_from_IPM.R
Rscript 00_data_prep/02_BBS_logistic_growth.R
Rscript 00_data_prep/03_download_weather_data.R
Rscript 00_data_prep/04_predict_recruitment_from_weather.R
Rscript 00_data_prep/05_fit_weather_recruitment_model.R
Rscript 00_data_prep/06_generate_april_predictions.R

# 1. Hunter preferences
Rscript 01_prep_choice_model/01_visualize_choice_model.R
Rscript 01_prep_choice_model/02_calculate_utility_weights.R

# 2. Population parameters (run in order 00-07)
Rscript 02_prep_population_parameters/07_create_population_scenarios.R

# 3. Decision model (MATLAB)
% In MATLAB:
cd 03_MDP_execution
run('01_main_decision_model.m')
run('04_calculate_VOI.m')
```

---

## Key Results

### Value of Information (VOI)

VOI quantifies utility gain from waiting for better information:
- **VOI(April)** = U(April) - U(January)
- **VOI(September)** = U(September) - U(April)

Results demonstrate utility improvements from delayed decisions, balanced against logistical constraints of regulation changes.

---

## About

Code to implement the analysis described in Winter et al. (20XX) "A Decision Framework for Balancing Hunting Opportunity and Population Abundance in Wild Turkey Management"

**Authors:** Veronica A. Winter, Duane R. Diefenbach, Mary Jo Casalena, Alexander R. Caraynoff, Andrew J. Mowen, Frances E. Buderman

## Contact

Veronica A. Winter - vaw5154@psu.edu