Analysis/
├── README.md                          # Overview of analysis workflow
├── 00_data/                          # Raw and processed data
│   ├── raw/                          # Original survey data, BBS data
│   └── processed/                    # Cleaned data ready for analysis
│
├── 01_prep_choice_model/             # Part 1 - R scripts
│   ├── 01_choice_model_analysis.R
│   ├── 02_utility_weights.R
│   └── README.md                     # What each script does
│
├── 02_prep_population_parameters/    # Density-dependent parameters
│   ├── 01_BBS_analysis.R             # BBS logistic growth model
│   ├── 02_inflection_point_calc.R
│   └── README.md
│
├── 03_MDP_execution/                 # Part 2 - MATLAB scripts
│   ├── 01_main_decision_model.m      # Your Monday_Testing script
│   ├── 02_helper_functions/          # Separate folder for functions
│   │   ├── PennTurkeyModel.m
│   │   ├── actionplot.m
│   │   ├── lrharvplot_fixed.m
│   │   └── getUtilityWeight.m
│   ├── 03_sensitivity_analysis.m     # If you have one
│   └── README.md
│
├── 04_output_formatting/             # Part 3 - Tables and figures
│   ├── 01_create_tables.R
│   ├── 02_create_figures.R
│   └── README.md
│
├── 05_results/                       # Generated outputs
│   ├── figures/
│   ├── tables/
│   └── model_outputs/
│
└── functions/                        # Shared functions (if any)
    └── utility_functions.R