##############################################################################
# Stitch Together Population Plots - Grid Layout
# Rows: January, April (Low), April (High), September
# Columns: Scenario A, B, C
##############################################################################
rm(list=ls())
gc()

# Load libraries
library(magick)
library(dplyr)
library(stringr)

# Directory containing your PNG files
input_dir <- "Results/Utility_Results"
output_dir <- "Stitched_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

trend <- "Stable"  # Options: "Increase", "Stable", "Decrease"

##############################################################################
# Get all PNG files containing the trend
##############################################################################
all_files <- list.files(input_dir, pattern = "\\.png$", full.names = FALSE)
trend_files <- all_files[str_detect(all_files, trend)]

cat(sprintf("Found %d %s files\n", length(trend_files), trend))

##############################################################################
# Parse file information
##############################################################################
file_info <- data.frame(
  filename = character(),
  month = character(),
  scenario = character(),
  type = character(),
  weather = character(),
  model_var = character(),
  stringsAsFactors = FALSE
)

for (fname in trend_files) {
  
  # Determine plot type
  if (str_detect(fname, "action_")) {
    plot_type <- "action"
  } else if (str_detect(fname, "season_dist_")) {
    plot_type <- "season"
  } else {
    next
  }
  
  # Extract month
  month <- case_when(
    str_detect(fname, "January") ~ "January",
    str_detect(fname, "April") ~ "April",
    str_detect(fname, "September") ~ "September",
    TRUE ~ ""
  )
  
  # Extract scenario - look at the end of filename before .png
  scenario <- case_when(
    str_detect(fname, "_A\\.png$") ~ "A",
    str_detect(fname, "_B\\.png$") ~ "B",
    str_detect(fname, "_C\\.png$") ~ "C",
    TRUE ~ ""
  )
  
  # Extract weather (for April)
  weather <- ""
  if (str_detect(fname, "April")) {
    weather <- case_when(
      str_detect(fname, "warm") ~ "warm",
      str_detect(fname, "cold") ~ "cold",
      TRUE ~ ""
    )
  }
  
  # Extract model variation (for April)
  model_var <- ""
  if (str_detect(fname, "April")) {
    model_var <- case_when(
      str_detect(fname, "Low") ~ "Low",
      str_detect(fname, "High") ~ "High",
      TRUE ~ ""
    )
  }
  
  file_info <- rbind(file_info, data.frame(
    filename = fname,
    month = month,
    scenario = scenario,
    type = plot_type,
    weather = weather,
    model_var = model_var,
    stringsAsFactors = FALSE
  ))
}

# Filter to only keep warm April files
file_info <- file_info %>%
  filter(!(month == "April" & weather == "cold"))

cat(sprintf("Filtered to %d files (warm April only)\n", nrow(file_info)))

# Print summary
cat("\nFile summary:\n")
print(file_info %>% 
        group_by(month, model_var, type) %>% 
        summarise(count = n(), .groups = "drop"))

##############################################################################
# Helper function to create grid
##############################################################################
create_grid <- function(files_df, plot_type_name) {
  
  # Define row structure: 4 rows
  rows_def <- list(
    list(month = "January", label = "January", model_var = ""),
    list(month = "April", label = "April (Low)", model_var = "Low"),
    list(month = "April", label = "April (High)", model_var = "High"),
    list(month = "September", label = "September", model_var = "")
  )
  
  scenarios <- c("A", "B", "C")
  
  # Create empty list for rows
  rows <- list()
  
  for (r in 1:length(rows_def)) {
    row_def <- rows_def[[r]]
    row_images <- list()
    
    for (scen in scenarios) {
      # Filter for this month, scenario, and model_var
      img_file <- files_df %>%
        filter(month == row_def$month, 
               scenario == scen,
               (model_var == row_def$model_var | row_def$model_var == ""))
      
      if (nrow(img_file) > 0) {
        # Load image
        img <- image_read(file.path(input_dir, img_file$filename[1]))
        
        # Get image dimensions
        info <- image_info(img)
        img_width <- info$width
        img_height <- info$height
        
        # Add top margin for column headers (first row only)
        if (r == 1) {
          img <- image_border(img, "white", geometry = "0x60x0x0")
          img <- image_annotate(img, paste0("Scenario ", scen),
                                size = 55,
                                gravity = "north",
                                color = "black",
                                weight = 700,
                                location = "+0+15")
        }
        
        # Add left margin for row headers (first column only)
        if (scen == "A") {
          img <- image_border(img, "white", geometry = "120x0x0x0")
          img <- image_annotate(img, row_def$label,
                                size = 50,
                                gravity = "west",
                                color = "black",
                                weight = 700,
                                degrees = 270,
                                location = "+30+0")
        }
        
        row_images[[length(row_images) + 1]] <- img
      } else {
        # Create blank placeholder if file missing
        blank <- image_blank(width = 750, height = 750, color = "white")
        blank <- image_annotate(blank, "Missing",
                                size = 30,
                                gravity = "center",
                                color = "red")
        row_images[[length(row_images) + 1]] <- blank
        
        cat(sprintf("WARNING: Missing file for %s - Scenario %s - %s\n", 
                    row_def$label, scen, plot_type_name))
      }
    }
    
    # Append row horizontally
    row <- image_append(do.call(c, row_images), stack = FALSE)
    rows[[length(rows) + 1]] <- row
  }
  
  # Stack all rows vertically
  grid <- image_append(do.call(c, rows), stack = TRUE)
  
  return(grid)
}

##############################################################################
# Create grids for action and season plots
##############################################################################

# Action plots
action_files <- file_info %>% filter(type == "action")
if (nrow(action_files) > 0) {
  cat("\nCreating action plot grid...\n")
  cat(sprintf("Using %d action files\n", nrow(action_files)))
  action_grid <- create_grid(action_files, "Action Plots")
  image_write(action_grid,
              path = file.path(output_dir, paste0("action_", trend, "_grid.png")),
              density = 300)
  image_write(action_grid,
              path = file.path(output_dir, paste0("action_", trend, "_grid.pdf")),
              format = "pdf")
  cat("Saved action plot grid\n")
}

# Season plots
season_files <- file_info %>% filter(type == "season")
if (nrow(season_files) > 0) {
  cat("\nCreating season plot grid...\n")
  cat(sprintf("Using %d season files\n", nrow(season_files)))
  season_grid <- create_grid(season_files, "Season Distribution")
  image_write(season_grid,
              path = file.path(output_dir, paste0("season_", trend, "_grid.png")),
              density = 300)
  image_write(season_grid,
              path = file.path(output_dir, paste0("season_", trend, "_grid.pdf")),
              format = "pdf")
  cat("Saved season plot grid\n")
}

cat(sprintf("\nDone! Grid plots saved in: %s\n", output_dir))