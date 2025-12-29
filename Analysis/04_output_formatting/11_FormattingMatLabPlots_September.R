##############################################################################
# Stitch Together September Plots - Compare Population Trends
# Rows: Increase, Stable, Decrease
# Columns: Scenario A, B, C
##############################################################################
rm(list=ls())
gc()

# Load libraries
library(magick)
library(dplyr)
library(stringr)

# Directory containing your PNG files
input_dir <- "Results/ResultsMastLowSept"
output_dir <- "Stitched_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

mast_level <- "Low"  # or "High" - choose which mast level to show

##############################################################################
# Get all PNG files containing September
##############################################################################
all_files <- list.files(input_dir, pattern = "\\.png$", full.names = FALSE)
sept_files <- all_files[str_detect(all_files, "September")]

cat(sprintf("Found %d September files\n", length(sept_files)))

##############################################################################
# Parse file information
##############################################################################
file_info <- data.frame(
  filename = character(),
  trend = character(),
  scenario = character(),
  type = character(),
  mast = character(),
  stringsAsFactors = FALSE
)

for (fname in sept_files) {
  
  # Determine plot type
  if (str_detect(fname, "action_")) {
    plot_type <- "action"
  } else if (str_detect(fname, "season_dist_")) {
    plot_type <- "season"
  } else {
    next
  }
  
  # Extract population trend
  trend <- case_when(
    str_detect(fname, "Increase") ~ "Increase",
    str_detect(fname, "Stable") ~ "Stable",
    str_detect(fname, "Decrease") ~ "Decrease",
    TRUE ~ ""
  )
  
  # Extract mast level (adjust pattern matching if needed)
  mast <- case_when(
    str_detect(tolower(fname), "low") ~ "Low",
    str_detect(tolower(fname), "high") ~ "High",
    TRUE ~ ""
  )
  
  # Extract scenario - look at the end of filename before .png
  scenario <- case_when(
    str_detect(fname, "_A\\.png$") ~ "A",
    str_detect(fname, "_B\\.png$") ~ "B",
    str_detect(fname, "_C\\.png$") ~ "C",
    TRUE ~ ""
  )
  
  file_info <- rbind(file_info, data.frame(
    filename = fname,
    trend = trend,
    scenario = scenario,
    type = plot_type,
    mast = mast,
    stringsAsFactors = FALSE
  ))
}

# Print what we found
cat("\nFiles by trend:\n")
print(table(file_info$trend, file_info$type))

##############################################################################
# Helper function to create grid
##############################################################################
create_grid <- function(files_df, plot_type_name) {
  
  trend_order <- c("Increase", "Stable", "Decrease")
  scenarios <- c("A", "B", "C")
  
  # Create empty list for rows
  rows <- list()
  
  for (pop_trend in trend_order) {
    row_images <- list()
    
    for (scen in scenarios) {
      # Filter for this trend and scenario
      img_file <- files_df %>%
        filter(trend == pop_trend, scenario == scen)
      
      if (nrow(img_file) > 0) {
        # Load image
        img <- image_read(file.path(input_dir, img_file$filename[1]))
        
        # Add top margin for column headers (first row only)
        if (pop_trend == "Increase") {
          img <- image_border(img, "white", geometry = "0x60x0x0")
          img <- image_annotate(img, paste0("Scenario ", scen),
                                size = 45,
                                gravity = "north",
                                color = "black",
                                weight = 700,
                                location = "+0+15")
        }
        
        # Add left margin for row headers (first column only)
        if (scen == "A") {
          img <- image_border(img, "white", geometry = "80x0x0x0")
          trend_label <- pop_trend
          img <- image_annotate(img, trend_label,
                                size = 55,
                                gravity = "west",
                                color = "black",
                                weight = 700,
                                degrees = 270,
                                location = "+25+0")
        }
        
        row_images[[length(row_images) + 1]] <- img
      } else {
        # Create blank placeholder if file missing
        blank <- image_blank(width = 750, height = 750, color = "white")
        blank <- image_annotate(blank, 
                                paste0("Missing: ", pop_trend, " ", scen),
                                size = 30,
                                gravity = "center",
                                color = "red")
        row_images[[length(row_images) + 1]] <- blank
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
              path = file.path(output_dir, paste0("Sept_action_trends_", mast_level, "Mast_grid.png")),
              density = 300)
  image_write(action_grid,
              path = file.path(output_dir, paste0("Sept_action_trends_", mast_level, "Mast_grid.pdf")),
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
              path = file.path(output_dir, paste0("Sept_season_trends_", mast_level, "Mast_grid.png")),
              density = 300)
  image_write(season_grid,
              path = file.path(output_dir, paste0("Sept_season_trends_", mast_level, "Mast_grid.pdf")),
              format = "pdf")
  cat("Saved season plot grid\n")
}

cat(sprintf("\nDone! Grid plots saved in: %s\n", output_dir))
