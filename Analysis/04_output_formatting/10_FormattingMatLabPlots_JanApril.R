##############################################################################
# Stitch Together Stable Population Plots - Grid Layout
# Rows: January, April (warm), September
# Columns: Scenario A, B, C
##############################################################################
rm(list=ls())
gc()

# Load libraries
library(magick)
library(dplyr)
library(stringr)

# Directory containing your PNG files
input_dir <- "Results/"
output_dir <- "Stitched_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

trend <- "increase"

##############################################################################
# Get all PNG files containing 'stable'
##############################################################################
all_files <- list.files(input_dir, pattern = "\\.png$", full.names = FALSE)
trend_files <- all_files[str_detect(tolower(all_files), trend)]

cat(sprintf(paste("Found %d", trend, "files\n", sep = " "), length(trend_files)))

##############################################################################
# Parse file information
##############################################################################
file_info <- data.frame(
  filename = character(),
  month = character(),
  scenario = character(),
  type = character(),
  weather = character(),
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
  
  # Extract weather (only warm for April)
  weather <- ""
  if (str_detect(fname, "April")) {
    weather <- case_when(
      str_detect(fname, "warm") ~ "warm",
      str_detect(fname, "cold") ~ "cold",
      TRUE ~ ""
    )
  }
  
  file_info <- rbind(file_info, data.frame(
    filename = fname,
    month = month,
    scenario = scenario,
    type = plot_type,
    weather = weather,
    stringsAsFactors = FALSE
  ))
}

# Filter to only keep warm April files
file_info <- file_info %>%
  filter(!(month == "April" & weather == "cold"))

cat(sprintf("Filtered to %d files (warm April only)\n", nrow(file_info)))

##############################################################################
# Helper function to create grid
##############################################################################
create_grid <- function(files_df, plot_type_name) {
  
  month_order <- c("January", "April", "September")
  scenarios <- c("A", "B", "C")
  
  # Create empty list for rows
  rows <- list()
  
  for (mon in month_order) {
    row_images <- list()
    
    for (scen in scenarios) {
      # Filter for this month and scenario
      img_file <- files_df %>%
        filter(month == mon, scenario == scen)
      
      if (nrow(img_file) > 0) {
        # Load image
        img <- image_read(file.path(input_dir, img_file$filename[1]))
        
        # Get image dimensions
        info <- image_info(img)
        img_width <- info$width
        img_height <- info$height
        
        # Add top margin for column headers (first row only)
        if (mon == "January") {
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
          img <- image_border(img, "white", geometry = "80x0x0x0")
          month_label <- mon
          img <- image_annotate(img, month_label,
                                size = 55,
                                gravity = "west",
                                color = "black",
                                weight = 700,
                                degrees = 270,  # Changed from 90 to 270
                                location = "+25+0")
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
      }
    }
    
    # Append row horizontally
    row <- image_append(do.call(c, row_images), stack = FALSE)
    rows[[length(rows) + 1]] <- row
  }
  
  # Stack all rows vertically
  grid <- image_append(do.call(c, rows), stack = TRUE)
  
  # Remove overall title - comment out these lines:
  # grid <- image_border(grid, "white", "50x120x50x50")
  # grid <- image_annotate(grid,
  #                        paste0(plot_type_name, ": Stable Population"),
  #                        size = 60,
  #                        gravity = "north",
  #                        color = "black",
  #                        weight = 700,
  #                        location = "+0+30")
  
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
              path = file.path(output_dir, paste0("action_",trend,"_grid.png")),
              density = 300)
  image_write(action_grid,
              path = file.path(output_dir, paste0("action_",trend,"_grid.pdf")),
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
              path = file.path(output_dir, paste0("season_",trend,"_grid.png")),
              density = 300)
  image_write(season_grid,
              path = file.path(output_dir, paste0("season_",trend,"_grid.pdf")),
              format = "pdf")
  cat("Saved season plot grid\n")
}

cat(sprintf("\nDone! Grid plots saved in: %s\n", output_dir))
