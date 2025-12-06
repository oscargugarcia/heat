# ===============================================================================
# CREATE TEST DATA FOR UNIT TESTS
# ===============================================================================
#
# Description:
#   This script creates minimal test data for unit testing the climate
#   aggregation pipeline. It generates:
#   - Environmental raster (3x3 grid, 30 daily layers)
#   - Secondary weight raster (9x9 grid, finer resolution)
#   - Two non-overlapping test polygons
#   - Two test points in different cells
#
# WARNING:
#   This script will OVERWRITE existing test data files!
#   Only run this when you intentionally want to recreate the test data.
#
# Last Updated: 2025-11-10 16:01:15 PST
# Author: Jonas Wallstein
#
# Output Structure:
#   unit_tests/data/
#     ├── env_rast/
#     │   └── env_rast.tif (30 layers, 2000-01-01 to 2000-01-30)
#     ├── sec_weight_rast/
#     │   └── sec_weight_rast.tif (1 layer, year 2000)
#     ├── polygons.gpkg (2 test polygons)
#     └── points.gpkg (2 test points)
#
# ===============================================================================

# ---- Safety Check: Get User Confirmation ----------------------------------------------------------

cat("\n")
cat(strrep("=", 80), "\n")
cat("WARNING: CREATE TEST DATA\n")
cat(strrep("=", 80), "\n")
cat("This script will create/overwrite test data files for unit tests.\n")
cat("This should only be run when you intentionally want to recreate the test data.\n")
cat("\nThe following files will be created/overwritten:\n")
cat("  - unit_tests/data/env_rast/env_rast.tif\n")
cat("  - unit_tests/data/sec_weight_rast/sec_weight_rast.tif\n")
cat("  - unit_tests/data/polygons.gpkg\n")
cat("  - unit_tests/data/points.gpkg\n")
cat(strrep("=", 80), "\n")

response <- readline(prompt = "Do you want to continue? (yes/no): ")

if (tolower(trimws(response)) != "yes") {
  stop("Script execution cancelled by user.")
}

cat("\nProceeding with test data creation...\n\n")

# ---- 1. Setup ----------------------------------------------------------

library(terra)
library(sf)

#' Create Simple Test Data
#' 
#' Creates a minimal example with two polygons and environmental/weight rasters
#' for testing the climate aggregation pipeline
create_test_data <- function() {
  
  # ---- 1. Create a simple 3x3 environmental raster ----
  
  # Define extent (simple coordinate system)
  ext <- ext(0, 3, 0, 3)
   
  # Create 100 layers (1999-12-01 to 2000-04-09) with different random values
  env_rast_list <- list()
  start_date <- as.Date("1999-12-01")
  
  for (i in 1:100) {
    set.seed(123 + i)  # Different seed for each layer
    layer_values <- runif(9, min = -10, max = 10)
    current_date <- start_date + (i - 1)
    
    env_rast_layer <- rast(ncols = 3, nrows = 3, ext = ext, 
                           vals = layer_values, 
                           crs = "EPSG:4326")
    names(env_rast_layer) <- as.character(current_date)
    env_rast_list[[i]] <- env_rast_layer
  }
  
  # Combine into a multi-layer raster
  env_rast <- rast(env_rast_list)
  
  
  # ---- 2. Create secondary weight raster (9x9, finer resolution) ----
  
  set.seed(789)
  sec_weight_values <- runif(81, min = 0.1, max = 10)  # Positive values only
  sec_weight_rast <- rast(ncols = 9, nrows = 9, ext = ext, 
                          vals = sec_weight_values, 
                          crs = "EPSG:4326")
  names(sec_weight_rast) <- "2000"
  
  
  # ---- 3. Create two neighboring non-intersecting polygons ----
  
  # Polygon 1: Small square in the lower-left area
  poly1_coords <- matrix(c(
    0.5, 0.5,
    1.2, 0.5,
    1.2, 1.2,
    0.5, 1.2,
    0.5, 0.5
  ), ncol = 2, byrow = TRUE)
  
  poly1 <- st_polygon(list(poly1_coords))
  
  # Polygon 2: Rectangle in the upper-right area
  poly2_coords <- matrix(c(
    1.5, 1.5,
    2.7, 1.5,
    2.7, 2.5,
    1.5, 2.5,
    1.5, 1.5
  ), ncol = 2, byrow = TRUE)
  
  poly2 <- st_polygon(list(poly2_coords))
  
  # Create sf object with both polygons
  polygons <- st_sf(
    poly_id = c("poly_1", "poly_2"),
    name = c("Polygon 1", "Polygon 2"),
    geometry = st_sfc(poly1, poly2, crs = "EPSG:4326")
  )
  
  
  # ---- 4. Create two points in different cells ----
  
  # Point 1: In the lower-left cell (centered)
  point1_coords <- matrix(c(0.5, 0.5), ncol = 2)
  point1 <- st_point(point1_coords)
  
  # Point 2: In the upper-right cell (centered)
  point2_coords <- matrix(c(2.5, 2.5), ncol = 2)
  point2 <- st_point(point2_coords)
  
  # Create sf object with both points
  points <- st_sf(
    point_id = c("point_1", "point_2"),
    name = c("Point 1", "Point 2"),
    geometry = st_sfc(point1, point2, crs = "EPSG:4326")
  )
  
  
  # ---- 5. Return all test data ----
  
  return(list(
    env_rast = env_rast,
    sec_weight_rast = sec_weight_rast,
    polygons = polygons,
    points = points
  ))
}


#' Plot Environmental Raster with Polygons
#' 
#' @param env_rast SpatRaster with environmental data
#' @param polygons sf object with polygon geometries
#' @param layer Layer number to plot (default = 1)
plot_env_rast_polygons <- function(env_rast, polygons, layer = 1) {
  
  # Set up plot area
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 6))
  
  # Plot the environmental raster (first layer)
  plot(env_rast[[layer]], 
       main = paste0("Environmental Raster (", names(env_rast)[layer], ")\nwith Test Polygons"),
       col = hcl.colors(20, "RdYlBu", rev = TRUE),
       legend = TRUE)
  
  # Add grid lines for cells
  grid_x <- seq(ext(env_rast)[1], ext(env_rast)[2], length.out = ncol(env_rast) + 1)
  grid_y <- seq(ext(env_rast)[3], ext(env_rast)[4], length.out = nrow(env_rast) + 1)
  abline(v = grid_x, col = "gray70", lty = 2, lwd = 0.5)
  abline(h = grid_y, col = "gray70", lty = 2, lwd = 0.5)
  
  # Overlay polygons
  plot(st_geometry(polygons), add = TRUE, border = "black", lwd = 2)
  
  # Add polygon labels
  centroids <- st_coordinates(st_centroid(polygons))
  text(centroids[, 1], centroids[, 2], 
       labels = polygons$poly_id, 
       col = "black", font = 2, cex = 0.8)
}


#' Plot Secondary Weight Raster with Polygons
#' 
#' @param sec_weight_rast SpatRaster with secondary weights
#' @param polygons sf object with polygon geometries
plot_sec_weight_polygons <- function(sec_weight_rast, polygons) {
  
  # Set up plot area
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 6))
  
  # Plot the secondary weight raster
  plot(sec_weight_rast, 
       main = paste0("Secondary Weight Raster (", names(sec_weight_rast), ")\nwith Test Polygons"),
       col = hcl.colors(20, "YlOrRd", rev = FALSE),
       legend = TRUE)
  
  # Add grid lines for cells (9x9)
  grid_x <- seq(ext(sec_weight_rast)[1], ext(sec_weight_rast)[2], 
                length.out = ncol(sec_weight_rast) + 1)
  grid_y <- seq(ext(sec_weight_rast)[3], ext(sec_weight_rast)[4], 
                length.out = nrow(sec_weight_rast) + 1)
  abline(v = grid_x, col = "gray70", lty = 2, lwd = 0.3)
  abline(h = grid_y, col = "gray70", lty = 2, lwd = 0.3)
  
  # Overlay polygons
  plot(st_geometry(polygons), add = TRUE, border = "blue", lwd = 2)
  
  # Add polygon labels
  centroids <- st_coordinates(st_centroid(polygons))
  text(centroids[, 1], centroids[, 2], 
       labels = polygons$poly_id, 
       col = "blue", font = 2, cex = 0.8)
}


# ---- Example usage ----

# Create test data
test_data <- create_test_data()

# Extract components
env_rast <- test_data$env_rast
sec_weight_rast <- test_data$sec_weight_rast
polygons <- test_data$polygons
points <- test_data$points

# save the components in the folder unit_tests/data/
# dir.create("unit_tests/test_data", showWarnings = FALSE)
writeRaster(env_rast, filename = "unit_tests/data/env_rast/env_rast.tif", overwrite = TRUE)
writeRaster(sec_weight_rast, filename = "unit_tests/data/sec_weight_rast/sec_weight_rast.tif", overwrite = TRUE)
st_write(polygons, dsn = "unit_tests/data/polygons.gpkg", delete_dsn = TRUE)
st_write(points, dsn = "unit_tests/data/points.gpkg", delete_dsn = TRUE)

message("\n")
message(strrep("=", 80))
message("TEST DATA CREATION COMPLETE")
message(strrep("=", 80))
message("Test data files have been saved:")
message("  - unit_tests/data/env_rast/env_rast.tif (30 daily layers)")
message("  - unit_tests/data/sec_weight_rast/sec_weight_rast.tif")
message("  - unit_tests/data/polygons.gpkg (2 polygons)")
message("  - unit_tests/data/points.gpkg (2 points)")
message(strrep("=", 80))

# ---- Update Timestamp ----------------------------------------------------------

# Record the creation timestamp
timestamp <- Sys.time()
cat("\n")
cat(strrep("=", 80), "\n")
cat("TEST DATA CREATED ON: ", format(timestamp, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat(strrep("=", 80), "\n")

# Update the timestamp in this script file
script_path <- "unit_tests/create_test_data.R"
script_lines <- readLines(script_path)

# Find the line with "Last Updated:" and replace it
timestamp_line_idx <- grep("^# Last Updated:", script_lines)
if (length(timestamp_line_idx) > 0) {
  script_lines[timestamp_line_idx] <- paste0("# Last Updated: ", format(timestamp, "%Y-%m-%d %H:%M:%S %Z"))
  writeLines(script_lines, script_path)
  message("\nTimestamp updated in script file.")
} else {
  message("\nWarning: Could not find timestamp line in script to update.")
}

# # Plot 1: Environmental raster with polygons (layer 1)
# plot_env_rast_polygons(env_rast, polygons, layer = 1)
# 
# # Plot 2: second layer of env_rast
# plot_env_rast_polygons(env_rast, polygons, layer = 2)
# 
# # Plot 3: Secondary weight raster with polygons
# plot_sec_weight_polygons(sec_weight_rast, polygons)

# Plot 4: Environmental raster with points
# plot_env_rast_polygons(env_rast, points, layer = 1)







