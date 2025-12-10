# ===============================================================================
# UNIT TEST FOR ALL TEMPORAL RESOLUTIONS
# ===============================================================================
#
# Description:
#   This script tests that r2e2() correctly handles all input temporal resolutions
#   from minutely to yearly, with and without secondary weight rasters. It verifies:
#   - Temporal resolution detection works correctly
#   - Aggregation to daily/monthly/yearly outputs works correctly
#   - Secondary weight rasters at different temporal resolutions work
#
# Test Strategy:
#   - Create synthetic raster data at different temporal resolutions
#   - Test r2e2() pipeline with each resolution
#   - Verify output structure and temporal aggregation
#   - Test secondary weight rasters at various resolutions
#
# Author: Jonas Wallstein
# Last Updated: 2025-12-09
#
# ===============================================================================

library(testthat)

# Load necessary libraries
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(data.table)
  library(dplyr)
  library(lubridate)
})

# Load the heat package
library(heat)

# ---- Helper Functions ----------------------------------------------------------

should_print_banners <- function() {
  nzchar(Sys.getenv("TESTTHAT_VERBOSE")) || isTRUE(getOption("testthat.verbose", FALSE))
}

# Create test geometry
create_test_geometry <- function() {
  poly_coords <- matrix(c(
    0.5, 0.5,
    2.5, 0.5,
    2.5, 2.5,
    0.5, 2.5,
    0.5, 0.5
  ), ncol = 2, byrow = TRUE)
  
  poly <- st_polygon(list(poly_coords))
  
  geometry <- st_sf(
    geom_id = "test_poly",
    geometry = st_sfc(poly, crs = "EPSG:4326")
  )
  
  return(geometry)
}

# Create raster at specified temporal resolution
create_temporal_raster <- function(resolution = "daily", n_periods = 5, seed = 123) {
  ext_obj <- ext(0, 3, 0, 3)
  
  # Generate dates based on resolution
  start_datetime <- as.POSIXct("2000-01-01 00:00", tz = "UTC")
  
  if (resolution == "minutely") {
    datetimes <- seq(start_datetime, by = "min", length.out = n_periods)
    date_format <- "%Y-%m-%d %H:%M"
  } else if (resolution == "hourly") {
    datetimes <- seq(start_datetime, by = "hour", length.out = n_periods)
    date_format <- "%Y-%m-%d %H:%M"
  } else if (resolution == "daily") {
    datetimes <- seq(as.Date("2000-01-01"), by = "day", length.out = n_periods)
    date_format <- "%Y-%m-%d"
  } else if (resolution == "monthly") {
    # For monthly, use YYYY-MM format
    datetimes <- seq(as.Date("2000-01-01"), by = "month", length.out = n_periods)
    date_format <- "%Y-%m"
  } else if (resolution == "yearly") {
    # For yearly, use YYYY format
    datetimes <- seq(as.Date("2000-01-01"), by = "year", length.out = n_periods)
    date_format <- "%Y"
  } else {
    stop("Unsupported resolution: ", resolution)
  }
  
  # Create raster layers
  rast_list <- list()
  for (i in 1:n_periods) {
    set.seed(seed + i)
    values <- runif(9, min = 10, max = 30)
    rast_layer <- rast(ncols = 3, nrows = 3, ext = ext_obj, 
                       vals = values, crs = "EPSG:4326")
    names(rast_layer) <- format(datetimes[i], date_format)
    rast_list[[i]] <- rast_layer
  }
  
  return(rast(rast_list))
}


# ===============================================================================
# TEST 1: Hourly input temporal resolution
# ===============================================================================

test_that("r2e2 handles hourly input data correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 1: Hourly input resolution ===\n")
  }
  
  # Create hourly raster (48 hours = 2 days)
  env_rast <- create_temporal_raster("hourly", n_periods = 48)
  geometry <- create_test_geometry()
  
  # Run r2e2 with hourly data, aggregate to daily
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "none",
    daily_agg_fun = "mean",
    out_temp_res = "daily",
    temp_agg_fun = "mean",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("daily_long" %in% names(result))
  
  # Check that we have 2 daily values
  expect_equal(nrow(result$daily_long), 2)
  
  # Check date column exists
  expect_true("date" %in% names(result$daily_long))
  
  if (should_print_banners()) {
    cat("✓ Hourly input processed successfully\n")
  }
})


# ===============================================================================
# TEST 2: Daily input temporal resolution
# ===============================================================================

test_that("r2e2 handles daily input data correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 2: Daily input resolution ===\n")
  }
  
  # Create daily raster (7 days)
  env_rast <- create_temporal_raster("daily", n_periods = 7)
  geometry <- create_test_geometry()
  
  # Run r2e2 with daily data, output monthly
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "polynomial",
    trans_args = list(degree = 2),
    out_temp_res = "daily",
    temp_agg_fun = "mean",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("daily_long" %in% names(result))
  
  # Check that we have 7 daily values
  expect_equal(nrow(result$daily_long), 7)
  
  # Check transformation columns exist
  expect_true("degree_1" %in% names(result$daily_long))
  expect_true("degree_2" %in% names(result$daily_long))
  
  if (should_print_banners()) {
    cat("✓ Daily input processed successfully\n")
  }
})


# ===============================================================================
# TEST 3: Monthly input temporal resolution
# ===============================================================================

test_that("r2e2 handles monthly input data correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 3: Monthly input resolution ===\n")
  }
  
  # Create monthly raster (12 months)
  env_rast <- create_temporal_raster("monthly", n_periods = 12)
  geometry <- create_test_geometry()
  
  # Run r2e2 with monthly data
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "none",
    out_temp_res = "monthly",
    temp_agg_fun = "mean",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("monthly_long" %in% names(result))
  
  # Check that we have 12 monthly values
  expect_equal(nrow(result$monthly_long), 12)
  
  # Check date columns exist
  expect_true("date" %in% names(result$monthly_long))
  expect_true("year" %in% names(result$monthly_long))
  expect_true("month" %in% names(result$monthly_long))
  
  if (should_print_banners()) {
    cat("✓ Monthly input processed successfully\n")
  }
})


# ===============================================================================
# TEST 4: Yearly input temporal resolution
# ===============================================================================

test_that("r2e2 handles yearly input data correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 4: Yearly input resolution ===\n")
  }
  
  # Create yearly raster (3 years)
  env_rast <- create_temporal_raster("yearly", n_periods = 3)
  geometry <- create_test_geometry()
  
  # Run r2e2 with yearly data
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "bin",
    trans_args = list(breaks = c(0, 15, 25, 40)),
    out_temp_res = "yearly",
    temp_agg_fun = "mean",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("yearly_long" %in% names(result))
  
  # Check that we have 3 yearly values
  expect_equal(nrow(result$yearly_long), 3)
  
  # Check year column exists
  expect_true("year" %in% names(result$yearly_long))
  
  # Check bin columns exist (bin names use underscore format)
  result_cols <- names(result$yearly_long)
  expect_true(any(grepl("^bin_", result_cols)))
  
  if (should_print_banners()) {
    cat("✓ Yearly input processed successfully\n")
  }
})


# ===============================================================================
# TEST 5: Secondary weight raster - Daily resolution
# ===============================================================================

test_that("r2e2 handles daily secondary weight raster correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 5: Daily secondary weight raster ===\n")
  }
  
  # Create daily env raster (7 days)
  env_rast <- create_temporal_raster("daily", n_periods = 7)
  
  # Create daily secondary weight raster (7 days, same dates)
  sec_weight_rast <- create_temporal_raster("daily", n_periods = 7, seed = 456)
  
  geometry <- create_test_geometry()
  
  # Run r2e2 with secondary weights
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "none",
    out_temp_res = "daily",
    temp_agg_fun = "mean",
    sec_weight_rast = sec_weight_rast,
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("daily_long" %in% names(result))
  
  # Check that we have 7 daily values
  expect_equal(nrow(result$daily_long), 7)
  
  if (should_print_banners()) {
    cat("✓ Daily secondary weight raster processed successfully\n")
  }
})


# ===============================================================================
# TEST 6: Secondary weight raster - Monthly resolution
# ===============================================================================

test_that("r2e2 handles monthly secondary weight raster correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 6: Monthly secondary weight raster ===\n")
  }
  
  # Create daily env raster (31 days = 1 month of Jan 2000)
  env_rast <- create_temporal_raster("daily", n_periods = 31)
  
  # Create monthly secondary weight raster (1 month)
  sec_weight_rast <- create_temporal_raster("monthly", n_periods = 1, seed = 789)
  
  geometry <- create_test_geometry()
  
  # Run r2e2 with monthly secondary weights
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "polynomial",
    trans_args = list(degree = 2),
    out_temp_res = "monthly",
    temp_agg_fun = "mean",
    sec_weight_rast = sec_weight_rast,
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("monthly_long" %in% names(result))
  
  # Check that we have 1 monthly value
  expect_equal(nrow(result$monthly_long), 1)
  
  if (should_print_banners()) {
    cat("✓ Monthly secondary weight raster processed successfully\n")
  }
})


# ===============================================================================
# TEST 7: Secondary weight raster - Yearly resolution
# ===============================================================================

test_that("r2e2 handles yearly secondary weight raster correctly", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 7: Yearly secondary weight raster ===\n")
  }
  
  # Create daily env raster for one year (365 days) to aggregate to yearly
  env_rast <- create_temporal_raster("daily", n_periods = 365)
  
  # Create yearly secondary weight raster - single layer with year format
  ext_obj <- ext(0, 3, 0, 3)
  set.seed(999)
  values <- runif(9, min = 0.1, max = 10)
  sec_weight_rast <- rast(ncols = 3, nrows = 3, ext = ext_obj, 
                          vals = values, crs = "EPSG:4326")
  names(sec_weight_rast) <- "2000"  # Use year format
  
  geometry <- create_test_geometry()
  
  # Run r2e2 with yearly secondary weights
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "none",
    out_temp_res = "yearly",
    temp_agg_fun = "mean",
    sec_weight_rast = sec_weight_rast,
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result exists
  expect_true("yearly_long" %in% names(result))
  
  # Check that we have 1 yearly value
  expect_equal(nrow(result$yearly_long), 1)
  
  if (should_print_banners()) {
    cat("✓ Yearly secondary weight raster processed successfully\n")
  }
})


# ===============================================================================
# TEST 8: Subdaily to daily aggregation (hourly -> daily)
# ===============================================================================

test_that("r2e2 correctly aggregates hourly to daily with mean", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 8: Hourly to daily aggregation ===\n")
  }
  
  # Create hourly raster with constant values per day for easier verification
  ext_obj <- ext(0, 3, 0, 3)
  rast_list <- list()
  
  # Create 48 hours (2 days) with known values
  # Day 1: all hours = 10, Day 2: all hours = 20
  for (i in 1:48) {
    day_value <- if (i <= 24) 10 else 20
    rast_layer <- rast(ncols = 3, nrows = 3, ext = ext_obj, 
                       vals = rep(day_value, 9), crs = "EPSG:4326")
    
    # Format as hourly datetime
    datetime <- as.POSIXct("2000-01-01 00:00", tz = "UTC") + (i - 1) * 3600
    names(rast_layer) <- format(datetime, "%Y-%m-%d %H:%M")
    rast_list[[i]] <- rast_layer
  }
  
  env_rast <- rast(rast_list)
  geometry <- create_test_geometry()
  
  # Run r2e2 with hourly data, aggregate to daily mean
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "none",
    daily_agg_fun = "mean",
    out_temp_res = "daily",
    temp_agg_fun = "mean",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result has 2 daily values
  expect_equal(nrow(result$daily_long), 2)
  
  # Check that the mean values are correct (10 for day 1, 20 for day 2)
  expect_equal(result$daily_long$value[1], 10, tolerance = 1e-10)
  expect_equal(result$daily_long$value[2], 20, tolerance = 1e-10)
  
  if (should_print_banners()) {
    cat("✓ Hourly to daily aggregation computed correctly\n")
  }
})


# ===============================================================================
# TEST 9: Subdaily to daily aggregation with sum
# ===============================================================================

test_that("r2e2 correctly aggregates hourly to daily with sum", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 9: Hourly to daily aggregation with sum ===\n")
  }
  
  # Create hourly raster with value = 1 for all hours
  ext_obj <- ext(0, 3, 0, 3)
  rast_list <- list()
  
  # Create 24 hours (1 day) with all values = 1
  for (i in 1:24) {
    rast_layer <- rast(ncols = 3, nrows = 3, ext = ext_obj, 
                       vals = rep(1, 9), crs = "EPSG:4326")
    
    datetime <- as.POSIXct("2000-01-01 00:00", tz = "UTC") + (i - 1) * 3600
    names(rast_layer) <- format(datetime, "%Y-%m-%d %H:%M")
    rast_list[[i]] <- rast_layer
  }
  
  env_rast <- rast(rast_list)
  geometry <- create_test_geometry()
  
  # Run r2e2 with hourly data, aggregate to daily sum
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "none",
    daily_agg_fun = "sum",
    out_temp_res = "daily",
    temp_agg_fun = "sum",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result has 1 daily value
  expect_equal(nrow(result$daily_long), 1)
  
  # Check that the sum is 24 (24 hours * 1)
  expect_equal(result$daily_long$value[1], 24, tolerance = 1e-10)
  
  if (should_print_banners()) {
    cat("✓ Hourly to daily sum aggregation computed correctly\n")
  }
})


# ===============================================================================
# TEST 10: Mixed temporal aggregations (daily -> monthly -> yearly)
# ===============================================================================

test_that("r2e2 correctly aggregates across multiple temporal scales", {
  
  if (should_print_banners()) {
    cat("\n=== TEST 10: Multi-scale temporal aggregation ===\n")
  }
  
  # Create daily raster for one year (365 days)
  env_rast <- create_temporal_raster("daily", n_periods = 365)
  geometry <- create_test_geometry()
  
  # Run r2e2 aggregating daily -> yearly
  result <- r2e2(
    env_rast = env_rast,
    geometry = geometry,
    geom_id_col = "geom_id",
    trans_type = "polynomial",
    trans_args = list(degree = 2),
    out_temp_res = "yearly",
    temp_agg_fun = "mean",
    verbose = 0,
    validation = FALSE
  )
  
  # Check that result has 1 yearly value
  expect_equal(nrow(result$yearly_long), 1)
  
  # Check that year column exists
  expect_true("year" %in% names(result$yearly_long))
  
  if (should_print_banners()) {
    cat("✓ Multi-scale temporal aggregation works correctly\n")
  }
})


# ===============================================================================
# Summary Banner
# ===============================================================================

if (should_print_banners()) {
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("TEMPORAL RESOLUTION TESTS COMPLETED\n")
  cat(strrep("=", 80), "\n")
}
