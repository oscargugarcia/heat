#' @importFrom terra rast window time nlyr
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr left_join rename select sym
#' @importFrom lubridate year
#' @importFrom arrow write_parquet
NULL


#' Run the Climate Data Pipeline
#'
#' This function runs an aggregation pipeline that processes environmental raster data
#' and polygon datasets to perform spatial and temporal aggregations. It loads the necessary
#' data, applies specified transformations, and outputs aggregated results in both wide
#' and long formats. The function also allows for the calculation of area weights if
#' specified.
#'
#' @param env_rast Either a character string specifying the directory path to environmental 
#'                  rasters, OR a SpatRaster object containing the environmental data. 
#'                  Layer names must be in date format (e.g., "YYYY-MM-DD").
#'
#' @param polygons Either a character string specifying the file path to a spatial polygon file
#'                 (e.g., shapefile, geopackage), OR an sf/SpatVector object containing the 
#'                 polygons for spatial aggregation.
#'
#' @param poly_id_col A character string representing the column name of the unique ID
#'                     in the polygon dataset (e.g., "GEOID"). If NULL, a row_index
#'                     column will be created and used as the polygon ID, 
#'                     and all other columns are kept in the output.
#'
#' @param trans_type A character string indicating the type of transformation to apply
#'                   to the raster data (e.g., "polynomial", "bin", etc.").
#'
#' @param trans_type A character string indicating the type of transformation to apply
#'                   to the raster data (e.g., "polynomial", "bin", etc.).
#'
#' @param trans_args A list of additional parameters to be passed to the transformation
#'                   function (e.g., polynomial degree).
#'
#' @param temp_agg_args A list containing parameters for temporal aggregation, including
#'                       the output temporal resolution (out_temp_res) and aggregation function.
#'
#' @param sec_weight_rast Either a character string specifying the directory path to secondary 
#'                        weight rasters, OR a SpatRaster object containing the weight data. 
#'                        If NULL (default), no secondary weighting is applied. 
#'                        Layer names must be in date format.
#'
#' @param boundary_dates A vector of Date objects defining the start and end dates for
#'                       the weighting periods.
#'
#' @param spatial_agg_args A list containing any specific arguments related to spatial
#'                         aggregation, including optional append_cols for additional columns
#'                         to include in output (default is NULL).
#'
#' @param save_path Optional character string specifying the output directory for storing 
#'                  results. If NULL (default), results will be returned without saving to disk.
#'                  When provided, area weights are automatically saved.
#'                       
#' @param out_format Character string specifying the output format. Options are:
#'                   - "all" (default): Returns and saves both wide and long formats
#'                   - "wide": Returns and saves only wide format (time steps as columns)
#'                   - "long": Returns and saves only long format (time steps as rows)
#'                   Wide format has one row per polygon with time steps as separate columns.
#'                   Long format has one row per polygon-time combination with time in rows.
#'                       
#' @param run_control_checks A logical value indicating whether to generate control check
#'                  reports for data validation (default is TRUE).
#'
#' @param run_control_checks A logical value indicating whether to generate control check
#'                  reports for data validation (default is TRUE).
#'                       
#' @param max_cells A numeric value representing the maximum number of raster cells that
#'                  can be processed at once (default is 3e7).
#'
#' @param metadata A list containing metadata related to the dataset and analysis,
#'                 such as variable names, descriptions, and additional identifiers (default is NULL).
#'
#' @param save_console_output A logical value indicating whether to save console output
#'                           to a log file (default is FALSE).
#'                           
#' @param save_batch_output Logical indicating whether to save intermediate batch outputs to disk for memory efficiency (default: TRUE).
#' 
#' @param overwrite_batch_output Logical indicating whether to overwrite existing batch files. When FALSE, existing files
#'                                are validated to ensure layer counts match current processing parameters (default: FALSE).
#'
#' @param compression A character string specifying the compression algorithm for parquet
#'                    files (default is 'zstd').
#'                       
#' @param control_checks_rmd A character string specifying the name of the RMD template
#'                          file to look for in the utils folder (default is "02_control_checks_source_template.Rmd")
#'                          
#' @param selected_trans_var The name of the transformation variable to be used in the report.
#'                           (default is 'degree_1' but if not found, it will use the first detected transformation variable)
#'                           
#' @param trans_var_description A character string describing the transformation variable
#'                             for plot labels (default is "Temperature (C°)")
#'
#'
#' @details The function performs the following main operations:
#' - Loads spatial and environmental data.
#' - Processes and applies transformations to raster data.
#' - Executes spatial and temporal aggregations with optional secondary weighting.
#' - Saves the results in both wide and long formats along with optional area weights.
#'
#' @return A list containing the processed data frames:
#'         - spatial_agg: Daily aggregated data in wide format (time steps as columns, NULL if out_format = "long")
#'         - spatial_agg_long: Daily aggregated data in long format (time steps as rows, NULL if out_format = "wide" or temporal resolution is daily)
#'         - temp_agg_wide: Temporally aggregated data in wide format (time steps as columns, NULL if out_format = "long" or temporal resolution is daily)
#'         - temp_agg_long: Temporally aggregated data in long format (time steps as rows, NULL if out_format = "wide")
#'         - area_weights: Area weights data frame (always included)
#'         Results are also saved directly to the specified output directory in parquet format when save_path is provided.
r2e2 <- function(env_rast,
                 polygons,
                 poly_id_col = NULL,
                 trans_type,
                 trans_args,
                 temp_agg_args,
                 sec_weight_rast = NULL,
                 boundary_dates = NULL,
                 spatial_agg_args = NULL,
                 save_path = NULL,
                 out_format = "all",
                 run_control_checks = TRUE,
                 max_cells = 3e7,
                 metadata = NULL,
                 save_console_output = FALSE,
                 save_batch_output = TRUE,
                 overwrite_batch_output = FALSE,
                 compression = 'zstd',
                 control_checks_rmd = "02_control_checks_source_template.Rmd",
                 selected_trans_var = 'degree_1',
                 trans_var_description = "Temperature (C°)") {
  
  # Check if output path is provided
  save_path_provided <- !is.null(save_path) && save_path != ""
  
  # Infer sec_weights from whether sec_weight_rast is provided
  sec_weights <- !is.null(sec_weight_rast)
  
  # Validate out_format parameter
  if (!out_format %in% c("all", "wide", "long")) {
    stop("out_format must be one of: 'all', 'wide', or 'long'")
  }
  
  # Determine what to save/return based on out_format
  save_wide <- out_format %in% c("all", "wide")
  save_long <- out_format %in% c("all", "long")
  
  # Ensure output directory exists, create it if it does not exist
  if (save_path_provided && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  if (!save_path_provided) {
    message("No save_path provided. Results will be returned without saving.")
    # Disable operations that require file output
    save_console_output <- FALSE
    run_control_checks <- FALSE
    save_batch_output <- FALSE
  }
  
  if (save_console_output && save_path_provided) {
    message("To view console output, see 'console_output.log' in the output folder.")
    # Create a log file to store the console output
    log_file_path <- file.path(save_path, "console_output.log")
    log_file <- file(log_file_path)
    
    # Start capturing the console output
    sink(log_file, append=TRUE, split = TRUE)
    sink(log_file, append=TRUE, type="message")
  }
  
  
  # ---- 1. Data --------------------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 1. Load Data",
          "\n===========================================================================\n")
  
  ## ---- 1.1 Create output folder and console output --------------------------------------------------------------------
  
  if (save_path_provided) {
    # Combine all global parameter objects into a single list.
    global_params <- list(
      env_rast_type = if (is.character(env_rast)) "path" else if (inherits(env_rast, "SpatRaster")) "SpatRaster" else "NULL",
      env_rast_path = if (is.character(env_rast)) env_rast else NULL,
      polygons_type = if (is.character(polygons)) "path" else if (inherits(polygons, c("sf", "SpatVector"))) class(polygons)[1] else "NULL",
      polygons_path = if (is.character(polygons)) polygons else NULL,
      sec_weight_rast_type = if (is.character(sec_weight_rast)) "path" else if (inherits(sec_weight_rast, "SpatRaster")) "SpatRaster" else "NULL",
      sec_weight_rast_path = if (is.character(sec_weight_rast)) sec_weight_rast else NULL,
      save_path = save_path,
      trans_type = trans_type,
      poly_id_col = poly_id_col,
      sec_weights = sec_weights,  # Inferred from sec_weight_rast
      out_format = out_format,
      boundary_dates = boundary_dates,
      trans_args = trans_args,
      spatial_agg_args = spatial_agg_args,
      temp_agg_args = temp_agg_args,
      metadata = metadata
    )
    
    # Capture the output of dput() to get a complete code representation.
    output <- capture.output(dput(global_params))
    
    # Write to file
    writeLines(output, file.path(save_path, "global_parameters.txt"))
    
    message("✓ Global parameters saved")
  }

  
  ## ---- 1.2 Import Polygons ---------------------------------------------------------
  
  # Handle polygons - either load from path or validate existing object
  if (is.character(polygons)) {
    polygons <- read_spatial_file(polygons)
    message("✓ Polygons loaded: ", nrow(polygons), " features")
  } else if (!inherits(polygons, c("sf", "SpatVector"))) {
    stop("polygons must be either a character path to a spatial file or an sf/SpatVector object")
  } else {
    # If a SpatVector is provided, convert to sf first
    if (inherits(polygons, "SpatVector")) {
      polygons <- sf::st_as_sf(polygons)
    }
    # Validate and clean the provided sf object
    polygons <- check_spatial_file(polygons)
    message("✓ Geometry provided and validated: ", nrow(polygons), " rows")
  }
  
  # Handle NULL poly_id_col by creating row_index and updating spatial_agg_args
  if (is.null(poly_id_col)) {
    message("No poly_id_col provided. Using row index as polygon ID column. All original columns will be kept in the output.")
    
    # Get column names of the polygons
    polygon_colnames <- polygons %>% 
      st_drop_geometry() %>% 
      colnames()
    
    # Add all original column names (except row_index) to append_cols
    spatial_agg_args$append_cols <- polygon_colnames
    
    # Add row_index column to polygons
    polygons$row_index <- seq_len(nrow(polygons))
  
    # Assign row_index as the polygon ID column
    poly_id_col <- "row_index"
    
    # Initialize spatial_agg_args if it's NULL
    if (is.null(spatial_agg_args)) {
      spatial_agg_args <- list()
    }
  }


  
  ## ---- 1.3 Import environmental Raster ---------------------------------------------------------
  
  # Handle env_rast - either load from path or validate existing raster
  if (is.character(env_rast)) {
    # Read the environmental raster files from the specified path
    env_rast_files <- filter_read_rast(env_rast)
    
    # Read the environmental rasters into a SpatRaster
    env_rast <- rast(file.path(env_rast, env_rast_files))
  } else if (!inherits(env_rast, "SpatRaster")) {
    stop("env_rast must be either a character path or a SpatRaster object")
  }
  
  # Validate raster layer names
  env_rast <- check_raster_layers(env_rast, "Environmental raster")
  
  # Create a buffered extent around the input polygons using the raster's resolution
  buffered_extent <- buffer_polygons(env_rast, polygons, buffer_factor = 1)
  
  # Crop the extent of the raster to the polygon buffer (terra::window is faster than terra::crop. To undo the crop, load the raster back in with terra::rast())
  terra::window(env_rast) <- buffered_extent
  # Clean up large objects no longer needed
  gc(verbose = FALSE)


  ## ---- 1.4 Handle boundary dates ---------------------------------------------------------
  
  if (is.null(boundary_dates)) {
    # No boundary dates provided at all
    boundary_dates <- parse_boundary_dates(env_rast)
    message("boundary_dates not provided. Using first and last dates from the environmental raster: ",
            paste(boundary_dates, collapse = " to "))
    
  } else {
    # Get dates from raster for potential replacement
    raster_dates <- parse_boundary_dates(env_rast)
    
    # Check and replace start date
    if (is.na(boundary_dates[1]) || isTRUE(boundary_dates[1] == "")) {
      boundary_dates[1] <- as.Date(raster_dates[1])
      message("Start date in boundary_dates missing, proceeding with raster start date: ", raster_dates[1])
    }
    
    # Check and replace end date
    if (is.na(boundary_dates[2]) || isTRUE(boundary_dates[2] == "")) {
      boundary_dates[2] <- as.Date(raster_dates[2])
      message("End date in boundary_dates missing, proceeding with raster end date: ", raster_dates[2])
    }
  }
  
  boundary_dates <- c(boundary_dates[1], boundary_dates[2])


  ## ---- 1.5 Import secondary Weight Raster ----------------------------------------------------------
  
  if (sec_weights) {
    # Handle sec_weight_rast - either load from path or validate existing raster
    if (is.character(sec_weight_rast)) {
      sec_weight_files <- filter_read_rast(sec_weight_rast)
      if (length(sec_weight_files) != 0) {
        sec_weight_rast <- rast(file.path(sec_weight_rast, sec_weight_files))
      } else {
        stop("No secondary weight rasters found in the specified path.")
      }
    } else if (inherits(sec_weight_rast, "SpatRaster")) {
      # Raster object provided directly - already loaded
      # No need to do anything, sec_weight_rast is already set
    } else {
      stop("sec_weight_rast must be either a character path or a SpatRaster object")
    }
    
    # Validate raster layer names
    sec_weight_rast <- check_raster_layers(sec_weight_rast, "Secondary weight raster")
        
    # get the boundary dates for the secondary weighting periods and relevant sec_weight_layers
    period_results <- get_period_boundaries(sec_weight_rast, boundary_dates)
    period_boundaries <- period_results$boundaries
    relevant_sec_weight_indices <- period_results$indices
    
    # Create a list of SpatRaster layers for each period (subset by indices)
    sec_weight_rast_list <- lapply(relevant_sec_weight_indices, function(idx) {
      sec_weight_rast[[idx]]
    })
    
    # For verification/display purposes, get layer names
    sec_weight_layer_names <- names(sec_weight_rast)[relevant_sec_weight_indices]
    
  } else {
    message("No secondary weight raster provided. Proceeding without secondary weighting.")
    sec_weight_rast_list = NULL
    sec_weight_layer_names = NULL
    period_boundaries <-  boundary_dates
  }
  

  
  # ---- 2. Prepare Aggregation ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 2. Prepare Aggregation",
          "\n===========================================================================\n")
  
  ## ---- 2.1 Assign secondary weighting periods ----------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   2.1 Assign Periods",
          "\n---------------------------------------------------\n")
  
  # Generate [start_date, end_date] pairs for each period
  period_defs <- create_periods(period_boundaries)
  
  # Read and process raster stacks for each period
  env_rast_list <- assign_weighting_periods(env_rast, period_defs)
 
  sec_weight_name <- ifelse(is.null(metadata), "name not provided", metadata$sec_weight_product)
  weighting_periods <- check_periods(env_rast_list, sec_weight_layer_names, sec_weights, sec_weight_name)

  ## ---- 2.2 Transformation functions and arguments ---------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   2.2 Select Transformation Functions",
          "\n---------------------------------------------------\n")
  
  # Select the transformation function (e.g., stats::poly for "polynomial")
  trans_fun <- select_trans_fun(trans_type)
  
  # Check the function and its arguments
  checked_trans_args <- check_trans(trans_fun, trans_args)
  
  ## ---- 2.3 Spatial aggregation arguments ---------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   2.3 Check Spatial Aggregation Arguments",
          "\n---------------------------------------------------\n")
  
  # Get appended columns argument if provided
  appended_cols <- spatial_agg_args$append_cols
  
  # Check spatial aggregation arguments
  checked_spatial_agg_args <- check_spatial_agg_args(spatial_agg_args)
  
  # ---- 3. Transformation and Spatial Aggregation ---------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 3. Transformation and Spatial Aggregation",
          "\n===========================================================================")
  
  message("\n---------------------------------------------------",
          "\n   3.1 Run Transformation and Spatial Aggregation",
          "\n---------------------------------------------------\n")
  
  ## ---- 3.1 Process all periods sequentially ----------------------------------------------------------
  
  spatial_agg <- trans_spatial_agg(env_rast_list = env_rast_list, 
                                   sec_weight_rast_list = sec_weight_rast_list, 
                                   polygons = polygons,
                                   buffered_extent = buffered_extent, 
                                   trans_fun = trans_fun, 
                                   trans_type = trans_type,
                                   checked_trans_args = checked_trans_args,
                                   spatial_agg_args = checked_spatial_agg_args, 
                                   poly_id_col = poly_id_col, 
                                   weighting_periods = weighting_periods, 
                                   save_path = save_path, 
                                   sec_weights = sec_weights,
                                   max_cells = max_cells, 
                                   save_batch_output =  save_batch_output,
                                   overwrite_batch_output = overwrite_batch_output)
  
  # Clean up large objects no longer needed
  rm(env_rast_list)  # Remove the list of raster periods
  gc(verbose = FALSE)
  
  ## ---- 3.2 Process spatial aggregation output ----------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   3.2 Check and Process Spatial Aggregation Output",
          "\n---------------------------------------------------\n")
  
  # Check for any polygons with missing values
  date_cols <- get_date_cols(spatial_agg)
  poly_ids_missing_vals <- get_missing_vals(spatial_agg, date_cols)
  
  # Rename transformation variable
  spatial_agg <- rename_trans_var(spatial_agg, trans_type, checked_trans_args)
  
  # Add metadata columns
  spatial_agg <- add_metadata_cols(spatial_agg, trans_type, metadata, trans_args)
  
  # Note: No geometry is joined here - spatial_agg remains a regular data frame
  # The poly_id column already exists from the spatial aggregation step
  
  ## ---- 3.3 Save output ----------------------------------------------------------
  
  if (save_wide && save_path_provided) {
    message("\n---------------------------------------------------",
            "\n   3.3 Save Aggregated Output",
            "\n---------------------------------------------------\n")
    
    # Save as parquet (no geometry in output)
    write_parquet(spatial_agg, file.path(save_path, "aggregation_output_daily_wide.parquet"), 
                  compression = compression)
    message("Daily aggregation output in wide format saved to: aggregation_output_daily_wide.parquet")
  }
  
  ## ---- 3.4 Get area weights ----------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   3.4 Calculate Area Weights",
          "\n---------------------------------------------------\n")
  
  area_weights <- get_area_weights(env_rast[[1]], polygons, poly_id_col)
  if (exists("env_rast")) rm(env_rast)  # Remove original raster

  if (save_path_provided) {
     write_parquet(area_weights, file.path(save_path, "area_weights.parquet"))
     message("Area weights saved to: area_weights.parquet")
  } else {
    message("Area weights calculated but not saved (no output path provided).")
  }


  # ---- 4. Temporal Aggregation ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 4. Temporal Aggregation",
          "\n===========================================================================\n")
  
  temp_agg_wide <- temp_agg(spatial_agg, temp_agg_args, keep_metadata = TRUE)
  
  if (save_wide && save_path_provided) {
    # Save as parquet (no geometry in output)
    filename_temp_agg_wide <- paste0("aggregation_output_", temp_agg_args$out_temp_res, "_wide.parquet")
    write_parquet(temp_agg_wide, file.path(save_path, filename_temp_agg_wide), 
                  compression = compression)
    message("Temporally aggregated output in wide format saved to: aggregation_output_", 
            temp_agg_args$out_temp_res, "_wide.parquet")
  }
  
  # ---- 5. Reshape to Long Format ----------------------------------------------------------
  
  # Only perform long format conversion if requested
  if (save_long) {
    message("\n===========================================================================",
            "\n 5. Reshape Data to Long Format",
            "\n===========================================================================")
    
    message("\n---------------------------------------------------",
            "\n   5.1 Reshape Temporally Aggregated Output to Long Format",
            "\n---------------------------------------------------\n")
    
    # Save long temporally aggregated output
    temp_agg_long <- reshape_to_long(temp_agg_wide, add_time_columns = TRUE)
    temp_agg_long <- add_appended_cols(temp_agg_long, polygons, poly_id_col, appended_cols)
    
    if (save_path_provided) {
      # Save as parquet (no geometry in output)
      filename_temp_agg_long <- paste0("aggregation_output_", temp_agg_args$out_temp_res, "_long.parquet")
      write_parquet(temp_agg_long, file.path(save_path, filename_temp_agg_long), 
                    compression = compression)
      message("Temporally aggregated output in long format saved to: aggregation_output_", 
              temp_agg_args$out_temp_res, "_long.parquet")
    }
    
    # Save long daily output
    if (temp_agg_args$out_temp_res != "daily") {
      message("\n---------------------------------------------------",
              "\n   5.2 Reshape Daily Output to Long Format",
              "\n---------------------------------------------------\n")
      
      spatial_agg_long <- reshape_to_long(spatial_agg, add_time_columns = TRUE)
      spatial_agg_long <- add_appended_cols(spatial_agg_long, polygons, poly_id_col, appended_cols)
      
      if (save_path_provided) {
        # Save as parquet (no geometry in output)
        filename_spatial_agg_long <- "aggregation_output_daily_long.parquet"
        write_parquet(spatial_agg_long, file.path(save_path, filename_spatial_agg_long), 
                      compression = compression)
        message("Daily aggregation output in long format saved to: aggregation_output_daily_long.parquet")
      }
    }
    
    # Clean up objects based on format after all processing is complete
    if (!save_wide) {
      # Remove wide format objects if only long format is needed
      rm(spatial_agg, temp_agg_wide)
      gc(verbose = FALSE)
    }
  } else {
    message("\nSkipping long format conversion (out_format = 'wide').\n")
    temp_agg_long <- NULL
    spatial_agg_long <- NULL
  }
  
  # Collect the relevant data frames based on out_format
  results <- list(
    spatial_agg = if (save_wide) spatial_agg else NULL,
    spatial_agg_long = if (save_long && temp_agg_args$out_temp_res != "daily") spatial_agg_long else NULL,
    temp_agg_wide = if (save_wide && temp_agg_args$out_temp_res != "daily") temp_agg_wide else NULL,
    temp_agg_long = if (save_long) temp_agg_long else NULL,
    area_weights = area_weights
  )
  
  # Clean up objects not included in results to save memory
  if (!save_wide) {
    if (exists("spatial_agg")) rm(spatial_agg)
    if (exists("temp_agg_wide")) rm(temp_agg_wide)
  }
  if (!save_long) {
    if (exists("spatial_agg_long")) rm(spatial_agg_long)
    if (exists("temp_agg_long")) rm(temp_agg_long)
  }
  gc(verbose = FALSE)
  
  # ---- 6. Control Checks ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 6. Control Checks",
          "\n===========================================================================\n")
  
  # Check if control checks should be run
  if (!run_control_checks) {
    if (!save_path_provided) {
      message("Control checks skipped (no output path provided).\n")
    } else {
      message("Control checks skipped. If control checks are desired, set run_control_checks = TRUE.\n")
    }
    return(results)
  } else {
    message("\nRunning control checks...\n")
  
  # Run control checks on the long format data
  control_checks(
    df_long = temp_agg_long,
    polygons = polygons,
    area_weights = area_weights,
    paths = paths,
    poly_id_col = poly_id_col,
    selected_trans_var = selected_trans_var,
    trans_var_description = trans_var_description,
    control_checks_rmd = control_checks_rmd)
  }
  
  # ---- 7. Save Console Output ----------------------------------------------------------
  
  if (save_console_output) {
    # Stop capturing the output
    sink() 
    sink(type="message")
  }
  
  # Return the results as a list
  return(results)

}