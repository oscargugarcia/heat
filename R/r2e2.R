#' Run the CLimate Data Pipeline
#'
#' This function runs an aggregation pipeline that processes environmental raster data
#' and polygon datasets to perform spatial and temporal aggregations. It loads the necessary
#' data, applies specified transformations, and outputs aggregated results in both wide
#' and long formats. The function also allows for the calculation of area weights if
#' specified.
#'
#' @param paths A list containing paths for input and output files, including:
#'              - path_utils: Directory path containing pipeline functions
#'              - path_out_folder: Output directory for storing results.
#'              - path_polygons: File path to the polygons (e.g., shapefile).
#'              - path_env_rast: Directory path for input environmental rasters.
#'              - path_sec_weight_rast: Directory path for secondary weight rasters.
#'              - path_tmin_rast: Directory path for minimum temperature rasters (alternative to path_env_rast).
#'              - path_tmax_rast: Directory path for maximum temperature rasters (alternative to path_env_rast).
#'
#' @param trans_type A character string indicating the type of transformation to apply
#'                   to the raster data (e.g., "polynomial", "bin", etc.).
#'
#' @param polygons Spatial polygon object (sf or SpatVector) containing the polygons
#'                 for spatial aggregation.
#'
#' @param poly_id_col A character string representing the column name of the unique ID
#'                     in the polygon dataset (e.g., "GEOID"). If NULL, a row_index
#'                     column will be created and used as the polygon ID, 
#'                     and all other columns are kept in the output.
#'
#' @param boundary_dates A vector of Date objects defining the start and end dates for
#'                       the weighting periods.
#'
#' @param sec_weights A logical value indicating whether to apply secondary weights during
#'                    the analysis (default is FALSE).
#'
#' @param daily_agg_fun A character string specifying the daily aggregation function to apply
#'                      to sub-daily data. Use "none" to skip daily aggregation (default is "none").
#'
#' @param interpol_fun A character string specifying the interpolation function for temperature
#'                     data when using tmin/tmax inputs. Use "none" to skip interpolation (default is "none").
#'
#' @param interpol_args A list of arguments to pass to the interpolation function when
#'                      processing tmin/tmax data (default is NULL).
#'
#' @param trans_args A list of additional parameters to be passed to the transformation
#'                   function (e.g., polynomial degree).
#'
#' @param spatial_agg_args A list containing any specific arguments related to spatial
#'                         aggregation, including optional append_cols for additional columns
#'                         to include in output (default is NULL).
#'
#' @param temp_agg_args A list containing parameters for temporal aggregation, including
#'                       the output temporal resolution (out_temp_res) and aggregation function.
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
#' @param save_area_weights A logical value indicating whether to calculate and save area
#'                         weights in the output (default is TRUE).
#'
#' @param save_wide_output A logical value indicating whether to save the aggregated
#'                         dataframes to parquet format in wide format (default is TRUE).
#'

#' @param save_long_daily_output A logical value indicating whether to save the daily
#'                         output to parquet format in long format (default is TRUE).
#'
#' @param save_interpol_rast A logical value indicating whether to save interpolated
#'                          raster data when processing tmin/tmax inputs (default is FALSE).
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
#' - Loads spatial and environmental data (either single environmental raster or tmin/tmax pairs).
#' - Processes and applies transformations to raster data.
#' - Executes spatial and temporal aggregations with optional secondary weighting.
#' - Saves the results in both wide and long formats along with optional area weights.
#' - Supports interpolation of sub-daily temperature data from daily tmin/tmax inputs.
#' - Provides flexible daily aggregation for sub-daily input data.
#'
#' @return A list containing the processed data frames:
#'         - spatial_agg: Daily aggregated data in wide format
#'         - spatial_agg_long: Daily aggregated data in long format (if save_long_daily_output is TRUE and temporal resolution is not daily)
#'         - temp_agg_wide: Temporally aggregated data in wide format (if temporal resolution is not daily)
#'         - temp_agg_long: Temporally aggregated data in long format
#'         - area_weights: Area weights data frame (if save_area_weights is TRUE)
#'         Results are also saved directly to the specified output directory in parquet format.
r2e2 <- function(paths, trans_type, polygons, poly_id_col = NULL, boundary_dates = NULL,
                         sec_weights = FALSE, 
                         daily_agg_fun = "none", interpol_fun = "none", interpol_args = NULL,
                         trans_args, spatial_agg_args = NULL, temp_agg_args,
                         run_control_checks = TRUE,
                         max_cells = 3e7, metadata = NULL, save_area_weights = TRUE, save_wide_output = TRUE, save_long_daily_output = TRUE, save_interpol_rast = FALSE, save_console_output = FALSE,  save_batch_output = TRUE,
                         overwrite_batch_output = FALSE, compression = 'zstd',
                         control_checks_rmd = "02_control_checks_source_template.Rmd", selected_trans_var = 'degree_1', trans_var_description = "Temperature (C°)") {
  
  # Check if output path is provided
  output_path_provided <- !is.null(paths$path_out_folder) && paths$path_out_folder != ""
  
  # Ensure output directory exists, create it if it does not exist
  if (output_path_provided && !dir.exists(paths$path_out_folder)) {
    dir.create(paths$path_out_folder, recursive = TRUE)
  }
  
  if (!output_path_provided) {
    message("No output path provided. Results will be returned without saving.")
    # Disable operations that require file output
    save_console_output <- FALSE
    run_control_checks <- FALSE
    save_batch_output <- FALSE
  }
  
  if (save_console_output && output_path_provided) {
    message("To view console output, see 'console_output.log' in the output folder.")
    # Create a log file to store the console output
    log_file_path <- file.path(paths$path_out_folder, "console_output.log")
    log_file <- file(log_file_path)
    
    # Start capturing the console output
    sink(log_file, append=TRUE, split = TRUE)
    sink(log_file, append=TRUE, type="message")
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


    
  # ---- 1. Data --------------------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 1. Load Data",
          "\n===========================================================================")
  
  ## ---- 1.1 Create output folder and console output --------------------------------------------------------------------
  
  if (output_path_provided) {
    message("\n---------------------------------------------------",
            "\n   1.1 Save global parameters txt file",
            "\n---------------------------------------------------\n")
    
    # Combine all global parameter objects into a single list.
    global_params <- list(
      paths = paths,
      trans_type = trans_type,
      poly_id_col = poly_id_col,
      sec_weights = sec_weights,
      save_area_weights = save_area_weights,
      boundary_dates = boundary_dates,
      trans_args = trans_args,
      spatial_agg_args = spatial_agg_args,
      temp_agg_args = temp_agg_args,
      metadata = metadata
    )
    
    # Capture the output of dput() to get a complete code representation.
    output <- capture.output(dput(global_params))
    
    # Write to file
    writeLines(output, file.path(paths$path_out_folder, "global_parameters.txt"))
    
    message("Global parameters saved to 'global_parameters.txt'.")
  }

  
  ## ---- 1.2 Import environmental Raster ---------------------------------------------------------
  
  
  if (!is.null(paths$path_env_rast)) {
    
    message("\n---------------------------------------------------",
            "\n   1.2 Load Environmental Rasters",
            "\n---------------------------------------------------\n")
    
    # Read the environmental raster files from the specified path, filtering by date range
    env_rast_files <- filter_read_rast(paths$path_env_rast)
    
    # Read the environmental rasters into a RasterStack
    env_rast <- rast(file.path(paths$path_env_rast, env_rast_files))
    
    # Get layer names and check whether they are in date format. Otherwise try using the time dimension of the raster
    layer_names <- names(env_rast)
    temp_res <- tryCatch({
      get_temp_res(layer_names)
    }, error = function(e) {
      if (grepl("Multiple temporal resolutions detected", e$message)) {
        time_dim <- terra::time(env_rast)
        if (!is.null(time_dim)) {
          message("The layer names of the raster seem to not be in a time format (e.g. 'YYYY-MM-DD'. Attempting to use time dimension of the raster")
          names(env_rast) <<- as.character(time_dim)
        } else stop(e$message)
      } else stop(e$message)
    })
    
    # Create a buffered extent around the input polygons using the raster's resolution
    buffered_extent <- buffer_polygons(env_rast, polygons, buffer_factor = 1)
    
    # Crop the extent of the raster to the polygon buffer (terra::window is faster than terra::crop. To undo the crop, load the raster back in with terra::rast())
    window(env_rast) <- buffered_extent
    # Clean up large objects no longer needed
    gc(verbose = FALSE)
    
    # Aggregation to daily if daily aggregation function is not "none"
    if (daily_agg_fun != "none") {
      env_rast <- aggregation_to_daily(env_rast, fun = daily_agg_fun)
    }
    
    message("Environmental raster read successfully.")
  } else if (!is.null(paths$path_tmin_rast) & !is.null(paths$path_tmax_rast)) {
    
    message("\n---------------------------------------------------",
            "\n   1.2 Load and Interpolate Tmin and Tmax Rasters ",
            "\n---------------------------------------------------\n")
    
   
    # Process the tmin and tmax rasters to get the subdaily raster
    env_rast <- tmin_tmax_interpolation(paths = paths, polygons = polygons, 
                                       boundary_dates = boundary_dates, 
                                       max_cells = max_cells, save_interpol_rast = save_interpol_rast,
                                       interpol_fun = interpol_fun,  interpol_args = interpol_args)

    
    # Create a buffered extent around the input polygons using the raster's resolution
    buffered_extent <- buffer_polygons(env_rast, polygons, buffer_factor = 1)
    
    message("Tmin and Tmax rasters processed successfully.")
    
  } else {
    stop("No environmental raster or tmin + tmax rasters provided.")
  }  
  
  # Handle boundary dates - use env_rast dates if any are missing
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
  
  
  
  
  
  ## ---- 1.3 Import secondary Weight Raster ----------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   1.3 Load Secondary Weight Rasters",
          "\n---------------------------------------------------\n")
  
  if (sec_weights) {
    sec_weight_rast_files <- filter_read_rast(paths$path_sec_weight_rast)
    if (length(sec_weight_rast_files) != 0) {
      sec_weight_rast <- rast(file.path(paths$path_sec_weight_rast, sec_weight_rast_files))
      
      # get the boundary dates for the secondary weighting periods and relevant sec_weight_layers
      period_results <- get_period_boundaries(sec_weight_rast, boundary_dates)
      period_boundaries <- period_results$boundaries
      relevant_sec_weight_indices <- period_results$indices
      
      # subset seconday weight raster files to indices
      sec_weight_rast_files <- sec_weight_rast_files[relevant_sec_weight_indices]
    } else {
      stop("No secondary weight rasters found in the specified path with the given year range.")
    }
    
  } else {
    message("Secondary weights set to FALSE.")
    sec_weight_rast_files = NULL
    period_boundaries <-  boundary_dates
  }
  

  
  # ---- 2. Prepare Aggregation ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 2. Prepare Aggregation",
          "\n===========================================================================")
  
  ## ---- 2.1 Assign secondary weighting periods ----------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   2.1 Assign Periods",
          "\n---------------------------------------------------\n")
  
  # Generate [start_date, end_date] pairs for each period
  period_defs <- create_periods(period_boundaries)
  
  # Read and process raster stacks for each period
  env_rast_list <- assign_weighting_periods(env_rast, period_defs)
 
  sec_weight_name <- ifelse(is.null(metadata), "name not provided", metadata$sec_weight_product)
  weighting_periods <- check_periods(env_rast_list, sec_weight_rast_files, sec_weights, sec_weight_name)

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
                                   sec_weight_rast_files = sec_weight_rast_files, 
                                   polygons = polygons,
                                   buffered_extent = buffered_extent, 
                                   trans_fun = trans_fun, 
                                   trans_type = trans_type,
                                   checked_trans_args = checked_trans_args,
                                   spatial_agg_args = checked_spatial_agg_args, 
                                   poly_id_col = poly_id_col, 
                                   weighting_periods = weighting_periods, 
                                   paths = paths, 
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
  
  # Join the geometry column from polygons to spatial_agg
  spatial_agg <- spatial_agg %>%
    left_join(polygons %>% rename(poly_id = !!sym(poly_id_col)) %>% dplyr::select(all_of("poly_id")), by = "poly_id") 
  
  ## ---- 3.3 Save output ----------------------------------------------------------
  
  if(save_wide_output && output_path_provided) {
    message("\n---------------------------------------------------",
            "\n   3.3 Save Aggregated Output",
            "\n---------------------------------------------------\n")
    
    # Save standardized dataframe to parquet with Zstandard (zstd) compression
    write_parquet(spatial_agg, file.path(paths$path_out_folder, "aggregation_output_daily_wide.parquet"), compression = compression)
    message("Daily aggregation output in wide format saved to: aggregation_output_daily_wide.parquet")
  }
  
  ## ---- 3.4 Get area weights ----------------------------------------------------------
  
  message("\n---------------------------------------------------",
          "\n   3.4 Calculate Area Weights",
          "\n---------------------------------------------------\n")
  
  area_weights <- get_area_weights(env_rast[[1]], polygons, poly_id_col)
  if (exists("env_rast")) rm(env_rast)  # Remove original raster

  if (save_area_weights && output_path_provided) {
     write_parquet(area_weights, file.path(paths$path_out_folder, "area_weights.parquet"))
     message("Area weights saved to: area_weights.parquet")
  } else if (!output_path_provided) {
    message("Area weights calculated but not saved (no output path provided).")
  } else {
    message("Area weights not saved as per user request.")
  }


  # ---- 4. Temporal Aggregation ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 4. Temporal Aggregation",
          "\n===========================================================================\n")
  
  temp_agg_wide <- temp_agg(spatial_agg, temp_agg_args, keep_metadata = TRUE)
  
  if(save_wide_output && output_path_provided) {
    filename_temp_agg_wide <- paste0("aggregation_output_", temp_agg_args$out_temp_res,  "_wide.parquet")
    write_parquet(temp_agg_wide, file.path(paths$path_out_folder, filename_temp_agg_wide), compression = compression)
    message("Temporally aggregated output in wide format saved to: aggregation_output_",temp_agg_args$out_temp_res, "_wide.parquet")
  }
  
  # ---- 5. Reshape to Long Format ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 5. Reshape Data to Long Format",
          "\n===========================================================================")
  
  message("\n---------------------------------------------------",
          "\n   5.1 Reshape Temporally Aggregated Output to Long Format",
          "\n---------------------------------------------------\n")
  
  # Save long temporally aggregated output
  temp_agg_long <- reshape_to_long(temp_agg_wide, add_time_columns = TRUE)
  temp_agg_long <- add_appended_cols(temp_agg_long, polygons, poly_id_col, appended_cols)
  
  if (output_path_provided) {
    filename_temp_agg_long <- paste0("aggregation_output_", temp_agg_args$out_temp_res,  "_long.parquet")
    write_parquet(temp_agg_long, file.path(paths$path_out_folder, filename_temp_agg_long), compression = compression)
    message("Temporally aggregated output in long format saved to: aggregation_output_",temp_agg_args$out_temp_res, "_long.parquet")
  }
  
  # Save long daily output
  if (temp_agg_args$out_temp_res != "daily" && save_long_daily_output) {
    message("\n---------------------------------------------------",
            "\n   5.2 Reshape Daily Output to Long Format",
            "\n---------------------------------------------------\n")
    
    spatial_agg_long <- reshape_to_long(spatial_agg, add_time_columns = TRUE)
    spatial_agg_long <- add_appended_cols(spatial_agg_long, polygons, poly_id_col, appended_cols)
    
    if (output_path_provided) {
      filename_spatial_agg_long <- paste0("aggregation_output_daily_long.parquet")
      write_parquet(spatial_agg_long, file.path(paths$path_out_folder, filename_spatial_agg_long), compression = compression)
      message("Daily aggregation output in long format saved to: aggregation_output_daily_long.parquet")
    }
  }
  
  # Collect the relevant data frames
  results <- list(
    spatial_agg = spatial_agg,
    spatial_agg_long = if (temp_agg_args$out_temp_res != "daily" && save_long_daily_output) spatial_agg_long else NULL,
    temp_agg_wide = if (temp_agg_args$out_temp_res != "daily") temp_agg_wide else NULL,
    temp_agg_long = temp_agg_long,
    area_weights = if (save_area_weights) area_weights else NULL
  )
  
  # ---- 6. Control Checks ----------------------------------------------------------
  
  message("\n===========================================================================",
          "\n 6. Control Checks",
          "\n===========================================================================\n")
  
  # Check if control checks should be run
  if (!run_control_checks) {
    if (!output_path_provided) {
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