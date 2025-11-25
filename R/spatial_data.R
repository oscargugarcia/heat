#' Get Valid Raster Files
#'
#' @Description:
#'   This function searches for raster files with specific extensions in a given folder,
#'   checking whether they can be read into R. It returns a list of accessible raster files 
#'   while printing warnings for any files that cannot be read.
#'
#' @Input:
#'   - path: A character string specifying the path to the folder containing raster files.
#'   - pattern: A character string specifying the file pattern to match (default is "\\.tif$|\\.nc$").
#'
#' @Output:
#'   - A character vector containing paths to the valid raster files that can be read.
#'
#' @Examples:
#'   valid_rasters <- get_valid_raster_files(paths$path_env_rast)
get_valid_raster_files <- function(path, pattern = "\\.tif$|\\.nc$") {
  # Identify all raster files in the folder with full paths
  rast_files <- list.files(
    path = path,
    pattern = pattern,
    full.names = FALSE,
    recursive = FALSE
  )
  
  # Initialize a vector to hold paths of readable files
  available_files <- character(0)
  
  # Check each file to confirm it can be read
  for (file in rast_files) {
    file_path <- file.path(path, file)  # Construct full file path
    if (file.exists(file_path)) {
      # Attempt to read the file and catch any errors
      try({
        # Attempt to create a SpatRaster object
        test_raster <- terra::rast(file_path)
        available_files <- c(available_files, file) # Add to available list if successful
      }, silent = TRUE)  # suppress error messages
      
    } else {
      warning("File does not exist or is not accessible: ", file)  # Print warning if file does not exist
    }
  }
  
  unavailable_files <- setdiff(rast_files, available_files)
  if (length(unavailable_files) > 0) {
    message("The following files could not be read. If the raster is located in Dropbox, make sure it is available offline before executing the script. Unavailable files: \n", paste(unavailable_files, collapse = ", "))
  }
  
  return(available_files)
}



#' Filter and Read Rasters
#'
#' @description
#' This function loads files from a specified directory. It allows users to filter the files based on a year range, with the years inferred from the filenames, which are expected to be in the format 'yyyy.tif' or 'yyyy.nc'. If no filtering range is provided, all available raster files will be read. The function prints informative messages detailing the loading process, including the path and names of the files being read.
#'
#' @param path A character string specifying the directory path containing raster files.
#' @param sec_weight_range A numeric vector of length 2 (optional). This vector specifies the year range to filter the raster files. If NULL (default), all files are loaded without filtering.
#'
#' @return A character vector of file paths corresponding to the raster files that were read. The files are loaded into memory as a list of raster objects.
#'
#' @examples
#' \dontrun{
#'   # Load all secondary weight rasters
#'   all_rasters <- filter_read_sec_weights("path/to/rasters/")
#'   
#'   # Load secondary weight rasters filtered by year range
#'   filtered_rasters <- filter_read_sec_weights("path/to/rasters/", c(1970, 2010))
#' }
filter_read_rast <- function(path, year_range = NULL) {
  
  # Identify all available rasters in the folder
  rast_files <- get_valid_raster_files(path, pattern = "\\.tif$|\\.nc$") 
  
  # Print the path from which files are being loaded
  message(sprintf("Loading files from path: %s", path))

  if (is.null(year_range)) {
    # No specific range provided, read all files
    filtered_files <- rast_files
    message("No filtering applied. Reading in all available files.")
  } else {
    
    # Ensure dates in year_range are correctly formatted
    
    
    # Extract years from the file names
    years <- as.numeric(gsub("\\.tif$|\\.nc$", "", basename(rast_files)))

    # Filter files based on the specified year range
    filtered_files <- rast_files[years >= year_range[1] & years <= year_range[2]]
    message(sprintf("Filtering files based on the date range provided: %d to %d", year_range[1], year_range[2]))
  }
  
  # Print the names of the files being read
  message("Files being read:")
  # print filenames separated by commas
  message(paste(basename(filtered_files), collapse = ", "))
  
  return(filtered_files)
}


#' Read Spatial File (Polygons or Points)
#'
#' @description
#' Reads spatial data from a file (Parquet, GPKG, Shapefile, GeoJSON, FGB, RDS), converts the geometry 
#' to an sf object, and ensures the CRS is set to WGS84 (EPSG:4326). If the geometry column is named "geom", 
#' it is renamed to "geometry". For point geometries, the points are retained as-is and will be processed
#' using optimized point extraction methods. For polygon geometries, st_make_valid() is applied to ensure 
#' validity. Additionally, the function checks for empty or invalid geometries and warns the user if any 
#' are found, listing the corresponding row numbers. Empty geometries are removed.
#'
#' @param file_path A string specifying the path to the spatial data file.
#'
#' @return An sf object with geometry in WGS84 (EPSG:4326) and a geometry column named "geometry".
#'   Point and polygon geometries are both supported and processed accordingly.
#'
#' @examples
#' \dontrun{
#'   polygons <- read_spatial_file("data/polygons/CONUS_consistent_counties.gpkg")
#'   points <- read_spatial_file("data/points/sampling_points.gpkg")
#'   rds_data <- read_spatial_file("data/spatial_data.rds")
#' }
read_spatial_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Error: File does not exist - ", file_path)
  }
  
  # Extract file extension (lowercase)
  file_ext <- tolower(tools::file_ext(file_path))
  
  # Read based on format
  if (file_ext %in% c("parquet", "pq")) {
    message("Reading Parquet file: ", file_path)
    polygons_df <- arrow::read_parquet(file_path)
    
    if (!("geometry" %in% names(polygons_df)) && !("geom" %in% names(polygons_df))) {
      stop("Error: No geometry column found in Parquet file")
    }
    
    if ("geometry" %in% names(polygons_df)) {
      polygons_df$geometry <- st_as_sfc(polygons_df$geometry)
    } else {
      polygons_df$geom <- st_as_sfc(polygons_df$geom)
    }
    
    polygons <- st_as_sf(polygons_df, crs = "WGS84")
    
  } else if (file_ext == "rds") {
    message("Reading RDS file: ", file_path)
    polygons <- readRDS(file_path)
    
    # Check if the loaded object is an sf object
    if (!inherits(polygons, "sf")) {
      stop("Error: RDS file does not contain an sf object")
    }
    
  } else if (file_ext %in% c("gpkg", "shp", "geojson", "json", "fgb")) {
    message("Reading spatial file: ", file_path)
    polygons <- st_read(file_path, quiet = TRUE)
    
  } else {
    stop("Error: Unsupported file format - ", file_ext)
  }
  
  # # Create a row_id column
  # polygons$row_id <- 1:nrow(polygons) 
  # polygons <- relocate(polygons, row_id)
  
  # Check the original CRS
  original_crs <- st_crs(polygons)
  original_epsg <- st_crs(polygons)$epsg
  
  # If CRS is NA, assume it's WGS 84 (EPSG:4326)
  if (is.na(original_crs)) {
    warning("The CRS is missing in spatial file ", file_path, "\n  Setting the CRS to WGS 84 (EPSG:4326).")
    st_crs(polygons) <- 4326  # Set the original CRS to EPSG:4326 if unknown
  }
  
  # Check current CRS and transform if necessary
  if (st_crs(polygons) != 4326) {
    message("Transforming EPSG from ", original_epsg, " to 4326 (WGS84).")
    # Use st_transform to convert to WGS84
    polygons <- suppressMessages(st_transform(polygons, crs = 4326, quiet = TRUE))
  } else {
    message("CRS is already WGS 84 (EPSG:4326). No transformation needed.")
  }

  
  # Rename geometry column 'geom' to 'geometry' if needed.
  if ("geom" %in% names(polygons)) {
    message("Renaming geometry column 'geom' to 'geometry'.")
    names(polygons)[names(polygons) == "geom"] <- "geometry"
    st_geometry(polygons) <- "geometry"
  }
  
  # Check for empty geometries first (before any processing)
  empty_idx <- which(st_is_empty(polygons))
  
  if (length(empty_idx) > 0) {
    warning("The following ", length(empty_idx), "/", nrow(polygons), " rows have empty geometries and will be removed:\n      ",
            paste(empty_idx, collapse = ", "), immediate. = TRUE)
    
    # Remove empty geometries
    polygons <- polygons[-empty_idx, ]
  }
  
  # Check geometry types and handle mixed geometries
  geom_types <- st_geometry_type(polygons)
  unique_geom_types <- unique(geom_types)
  
  # Check if we have mixed geometry types
  has_points <- any(grepl("POINT", unique_geom_types))
  has_polygons <- any(grepl("POLYGON", unique_geom_types))
  
  if (has_points && has_polygons) {
    warning("Mixed geometry types detected (both points and polygons). This may cause issues in downstream analysis.")
    message("Geometry types found: ", paste(unique_geom_types, collapse = ", "))
  }
  
  # Report point geometries
  if (has_points) {
    point_indices <- grepl("POINT", geom_types)
    n_points <- sum(point_indices)
    message("Point geometries detected: ", n_points, " point(s) will be processed using optimized point extraction.")
  }
  
  # Handle polygon geometries (apply st_make_valid)
  if (has_polygons) {
    polygon_indices <- grepl("POLYGON", geom_types)
    n_polygons <- sum(polygon_indices)
    
    # Check for invalid polygons first
    invalid_polygon_idx <- polygon_indices & !st_is_valid(polygons)
    n_invalid_polygons <- sum(invalid_polygon_idx)
    
    if (n_invalid_polygons > 0) {
      invalid_rows <- which(invalid_polygon_idx)
      warning("The following ", n_invalid_polygons, "/", n_polygons, " polygon rows have invalid geometries that will be fixed using st_make_valid():\n      ",
              paste(invalid_rows, collapse = ", "), immediate. = TRUE)
    }
    
    # Apply st_make_valid() only to polygon geometries
    message("Applying st_make_valid() to ", n_polygons, " polygon geometries to ensure validity.")
    
    # Create a copy to work with
    temp_polygons <- polygons
    temp_polygons[polygon_indices, ] <- st_make_valid(temp_polygons[polygon_indices, ])
    
    # Check if any polygons became empty after st_make_valid()
    post_valid_empty_idx <- which(polygon_indices & st_is_empty(temp_polygons))
    
    if (length(post_valid_empty_idx) > 0) {
      warning("The following ", length(post_valid_empty_idx), "/", n_polygons, " polygon rows became empty after st_make_valid() and will be removed:\n      ",
              paste(post_valid_empty_idx, collapse = ", "), immediate. = TRUE)
      
      # Remove geometries that became empty after validation
      temp_polygons <- temp_polygons[-post_valid_empty_idx, ]
    }
    
    polygons <- temp_polygons
    
    # Handle geometries crossing the dateline
    polygons <- st_wrap_dateline(polygons, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  }
  
  # Final check for any remaining empty geometries
  final_empty_idx <- which(st_is_empty(polygons))
  
  if (length(final_empty_idx) > 0) {
    warning("Final check: The following ", length(final_empty_idx), "/", nrow(polygons), " rows have empty geometries and will be removed:\n      ",
            paste(final_empty_idx, collapse = ", "), immediate. = TRUE)
    
    # Remove empty geometries
    polygons <- polygons[-final_empty_idx, ]
  }
  
  # Final validation check
  final_invalid_idx <- which(!st_is_valid(polygons))
  
  if (length(final_invalid_idx) > 0) {
    warning("Final check: The following ", length(final_invalid_idx), "/", nrow(polygons), " rows still have invalid geometries:\n      ",
            paste(final_invalid_idx, collapse = ", "), immediate. = TRUE)
  }
  
  message("Final result: ", nrow(polygons), " valid geometries loaded.")
  
  return(polygons)
}

#' Buffer Polygons Based on Raster Resolution
#'
#' @Description:
#'   This function computes a buffered extent around the input polygons using the raster's resolution.
#'   The extent is expanded by a factor specified by the user (buffer_factor). The buffered extent is
#'   returned as a terra::ext object.
#'
#' @Input:
#'   - raster: (SpatRaster) A raster object from which the resolution is extracted.
#'   - polygons: (sf or SpatVector) A spatial object containing the polygons.
#'   - buffer_factor: (numeric) A factor by which to multiply the raster resolution when buffering
#'         the polygons' extent.
#'
#' @Output:
#'   - A buffered extent (terra::ext object) computed as:
#'         xmin = polygons_extent$xmin - buffer_factor * x_res,
#'         xmax = polygons_extent$xmax + buffer_factor * x_res,
#'         ymin = polygons_extent$ymin - buffer_factor * y_res,
#'         ymax = polygons_extent$ymax + buffer_factor * y_res.
#'
#' @Examples:
#'   buffered_extent <- buffer_polygons(environmental_raster, polygons, buffer_factor = 1)
buffer_polygons <- function(raster, polygons, buffer_factor) { 
  # Extract raster resolution
  x_res <- res(raster)[1]
  y_res <- res(raster)[2]
  
  # Get polygon bounding box as a terra::ext object.
  polygons_extent <- ext(polygons)
  
  # Create a buffered extent using the buffer_factor.
  buffered_extent <- ext(c(
    xmin = polygons_extent$xmin - buffer_factor * x_res,
    xmax = polygons_extent$xmax + buffer_factor * x_res,
    ymin = polygons_extent$ymin - buffer_factor * y_res,
    ymax = polygons_extent$ymax + buffer_factor * y_res
  ))
  
  return(buffered_extent)
}

#' Filter Raster by Date Range
#'
#' @Description:
#'   This function filters a multi-layer raster based on a specified date range. The raster should have
#'   layer names corresponding to date strings in 'YYYY-MM-DD' format. The function checks whether the
#'   specified date range falls within the range of available raster dates and returns the filtered
#'   raster containing only the layers that match the specified dates.
#'
#' @Input:
#'   - env_rast: A multi-layer raster object (e.g., from the 'raster' or 'terra' package) with layer names
#'     corresponding to valid date strings (e.g., 2017-01-01).
#'   - start_date: A string or Date object representing the start date of the desired date range.
#'   - end_date: A string or Date object representing the end date of the desired date range.
#'
#' @Output:
#'   - A filtered raster object containing only the layers corresponding to the specified date range.
#'
#' @Examples:
#'   # Assuming 'my_raster' is a multi-layer raster and you wish to filter
#'   # it for dates between January 1, 2018, and December 31, 2018:
#'   filtered_raster <- filter_env_rast(my_raster, start_date = "2018-01-01", end_date = "2018-12-31")
filter_env_rast <- function(env_rast, start_date, end_date) {
  # Check if raster layer names exist
  layer_names <- names(env_rast)
  if (is.null(layer_names)) {
    stop("Error: Raster layers have no names. Ensure layers are named with valid date strings (e.g., 'YYYY-MM-DD').")
  }
  
  # Try converting layer names to Date format and handle errors
  raster_dates <- tryCatch(
    as.Date(layer_names),
    error = function(e) {
      stop("Error: Raster layer names could not be converted to valid dates. ",
           "Ensure they are in a recognizable date format such as 'YYYY-MM-DD'.")
    }
  )
  
  # Ensure that at least some raster layer names were successfully converted
  if (all(is.na(raster_dates))) {
    stop("Error: All raster layer names are NA after attempting date conversion. ",
         "Please check the naming format.")
  }
  
  # Get the overall raster date range
  raster_start <- min(raster_dates, na.rm = TRUE)
  raster_end <- max(raster_dates, na.rm = TRUE)
  
  # Convert input dates to Date format
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  # Check if period definitions exceed raster date range
  if (start_date < raster_start || end_date > raster_end) {
    stop("Error: The specified date range exceeds the raster date range.\n",
         "   Raster covers: ", raster_start, " to ", raster_end, "\n",
         "   Periods cover: ", start_date, " to ", end_date, "\n",
         "   Ensure that period definitions are within the available raster time range.")
  }
  
  if (start_date > raster_start || end_date < raster_end) {
    message("Selected period date range is shorter than raster date range:\n",
            "   Raster covers: ", raster_start, " to ", raster_end, "\n",
            "   Periods cover: ", start_date, " to ", end_date, "\n",
            "Filtering to period date range...")
  }
  
  
  # Filter the raster based on the specified date range
  date_mask <- raster_dates >= start_date & raster_dates <= end_date
  filtered_raster <- env_rast[[date_mask]]
  message("Filtered Raster date range: ", start_date, " to ", end_date, "\n")
  
  
  return(filtered_raster)
}


#' Parse Boundary Dates from Raster Layer Names
#'
#' @Description:
#'   Determines temporal resolution and converts first/last layer names to date boundaries.
#'   For hourly/daily: returns dates as-is
#'   For monthly: returns first day of first month and last day of last month
#'   For yearly: returns Jan 1st of first year and Dec 31st of last year
#'
#' @Input:
#'   - env_rast: Terra raster with date-formatted layer names
#'
#' @Output:
#'   - Vector of two Date objects representing start and end boundaries
#'
parse_boundary_dates <- function(env_rast) {
  
  # Get layer names and determine temporal resolution with fallback
  layer_names <- names(env_rast)  
  temp_res <- get_temp_res(layer_names)
  
  first_name <- layer_names[1]
  last_name <- layer_names[length(layer_names)]
  
  # Convert based on resolution
  if (temp_res == "hourly" || temp_res == "daily") {
    # Extract just the date part (YYYY-MM-DD) from YYYY-MM-DD HH:MM
    first_date <- as.Date(first_name)
    last_date <- as.Date(last_name)
    
  } else if (temp_res == "monthly") {
    # Convert YYYY-MM to first day of first month and last day of last month
    first_date <- as.Date(paste0(first_name, "-01"))
    
    last_month_first <- as.Date(paste0(last_name, "-01"))
    days_in_last_month <- lubridate::days_in_month(last_month_first)
    last_date <- as.Date(paste0(last_name, "-", sprintf("%02d", days_in_last_month)))
    
  } else if (temp_res == "yearly") {
    # Convert YYYY to Jan 1st and Dec 31st
    first_date <- as.Date(paste0(first_name, "-01-01"))
    last_date <- as.Date(paste0(last_name, "-12-31"))
    
  } else {
    stop("Unsupported temporal resolution: ", temp_res)
  }
  
  boundary_dates <- c(first_date, last_date)

  return(boundary_dates)
}
