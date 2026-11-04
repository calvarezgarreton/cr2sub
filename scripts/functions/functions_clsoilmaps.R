# -----------------------------------------------------------------------------
# CLSoilMaps extraction helpers.
# Adapted from https://github.com/diegodinamarca/clsoilmap_extractor
# (Diego Dinamarca, MIT). Requires: terra, sf, dplyr, tibble, readr.
# -----------------------------------------------------------------------------


# Shared constants and helpers --------------------------------------------

CLSOILMAP_DEPTH_WEIGHTS <- c(
  "0-5cm" = 5,
  "5-15cm" = 10,
  "15-30cm" = 15,
  "30-60cm" = 30,
  "60-100cm" = 40
)

CLSOILMAP_OUTPUT_PREFIX <- "cr2sub_clsoilmap"

REQUIRED_WELL_ID <- "well_id"

PHYSICAL_PROPERTIES <- c("Clay", "Sand", "Bulkd")
HYDRAULIC_PROPERTIES <- c("AWC", "ksat")

SUBSOIL_DEPTH <- "100-200cm"

stop_if_missing_column <- function(data, column) {
  if (!column %in% names(data)) {
    stop(sprintf("Input wells layer must contain a '%s' column.", column), call. = FALSE)
  }
}

safe_list_tifs <- function(directory) {
  files <- list.files(directory, full.names = TRUE, pattern = "\\.tif$")
  if (length(files) == 0) {
    stop(sprintf("No .tif files found in directory: %s", directory), call. = FALSE)
  }
  files
}

# Input / output ----------------------------------------------------------

read_wells <- function(wells_path) {
  wells <- sf::read_sf(wells_path) |>
    sf::st_transform(4326)
  
  stop_if_missing_column(wells, REQUIRED_WELL_ID)
  wells
}

list_soil_rasters <- function(phys_dir, hydro_dir) {
  list(
    physical = safe_list_tifs(phys_dir),
    hydraulic = safe_list_tifs(hydro_dir)
  )
}

write_output_csv <- function(data, out_path) {
  readr::write_csv(data, out_path)
}

# Raster selection --------------------------------------------------------

select_physical_rasters <- function(files) {
  files[!grepl("Silt", basename(files), ignore.case = TRUE)]
}

select_hydraulic_rasters <- function(files) {
  files <- files[!grepl("Total_AWC", basename(files), ignore.case = TRUE)]
  
  keep <- grepl("AWC", basename(files), ignore.case = TRUE) |
    grepl("ksat", basename(files), ignore.case = TRUE)
  
  files[keep]
}

# Raster extraction -------------------------------------------------------

extract_mean_by_feature <- function(raster_path, wells) {
  r <- terra::rast(raster_path)
  
  extracted <- terra::extract(r, terra::vect(wells))
  
  extracted |>
    dplyr::group_by(ID) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) |>
    dplyr::ungroup() |>
    dplyr::select(-ID)
}

extract_raster_stack_table <- function(raster_files, wells, id_values, verbose = TRUE) {
  out <- tibble::tibble(well_id = id_values)
  
  for (i in seq_along(raster_files)) {
    if (verbose) {
      message(sprintf("%s / %s : %s", i, length(raster_files), basename(raster_files[i])))
    }
    
    values <- extract_mean_by_feature(raster_files[i], wells)
    out <- dplyr::bind_cols(out, values)
  }
  
  out
}

# Aggregation and final formatting ---------------------------------------

build_weighted_mean <- function(data, property, weights = CLSOILMAP_DEPTH_WEIGHTS) {
  layer_names <- names(weights)
  
  column_names <- if (property %in% c("Clay", "Sand", "Bulkd")) {
    paste0(property, ".", layer_names)
  } else {
    paste0(property, "_", layer_names)
  }
  
  missing_cols <- setdiff(column_names, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required columns for weighted mean of '%s': %s",
        property,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  weighted_sum <- Reduce(
    `+`,
    Map(function(col, w) data[[col]] * w, column_names, unname(weights))
  )
  
  weighted_sum / sum(weights)
}

build_awc_sum_0_100 <- function(data) {
  column_names <- paste0("AWC_", names(CLSOILMAP_DEPTH_WEIGHTS))
  
  missing_cols <- setdiff(column_names, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required columns for AWC 0-100cm: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  rowSums(data[, column_names, drop = FALSE], na.rm = FALSE)
}

build_subsoil_column_name <- function(property, depth = SUBSOIL_DEPTH) {
  if (property %in% c("Clay", "Sand", "Bulkd")) {
    paste0(property, ".", depth)
  } else {
    paste0(property, "_", depth)
  }
}

format_output_table <- function(full_df) {
  out <- full_df |>
    dplyr::transmute(
      well_id = .data$well_id,
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_ksat_0_100cm") := build_weighted_mean(full_df, "ksat"),
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_clay_0_100cm") := build_weighted_mean(full_df, "Clay"),
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_sand_0_100cm") := build_weighted_mean(full_df, "Sand"),
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_bulkd_0_100cm") := build_weighted_mean(full_df, "Bulkd"),
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_awc_0_100cm") := build_awc_sum_0_100(full_df),
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_ksat_100_200cm") := .data[[build_subsoil_column_name("ksat")]],
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_clay_100_200cm") := .data[[build_subsoil_column_name("Clay")]],
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_sand_100_200cm") := .data[[build_subsoil_column_name("Sand")]],
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_bulkd_100_200cm") := .data[[build_subsoil_column_name("Bulkd")]],
      !!paste0(CLSOILMAP_OUTPUT_PREFIX, "_awc_100_200cm") := .data[[build_subsoil_column_name("AWC")]]
    ) |>
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ round(.x, 3))
    )
  
  ordered_cols <- c("well_id", sort(setdiff(names(out), "well_id")))
  out[, ordered_cols]
}

#' Extract Soil Attributes from CLSoilMaps for Well Locations
#'
#' @param wells.shp Character. Path to wells vector file (.shp or .geojson).
#' @param phys.dir Character. Directory containing physical soil rasters.
#' @param hydro.dir Character. Directory containing hydraulic soil rasters.
#' @param out Character. Output CSV path.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return Invisibly returns the output data frame and writes a CSV file.
#' @export
extract_clsoilmaps_attr <- function(wells.shp, phys.dir, hydro.dir, out, verbose = TRUE) {
  wells <- read_wells(wells.shp)
  rasters <- list_soil_rasters(phys.dir, hydro.dir)
  
  physical_files <- select_physical_rasters(rasters$physical)
  hydraulic_files <- select_hydraulic_rasters(rasters$hydraulic)
  
  physical_df <- extract_raster_stack_table(
    raster_files = physical_files,
    wells = wells,
    id_values = wells[[REQUIRED_WELL_ID]],
    verbose = verbose
  )
  
  hydraulic_df <- extract_raster_stack_table(
    raster_files = hydraulic_files,
    wells = wells,
    id_values = wells[[REQUIRED_WELL_ID]],
    verbose = verbose
  )
  
  full_df <- physical_df |>
    dplyr::left_join(hydraulic_df, by = "well_id")
  
  result <- format_output_table(full_df)
  
  write_output_csv(result, out)
  
  invisible(result)
}
