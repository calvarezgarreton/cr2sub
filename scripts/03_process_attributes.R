# -----------------------------------------------------------------------------
# Script: 03_process_attributes.R
# Purpose: Derive CR2SUB well attributes by merging metadata, terrain rasters,
#   and hydroclimate summaries to produce curated attributes.
# Author: Rodrigo Marinao Rivas, Camila Alvarez Garretón
# Date: [2025-09-16]
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section: Dependencies
# -----------------------------------------------------------------------------

library(terra)
library(zoo)

source("scripts/functions/functions_pipeline_config.R")
source("scripts/functions/functions_process_cr2sub_attributes.R")

# -----------------------------------------------------------------------------
# Section: Configuration
# -----------------------------------------------------------------------------

tag <- pipeline_config("tag")
version <- pipeline_config("version")
project_tmp_dir <- pipeline_config("tmp_dir")
out_folder <- pipeline_config("output_dir")
ini_mean_values <- pipeline_config_integer("ini_mean_values")
end_mean_values <- pipeline_config_integer("end_mean_values")

if (ini_mean_values > end_mean_values) {
  stop("ini_mean_values debe ser menor o igual que end_mean_values.")
}
period_column_name <- function(prefix) {
  paste0(prefix, ini_mean_values, "_", end_mean_values)
}

# -----------------------------------------------------------------------------
# Section: Input Files
# -----------------------------------------------------------------------------

metadata_file <- pipeline_config("metadata_file")

if (!file.exists(metadata_file)) {
  stop(
    "No se encontró el archivo de metadata esperado en ",
    dirname(metadata_file), "."
  )
}

dem_file <- pipeline_config("dem_file")
slope_file <- pipeline_config("slope_file")

pr_file <- pipeline_config("cr2met_pr_file")
et0_file <- pipeline_config("cr2met_et0_file")
sf_file <- pipeline_config("cr2met_sf_file")
tmean_file <- pipeline_config("cr2met_tmean_file")

clsoil_file <- pipeline_config("clsoil_data_file")

basin_camels <- vect(pipeline_config("camels_basins_file"))
basin_bna <- vect(pipeline_config("bna_basins_file")) |>
  project(y = "epsg:4326")

join_table_camels_basins <- pipeline_config("camels_join_file")
join_table_bna_basins <- pipeline_config("bna_join_file")

gwl_ts_file <- pipeline_config("gwl_ts_file")
gwl_ts_clean_file <- pipeline_config("gwl_ts_clean_file")
output_file <- pipeline_config("attributes_file")

# -----------------------------------------------------------------------------
# Section: Load Metadata And Time Series
# -----------------------------------------------------------------------------

metadata_raw <- read.csv(metadata_file, stringsAsFactors = FALSE)
metadata_df <- data.frame(
  cr2sub_id = metadata_raw$cr2sub_id,
  dga_well_code = metadata_raw$dga_well_code,
  dga_well_name = metadata_raw$dga_well_name,
  dga_well_basin = metadata_raw$dga_well_basin,
  dga_well_subbasin = metadata_raw$dga_well_subbasin,
  dga_well_lat = as.numeric(metadata_raw$dga_well_lat),
  dga_well_lon = as.numeric(metadata_raw$dga_well_lon),
  dga_well_utm_north = as.numeric(metadata_raw$dga_well_utm_north),
  dga_well_utm_east = as.numeric(metadata_raw$dga_well_utm_east),
  dga_well_elev = as.numeric(metadata_raw$dga_well_elev)
)

# -----------------------------------------------------------------------------
# Section: Derive Core Spatial Attributes
# -----------------------------------------------------------------------------

metadata_df$cr2sub_well_lat <- metadata_df$dga_well_lat
metadata_df$cr2sub_well_lon <- metadata_df$dga_well_lon
metadata_df$cr2sub_well_utm_north_h19 <- metadata_df$dga_well_utm_north
metadata_df$cr2sub_well_utm_east_h19 <- metadata_df$dga_well_utm_east

well_vect <- vect(
  metadata_df,
  geom = c("dga_well_lon", "dga_well_lat"),
  crs = "EPSG:4326",
  keepgeom = TRUE
)

# -----------------------------------------------------------------------------
# Section: Raster Resources And Helpers
# -----------------------------------------------------------------------------

dem_raster <- terra::rast(dem_file)
slope_raster <- terra::rast(slope_file)

filter_annual_range <- function(rast_obj, start_year, end_year) {
  time_vals <- terra::time(rast_obj)
  if (is.null(time_vals)) {
    stop("Raster lacks temporal metadata for filtering.")
  }

  if (inherits(time_vals, "Date") || inherits(time_vals, "POSIXt")) {
    years <- as.integer(format(time_vals, "%Y"))
  } else {
    years <- as.integer(time_vals)
  }

  keep <- !is.na(years) & years >= start_year & years <= end_year
  if (!any(keep)) {
    stop("No raster layers within requested year range.")
  }

  rast_obj[[which(keep)]]
}

mean_over_years <- function(file_path, start_year, end_year) {
  rast_obj <- terra::rast(file_path)
  yearly_subset <- filter_annual_range(rast_obj, start_year, end_year)
  terra::app(yearly_subset, mean, na.rm = TRUE)
}

mean_annual_basin_ratio <- function(
  numerator_file,
  denominator_file,
  basins,
  basin_ids,
  start_year,
  end_year
) {
  numerator <- filter_annual_range(
    terra::rast(numerator_file), start_year, end_year
  )
  denominator <- filter_annual_range(
    terra::rast(denominator_file), start_year, end_year
  )

  if (!identical(
    as.character(terra::time(numerator)),
    as.character(terra::time(denominator))
  )) {
    stop("Annual ratio rasters must have matching years.")
  }

  numerator_by_year <- terra::extract(
    numerator, basins, fun = mean, na.rm = TRUE
  )[, -1, drop = FALSE] |> as.matrix()
  denominator_by_year <- terra::extract(
    denominator, basins, fun = mean, na.rm = TRUE
  )[, -1, drop = FALSE] |> as.matrix()

  annual_ratios <- numerator_by_year / denominator_by_year
  annual_ratios[!is.finite(annual_ratios)] <- NA_real_
  basin_means <- apply(annual_ratios, 1, function(ratios) {
    if (!any(!is.na(ratios))) {
      return(NA_real_)
    }
    mean(ratios, na.rm = TRUE)
  })

  stats::setNames(basin_means, as.character(basin_ids))
}

pr_mean_rast <- mean_over_years(pr_file, ini_mean_values, end_mean_values)
et0_mean_rast <- mean_over_years(et0_file, ini_mean_values, end_mean_values)
tmean_mean_rast <- mean_over_years(tmean_file, ini_mean_values, end_mean_values)

# -----------------------------------------------------------------------------
# Section: Topographic Attributes At Wells
# -----------------------------------------------------------------------------

dem_vals <- terra::extract(dem_raster, well_vect)[, 2]
slope_vals <- terra::extract(slope_raster, well_vect)[, 2]

metadata_df$cr2sub_well_elev <- round(dem_vals, 2)
metadata_df$cr2sub_well_slp <- round(slope_vals, 2)

# -----------------------------------------------------------------------------
# Section: Groundwater Time-Series Statistics
# -----------------------------------------------------------------------------

gwl_ts <- read.csv.zoo(gwl_ts_file, check.names = FALSE)

gwl_ts_clean <- read.csv.zoo(gwl_ts_clean_file, check.names = FALSE)

gwl_years <- as.integer(format(zoo::index(gwl_ts), "%Y"))
gwl_period <- !is.na(gwl_years) &
  gwl_years >= ini_mean_values &
  gwl_years <= end_mean_values
gwl_ts_period <- gwl_ts[gwl_period, ]

gwl_clean_years <- as.integer(format(zoo::index(gwl_ts_clean), "%Y"))
gwl_clean_period <- !is.na(gwl_clean_years) &
  gwl_clean_years >= ini_mean_values &
  gwl_clean_years <= end_mean_values
gwl_ts_clean_period <- gwl_ts_clean[gwl_clean_period, ]

calculate_gwl_stat <- function(gwl_period, stat_fun) {
  values <- zoo::coredata(gwl_period)
  if (nrow(values) == 0) {
    return(setNames(rep(NA_real_, ncol(values)), colnames(values)))
  }

  apply(values, 2, function(well_values) {
    if (!any(!is.na(well_values))) {
      return(NA_real_)
    }
    stat_fun(well_values, na.rm = TRUE)
  })
}

well_ids <- as.character(metadata_df$cr2sub_id)
mean_gwl <- calculate_gwl_stat(gwl_ts_period, mean)[well_ids]
sd_gwl <- calculate_gwl_stat(gwl_ts_period, sd)[well_ids]
mean_gwl_clean <- calculate_gwl_stat(gwl_ts_clean_period, mean)[well_ids]
sd_gwl_clean <- calculate_gwl_stat(gwl_ts_clean_period, sd)[well_ids]

metadata_df[[period_column_name("cr2sub_well_gwl_mon_mean_")]] <-
  round(mean_gwl, 2)
metadata_df[[period_column_name("cr2sub_well_gwl_mon_sd_")]] <-
  round(sd_gwl, 2)
metadata_df[[period_column_name("cr2sub_well_gwl_mon_cv_")]] <-
  round(sd_gwl / abs(mean_gwl), 2)
metadata_df[[period_column_name("cr2sub_well_gwl_mon_clean_mean_")]] <-
  round(mean_gwl_clean, 2)
metadata_df[[period_column_name("cr2sub_well_gwl_mon_clean_sd_")]] <-
  round(sd_gwl_clean, 2)
metadata_df[[period_column_name("cr2sub_well_gwl_mon_clean_cv_")]] <-
  round(sd_gwl_clean / abs(mean_gwl_clean), 2)

# -----------------------------------------------------------------------------
# Section: Basin Attributes – CAMELS-CL
# -----------------------------------------------------------------------------

join_table_camels <- read.csv(join_table_camels_basins)

camels_ids <- suppressWarnings(sapply(
  metadata_raw$cr2sub_id,
  select_basin_from_join_table,
  join_table = join_table_camels,
  fun = min
))

metadata_df$cr2sub_well_in_basin_camels <- camels_ids

camels_pr_vals <- terra::extract(pr_mean_rast, basin_camels,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(camels_pr_vals) <-
  as.data.frame(basin_camels)[, "gauge_id"] |> as.character()

camels_tmean_vals <- terra::extract(tmean_mean_rast, basin_camels,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(camels_tmean_vals) <-
  as.data.frame(basin_camels)[, "gauge_id"] |> as.character()

camels_et0_vals <-
  terra::extract(et0_mean_rast, basin_camels,
    fun = mean, na.rm = TRUE
  )[, 2] |> as.numeric()

names(camels_et0_vals) <-
  as.data.frame(basin_camels)[, "gauge_id"] |> as.character()

camels_elev_vals <- terra::extract(dem_raster,
  basin_camels,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(camels_elev_vals) <-
  as.data.frame(basin_camels)[, "gauge_id"] |> as.character()

camels_slope_vals <- terra::extract(slope_raster,
  basin_camels,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(camels_slope_vals) <-
  as.data.frame(basin_camels)[, "gauge_id"] |> as.character()

metadata_df[[period_column_name("cr2sub_camels_pr_yr_mean_")]] <-
  round(camels_pr_vals[camels_ids], 2)
metadata_df[[period_column_name("cr2sub_camels_tmean_yr_mean_")]] <-
  round(camels_tmean_vals[camels_ids], 2)
metadata_df[[period_column_name("cr2sub_camels_et0_yr_mean_")]] <-
  round(camels_et0_vals[camels_ids], 2)
camels_basin_ids <- as.data.frame(basin_camels)[, "gauge_id"]
camels_snow_fraction_vals <- mean_annual_basin_ratio(
  sf_file, pr_file, basin_camels, camels_basin_ids,
  ini_mean_values, end_mean_values
)
camels_aridity_vals <- mean_annual_basin_ratio(
  et0_file, pr_file, basin_camels, camels_basin_ids,
  ini_mean_values, end_mean_values
)
metadata_df[[period_column_name("cr2sub_camels_sf_")]] <-
  round(camels_snow_fraction_vals[camels_ids], 2)
metadata_df[[period_column_name("cr2sub_camels_aridity_")]] <-
  round(camels_aridity_vals[camels_ids], 2)
metadata_df$cr2sub_camels_elev <- round(camels_elev_vals[camels_ids], 2)
metadata_df$cr2sub_camels_slp <- round(camels_slope_vals[camels_ids], 2)

# -----------------------------------------------------------------------------
# Section: Basin Attributes – BNA
# -----------------------------------------------------------------------------

join_table_bna <- read.csv(join_table_bna_basins)

bna_ids <- suppressWarnings(sapply(
  metadata_raw$cr2sub_id,
  select_basin_from_join_table,
  join_table = join_table_bna,
  fun = min
))

metadata_df$cr2sub_well_in_basin_bna <- bna_ids

bna_pr_vals <- terra::extract(pr_mean_rast,
  basin_bna,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()
names(bna_pr_vals) <-
  as.data.frame(basin_bna)[, "COD_CUEN"] |> as.numeric() |> as.character()

bna_tmean_vals <- terra::extract(tmean_mean_rast, basin_bna,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(bna_tmean_vals) <-
  as.data.frame(basin_bna)[, "COD_CUEN"] |> as.numeric() |> as.character()

bna_et0_vals <- terra::extract(et0_mean_rast,
  basin_bna,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(bna_et0_vals) <-
  as.data.frame(basin_bna)[, "COD_CUEN"] |>
  as.numeric()

bna_elev_vals <- terra::extract(dem_raster,
  basin_bna,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(bna_elev_vals) <-
  as.data.frame(basin_bna)[, "COD_CUEN"] |> as.numeric()

bna_slope_vals <- terra::extract(slope_raster,
  basin_bna,
  fun = mean, na.rm = TRUE
)[, 2] |> as.numeric()

names(bna_slope_vals) <-
  as.data.frame(basin_bna)[, "COD_CUEN"] |> as.numeric()

metadata_df[[period_column_name("cr2sub_bna_pr_yr_mean_")]] <-
  round(bna_pr_vals[bna_ids], 2)
metadata_df[[period_column_name("cr2sub_bna_tmean_yr_mean_")]] <-
  round(bna_tmean_vals[bna_ids], 2)
metadata_df[[period_column_name("cr2sub_bna_et0_yr_mean_")]] <-
  round(bna_et0_vals[bna_ids], 2)
bna_basin_ids <- as.data.frame(basin_bna)[, "COD_CUEN"]
bna_snow_fraction_vals <- mean_annual_basin_ratio(
  sf_file, pr_file, basin_bna, bna_basin_ids,
  ini_mean_values, end_mean_values
)
bna_aridity_vals <- mean_annual_basin_ratio(
  et0_file, pr_file, basin_bna, bna_basin_ids,
  ini_mean_values, end_mean_values
)
metadata_df[[period_column_name("cr2sub_bna_sf_")]] <-
  round(bna_snow_fraction_vals[bna_ids], 2)
metadata_df[[period_column_name("cr2sub_bna_aridity_")]] <-
  round(bna_aridity_vals[bna_ids], 2)
metadata_df$cr2sub_bna_elev <- round(bna_elev_vals[bna_ids], 2)
metadata_df$cr2sub_bna_slp <- round(bna_slope_vals[bna_ids], 2)

# -----------------------------------------------------------------------------
# Section: Soil Attributes – CLSoilMaps
# -----------------------------------------------------------------------------

clsmap <- read.csv(
  clsoil_file,
  colClasses = c(cr2sub_id = "character"),
  row.names = 1
)
clsmap <- clsmap[sprintf("%08d", as.integer(metadata_df$cr2sub_id)), ]

metadata_df <- cbind.data.frame(metadata_df, clsmap)

# -----------------------------------------------------------------------------
# Section: Consolidation And Export
# -----------------------------------------------------------------------------

required_order <- c(
  "cr2sub_id",
  "dga_well_code",
  "dga_well_name",
  "dga_well_basin",
  "dga_well_subbasin",
  "dga_well_lat",
  "dga_well_lon",
  "dga_well_utm_north",
  "dga_well_utm_east",
  "dga_well_elev",
  "cr2sub_well_lat",
  "cr2sub_well_lon",
  "cr2sub_well_utm_north_h19",
  "cr2sub_well_utm_east_h19",
  "cr2sub_well_elev",
  "cr2sub_well_slp",
  period_column_name("cr2sub_well_gwl_mon_mean_"),
  period_column_name("cr2sub_well_gwl_mon_sd_"),
  period_column_name("cr2sub_well_gwl_mon_cv_"),
  period_column_name("cr2sub_well_gwl_mon_clean_mean_"),
  period_column_name("cr2sub_well_gwl_mon_clean_sd_"),
  period_column_name("cr2sub_well_gwl_mon_clean_cv_"),
  "cr2sub_well_in_basin_camels",
  period_column_name("cr2sub_camels_pr_yr_mean_"),
  period_column_name("cr2sub_camels_tmean_yr_mean_"),
  period_column_name("cr2sub_camels_et0_yr_mean_"),
  period_column_name("cr2sub_camels_aridity_"),
  period_column_name("cr2sub_camels_sf_"),
  "cr2sub_camels_elev",
  "cr2sub_camels_slp",
  "cr2sub_well_in_basin_bna",
  period_column_name("cr2sub_bna_pr_yr_mean_"),
  period_column_name("cr2sub_bna_tmean_yr_mean_"),
  period_column_name("cr2sub_bna_et0_yr_mean_"),
  period_column_name("cr2sub_bna_aridity_"),
  period_column_name("cr2sub_bna_sf_"),
  "cr2sub_bna_elev",
  "cr2sub_bna_slp",
  "cr2sub_clsoilmap_awc_0_100cm",
  "cr2sub_clsoilmap_awc_100_200cm",
  "cr2sub_clsoilmap_bulkd_0_100cm",
  "cr2sub_clsoilmap_bulkd_100_200cm",
  "cr2sub_clsoilmap_clay_0_100cm",
  "cr2sub_clsoilmap_clay_100_200cm",
  "cr2sub_clsoilmap_ksat_0_100cm",
  "cr2sub_clsoilmap_ksat_100_200cm",
  "cr2sub_clsoilmap_sand_0_100cm",
  "cr2sub_clsoilmap_sand_100_200cm"
)

cr2sub_attributes_df <- metadata_df[, required_order]

write.csv(cr2sub_attributes_df, output_file, row.names = FALSE)

message("Atributos exportados en ", output_file)
message(
  "Dimensiones finales: ", nrow(cr2sub_attributes_df), " pozos x ",
  ncol(cr2sub_attributes_df), " atributos"
)

if (file.exists(output_file)) {
  file_info <- file.info(output_file)
  if (!is.na(file_info$size) && file_info$size > 0 && dir.exists(project_tmp_dir)) {
    unlink(project_tmp_dir, recursive = TRUE)
  }
}
