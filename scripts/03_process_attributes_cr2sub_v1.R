rm(list = ls())

library(terra)

tag <- "cr2sub"
version <- "v1"

# -----------------------------------------------------------------------------
# archivos de entrada
# -----------------------------------------------------------------------------
metadata_file <- "cr2sub/cr2sub_v1_metadata.csv"
well_ts_file <- "cr2sub/cr2sub_v1_mon_clean.csv"

dem_file <- "input/other_data/dem_masl_fabdemv1.2_2015_300m_epsg4326.tif"
slope_file <- "input/other_data/slope_deg_fabdemv1.2_2015_300m_epsg4326.tif"

pcp_file <- "inpit/other_data/pr_mm_cr2metv2.5_ann_1960_2024_0.05deg_epsg4326.nc"

# -----------------------------------------------------------------------------
# utilidades
# -----------------------------------------------------------------------------
safe_cv <- function(mean_val, sd_val) {
  if (!is.finite(mean_val) || abs(mean_val) < .Machine$double.eps) {
    return(NA_real_)
  }
  sd_val / mean_val
}

clean_stats <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(list(
      mean = NA_real_, sd = NA_real_, cv = NA_real_,
      clean_mean = NA_real_, clean_sd = NA_real_, clean_cv = NA_real_
    ))
  }
  mean_val <- mean(x)
  sd_val <- sd(x)
  cv_val <- safe_cv(mean_val, sd_val)

  if (!is.finite(sd_val) || sd_val == 0) {
    clean <- x
  } else {
    clean <- x[abs(x - mean_val) <= (2 * sd_val)]
    if (!length(clean)) clean <- x
  }

  clean_mean <- mean(clean)
  clean_sd <- sd(clean)
  clean_cv <- safe_cv(clean_mean, clean_sd)

  list(
    mean = mean_val, sd = sd_val, cv = cv_val,
    clean_mean = clean_mean, clean_sd = clean_sd, clean_cv = clean_cv
  )
}

compute_gwl_stats <- function(ts_wide, metadata_ref) {
  well_ids <- intersect(metadata_ref$well_id, names(ts_wide))
  stats_list <- lapply(well_ids, function(wid) {
    clean_stats(ts_wide[[wid]])
  })
  stats_df <- data.frame(
    well_id = well_ids,
    cr2sub_mean_gwl = vapply(stats_list, function(x) x$mean, numeric(1)),
    cr2sub_sd_gwl = vapply(stats_list, function(x) x$sd, numeric(1)),
    cr2sub_cv_gwl = vapply(stats_list, function(x) x$cv, numeric(1)),
    cr2sub_clean_mean_gwl = vapply(stats_list, function(x) x$clean_mean, numeric(1)),
    cr2sub_clean_sd_gwl = vapply(stats_list, function(x) x$clean_sd, numeric(1)),
    cr2sub_clean_cv_gwl = vapply(stats_list, function(x) x$clean_cv, numeric(1)),
    stringsAsFactors = FALSE
  )
  merge(metadata_ref[c("cr2sub_id", "well_id")], stats_df, by = "well_id", all.x = TRUE)
}

safe_raster <- function(path) {
  tryCatch(
    terra::rast(path),
    error = function(e) {
      warning("No se pudo leer el ráster ", path, ": ", conditionMessage(e))
      NULL
    }
  )
}

extract_or_na <- function(raster_obj, vect_obj) {
  if (is.null(raster_obj)) {
    return(rep(NA_real_, length(vect_obj)))
  }
  vals <- terra::extract(raster_obj, vect_obj)[, 2]
  vals
}

# -----------------------------------------------------------------------------
# lectura de metadatos y series
# -----------------------------------------------------------------------------
metadata_raw <- read.csv(metadata_file, stringsAsFactors = FALSE)

metadata_df <- data.frame(
  well_id = metadata_raw$well_id,
  dga_well_code = metadata_raw$dga_code,
  dga_well_name = metadata_raw$well_name,
  dga_well_basin = metadata_raw$well_basin,
  dga_well_subbasin = metadata_raw$well_subbasin,
  dga_well_lat = as.numeric(metadata_raw$well_lat),
  dga_well_lon = as.numeric(metadata_raw$well_lon),
  dga_well_utm_north = as.numeric(metadata_raw$well_north),
  dga_well_utm_east = as.numeric(metadata_raw$well_east),
  dga_well_elev = as.numeric(metadata_raw$well_elev),
  well_id = metadata_raw$well_id,
  stringsAsFactors = FALSE
)

well_ts <- read.csv(well_ts_file, stringsAsFactors = FALSE)
if (!"date" %in% names(well_ts)) stop("La columna 'date' es requerida en ", well_ts_file)
well_ts$date <- as.Date(well_ts$date)

# -----------------------------------------------------------------------------
# atributos derivados
# -----------------------------------------------------------------------------
coordinates_df <- data.frame(
  cr2sub_id = metadata_df$cr2sub_id,
  cr2sub_lat = metadata_df$dga_well_lat,
  cr2sub_lon = metadata_df$dga_well_lon,
  cr2sub_utm_north_h19 = metadata_df$dga_well_utm_north,
  cr2sub_utm_south_h19 = metadata_df$dga_well_utm_east,
  stringsAsFactors = FALSE
)

gwl_stats_df <- compute_gwl_stats(well_ts, metadata_df)

# --- topografía real --------------------------------------------------------
well_vect <- terra::vect(metadata_df,
  geom = c("dga_well_lon", "dga_well_lat"),
  crs = "EPSG:4326", keepgeom = TRUE
)

dem_raster <- safe_raster(dem_file)
slope_raster <- safe_raster(slope_file)

dem_vals <- extract_or_na(dem_raster, well_vect)
slope_vals <- extract_or_na(slope_raster, well_vect)

topography_df <- data.frame(
  cr2sub_id = metadata_df$cr2sub_id,
  cr2sub_elev = ifelse(is.na(dem_vals), metadata_df$dga_well_elev, dem_vals),
  cr2sub_slp = slope_vals,
  stringsAsFactors = FALSE
)

# --- placeholders para datos pendientes ------------------------------------
camels_df <- data.frame(
  cr2sub_id = metadata_df$cr2sub_id,
  cr2sub_in_basin_camels = NA_character_,
  cr2sub_camels_pr_yr = NA_real_,
  cr2sub_camels_aridity = NA_real_,
  cr2sub_camels_snowf = NA_real_,
  cr2sub_camels_elev = NA_real_,
  cr2sub_camels_slp = NA_real_,
  stringsAsFactors = FALSE
)

bna_df <- data.frame(
  cr2sub_id = metadata_df$cr2sub_id,
  cr2sub_in_basin_bna = NA_character_,
  cr2sub_bna_pr_yr = NA_real_,
  cr2sub_bna_aridity = NA_real_,
  cr2sub_bna_snowf = NA_real_,
  cr2sub_bna_elev = NA_real_,
  cr2sub_bna_slp = NA_real_,
  stringsAsFactors = FALSE
)

soil_df <- data.frame(
  cr2sub_id = metadata_df$cr2sub_id,
  cr2sub_clsoilmap_awc_0_100cm = NA_real_,
  cr2sub_clsoilmap_awc_100_200cm = NA_real_,
  cr2sub_clsoilmap_bulkd_0_100cm = NA_real_,
  cr2sub_clsoilmap_bulkd_100_200cm = NA_real_,
  cr2sub_clsoilmap_clay_0_100cm = NA_real_,
  cr2sub_clsoilmap_clay_100_200cm = NA_real_,
  cr2sub_clsoilmap_ksat_0_100cm = NA_real_,
  cr2sub_clsoilmap_ksat_100_200cm = NA_real_,
  cr2sub_clsoilmap_sand_0_100cm = NA_real_,
  cr2sub_clsoilmap_sand_100_200cm = NA_real_,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# consolidación
# -----------------------------------------------------------------------------
merge_list <- list(
  metadata_df[, c(
    "cr2sub_id", "dga_well_code", "dga_well_name", "dga_well_basin",
    "dga_well_subbasin", "dga_well_lat", "dga_well_lon",
    "dga_well_utm_north", "dga_well_utm_east", "dga_well_elev"
  )],
  coordinates_df,
  topography_df,
  gwl_stats_df[, !(names(gwl_stats_df) %in% c("well_id"))],
  camels_df,
  bna_df,
  soil_df
)

cr2sub_attributes_df <- Reduce(
  function(x, y) merge(x, y, by = "cr2sub_id", all.x = TRUE),
  merge_list
)

required_order <- c(
  "cr2sub_id",
  "dga_well_code", "dga_well_name", "dga_well_basin", "dga_well_subbasin",
  "dga_well_lat", "dga_well_lon", "dga_well_utm_north", "dga_well_utm_east",
  "dga_well_elev",
  "cr2sub_lat", "cr2sub_lon", "cr2sub_utm_north_h19", "cr2sub_utm_south_h19",
  "cr2sub_elev", "cr2sub_slp",
  "cr2sub_mean_gwl", "cr2sub_sd_gwl", "cr2sub_cv_gwl",
  "cr2sub_clean_mean_gwl", "cr2sub_clean_sd_gwl", "cr2sub_clean_cv_gwl",
  "cr2sub_in_basin_camels", "cr2sub_camels_pr_yr", "cr2sub_camels_aridity",
  "cr2sub_camels_snowf", "cr2sub_camels_elev", "cr2sub_camels_slp",
  "cr2sub_in_basin_bna", "cr2sub_bna_pr_yr", "cr2sub_bna_aridity",
  "cr2sub_bna_snowf", "cr2sub_bna_elev", "cr2sub_bna_slp",
  "cr2sub_clsoilmap_awc_0_100cm", "cr2sub_clsoilmap_awc_100_200cm",
  "cr2sub_clsoilmap_bulkd_0_100cm", "cr2sub_clsoilmap_bulkd_100_200cm",
  "cr2sub_clsoilmap_clay_0_100cm", "cr2sub_clsoilmap_clay_100_200cm",
  "cr2sub_clsoilmap_ksat_0_100cm", "cr2sub_clsoilmap_ksat_100_200cm",
  "cr2sub_clsoilmap_sand_0_100cm", "cr2sub_clsoilmap_sand_100_200cm"
)

missing_cols <- setdiff(required_order, names(cr2sub_attributes_df))
if (length(missing_cols)) {
  stop("Columnas faltantes en la tabla final: ", paste(missing_cols, collapse = ", "))
}

cr2sub_attributes_df <- cr2sub_attributes_df[, required_order]

output_dir <- "cr2sub"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "cr2sub_v1_attributes.csv")
write.csv(cr2sub_attributes_df, output_file, row.names = FALSE)

message("Atributos exportados en ", output_file)
message(
  "Dimensiones finales: ", nrow(cr2sub_attributes_df), " pozos x ",
  ncol(cr2sub_attributes_df), " atributos"
)

str(cr2sub_attributes_df)
