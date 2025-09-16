# -----------------------------------------------------------------------------
# Script: 03_process_attributes_cr2sub_v1.R
# Purpose: Derive CR2SUB well attributes by merging metadata, terrain rasters,
#   and hydroclimate summaries to produce curated attributes.
# Author: Rodrigo Marinao Rivas, Camila Alvarez Garretón
# Date: [2025-09-16]
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# lectura de metadatos y series
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
  dga_well_utm_north = as.numeric(metadata_raw$dga_well_north),
  dga_well_utm_east = as.numeric(metadata_raw$dga_well_east),
  dga_well_elev = as.numeric(metadata_raw$dga_well_elev),
  stringsAsFactors = FALSE
)

well_ts <- read.csv(well_ts_file, stringsAsFactors = FALSE)
well_ts$date <- as.Date(well_ts$date)

# -----------------------------------------------------------------------------
# atributos derivados
# -----------------------------------------------------------------------------
metadata_df$cr2sub_lat <- metadata_df$dga_well_lat
metadata_df$cr2sub_lon <- metadata_df$dga_well_lon
metadata_df$cr2sub_utm_north_h19 <- metadata_df$dga_well_utm_north
metadata_df$cr2sub_utm_south_h19 <- metadata_df$dga_well_utm_east

gwl_ids <- intersect(metadata_df$cr2sub_id, names(well_ts))
if (length(gwl_ids)) {
  stats_list <- lapply(gwl_ids, function(id) {
    series_vals <- as.numeric(well_ts[[id]])
    series_vals <- series_vals[is.finite(series_vals)]
    if (!length(series_vals)) {
      return(rep(NA_real_, 6))
    }
    mean_val <- mean(series_vals)
    sd_val <- sd(series_vals)
    clean_vals <- if (is.finite(sd_val) && sd_val > 0) {
      kept <- series_vals[abs(series_vals - mean_val) <= (2 * sd_val)]
      if (length(kept)) kept else series_vals
    } else {
      series_vals
    }
    clean_mean <- mean(clean_vals)
    clean_sd <- sd(clean_vals)
    c(
      mean_val,
      sd_val,
      if (is.finite(mean_val) && abs(mean_val) > .Machine$double.eps) sd_val / mean_val else NA_real_,
      clean_mean,
      clean_sd,
      if (is.finite(clean_mean) && abs(clean_mean) > .Machine$double.eps) clean_sd / clean_mean else NA_real_
    )
  })
  stats_mat <- do.call(rbind, stats_list)
  colnames(stats_mat) <- c(
    "cr2sub_mean_gwl", "cr2sub_sd_gwl", "cr2sub_cv_gwl",
    "cr2sub_clean_mean_gwl", "cr2sub_clean_sd_gwl", "cr2sub_clean_cv_gwl"
  )
  stats_df <- data.frame(cr2sub_id = gwl_ids, stats_mat, stringsAsFactors = FALSE)
  idx <- match(metadata_df$cr2sub_id, stats_df$cr2sub_id)
  metadata_df$cr2sub_mean_gwl <- stats_df$cr2sub_mean_gwl[idx]
  metadata_df$cr2sub_sd_gwl <- stats_df$cr2sub_sd_gwl[idx]
  metadata_df$cr2sub_cv_gwl <- stats_df$cr2sub_cv_gwl[idx]
  metadata_df$cr2sub_clean_mean_gwl <- stats_df$cr2sub_clean_mean_gwl[idx]
  metadata_df$cr2sub_clean_sd_gwl <- stats_df$cr2sub_clean_sd_gwl[idx]
  metadata_df$cr2sub_clean_cv_gwl <- stats_df$cr2sub_clean_cv_gwl[idx]
} else {
  metadata_df$cr2sub_mean_gwl <- NA_real_
  metadata_df$cr2sub_sd_gwl <- NA_real_
  metadata_df$cr2sub_cv_gwl <- NA_real_
  metadata_df$cr2sub_clean_mean_gwl <- NA_real_
  metadata_df$cr2sub_clean_sd_gwl <- NA_real_
  metadata_df$cr2sub_clean_cv_gwl <- NA_real_
}

well_vect <- terra::vect(
  metadata_df,
  geom = c("dga_well_lon", "dga_well_lat"),
  crs = "EPSG:4326",
  keepgeom = TRUE
)

dem_raster <- terra::rast(dem_file)
slope_raster <- terra::rast(slope_file)

dem_vals <- terra::extract(dem_raster, well_vect)[, 2]
slope_vals <- terra::extract(slope_raster, well_vect)[, 2]

metadata_df$cr2sub_elev <- ifelse(is.na(dem_vals), metadata_df$dga_well_elev, dem_vals)
metadata_df$cr2sub_slp <- slope_vals

metadata_df$cr2sub_in_basin_camels <- NA_character_
metadata_df$cr2sub_camels_pr_yr <- NA_real_
metadata_df$cr2sub_camels_aridity <- NA_real_
metadata_df$cr2sub_camels_snowf <- NA_real_
metadata_df$cr2sub_camels_elev <- NA_real_
metadata_df$cr2sub_camels_slp <- NA_real_

metadata_df$cr2sub_in_basin_bna <- NA_character_
metadata_df$cr2sub_bna_pr_yr <- NA_real_
metadata_df$cr2sub_bna_aridity <- NA_real_
metadata_df$cr2sub_bna_snowf <- NA_real_
metadata_df$cr2sub_bna_elev <- NA_real_
metadata_df$cr2sub_bna_slp <- NA_real_

metadata_df$cr2sub_clsoilmap_awc_0_100cm <- NA_real_
metadata_df$cr2sub_clsoilmap_awc_100_200cm <- NA_real_
metadata_df$cr2sub_clsoilmap_bulkd_0_100cm <- NA_real_
metadata_df$cr2sub_clsoilmap_bulkd_100_200cm <- NA_real_
metadata_df$cr2sub_clsoilmap_clay_0_100cm <- NA_real_
metadata_df$cr2sub_clsoilmap_clay_100_200cm <- NA_real_
metadata_df$cr2sub_clsoilmap_ksat_0_100cm <- NA_real_
metadata_df$cr2sub_clsoilmap_ksat_100_200cm <- NA_real_
metadata_df$cr2sub_clsoilmap_sand_0_100cm <- NA_real_
metadata_df$cr2sub_clsoilmap_sand_100_200cm <- NA_real_

# -----------------------------------------------------------------------------
# consolidación
# -----------------------------------------------------------------------------
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

cr2sub_attributes_df <- metadata_df[, required_order]

output_dir <- "cr2sub"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "cr2sub_v1_attributes.csv")
write.csv(cr2sub_attributes_df, output_file, row.names = FALSE)

message("Atributos exportados en ", output_file)
message(
  "Dimensiones finales: ", nrow(cr2sub_attributes_df), " pozos x ",
  ncol(cr2sub_attributes_df), " atributos"
)
