# -----------------------------------------------------------------------------
# Script: 03_process_attributes_cr2sub_v1.R
# Purpose: Derive CR2SUB well attributes by merging metadata, terrain rasters,
#   and hydroclimate summaries to produce curated attributes.
# Author: Rodrigo Marinao Rivas, Camila Alvarez Garretón
# Date: [2025-09-16]
# -----------------------------------------------------------------------------

library(terra)
library(zoo)

source("scripts/functions/functions_process_attributes.R")

tag <- "cr2sub"
version <- "v1"

# -----------------------------------------------------------------------------
# archivos de entrada
# -----------------------------------------------------------------------------
metadata_file <- "cr2sub/cr2sub_v1_metadata.csv"
well_ts_file <- "cr2sub/cr2sub_v1_mon_clean.csv"
dem_file <- "input/other_data/dem_masl_fabdemv1.2_2015_300m_epsg4326.tif"
slope_file <- "input/other_data/slope_deg_fabdemv1.2_2015_300m_epsg4326.tif"
pr_file <- "input/other_data/pr_mm_cr2metv2.5_ann_1960_2024_0.05deg_epsg4326.nc"
pet_file <- "input/other_data/et0_mm_cr2met_2.5_ann_1960_2024_0.05deg_epsg4326.nc"
snow_file <- "input/other_data/snow_mm_cr2met_2.5_ann_1960_2024_0.05deg_epsg4326.nc"
clsoil_file <- "input/CLSoilMaps/cr2sub_CLSoilMap_data.csv"

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

pr_mean_rast <- mean_over_years(pr_file, 1980L, 2010L)
pet_mean_rast <- mean_over_years(pet_file, 1980L, 2010L)
snow_mean_rast <- mean_over_years(snow_file, 1980L, 2010L)





# topography


dem_vals <- terra::extract(dem_raster, well_vect)[, 2]
slope_vals <- terra::extract(slope_raster, well_vect)[, 2]

metadata_df$cr2sub_elev <- dem_vals
metadata_df$cr2sub_slp <- slope_vals


# mean, sd and cv

gwl_ts_file <- "cr2sub/cr2sub_v1_mon_clean.csv"
gwl_ts_clean_file <- "cr2sub/cr2sub_v1_mon.csv"

gwl_ts <- read.csv.zoo(gwl_ts_file, check.names = FALSE)
gwl_ts_clean <- read.csv.zoo(gwl_ts_clean_file, check.names = FALSE)


avg_gwl <- apply(coredata(gwl_ts), 2, mean, na.rm = TRUE) [as.character(metadata_df$cr2sub_id)]
sd_gwl <- apply(coredata(gwl_ts), 2, sd, na.rm = TRUE) [as.character(metadata_df$cr2sub_id)]
cv_gwl <- (avg_gwl / sd_gwl)[as.character(metadata_df$cr2sub_id)]

avg_gwl_clean <- apply(coredata(gwl_ts_clean), 2, mean, na.rm = TRUE)
sd_gwl_clean <- apply(coredata(gwl_ts_clean), 2, sd, na.rm = TRUE)
cv_gwl_clean <- avg_gwl_clean / sd_gwl_clean

metadata_df$cr2sub_avg_gwl <- round(avg_gwl, 2)
metadata_df$cr2sub_sd_gwl <- round(sd_gwl, 2)
metadata_df$cr2sub_cv_gwl <- round(cv_gwl, 2)

metadata_df$cr2sub_avg_gwl_clean <- round(avg_gwl_clean, 2)
metadata_df$cr2sub_sd_gwl_clean <- round(sd_gwl_clean, 2)
metadata_df$cr2sub_cv_gwl_clean <- round(cv_gwl_clean, 2)

# camels-cl related

basin_camels <- vect("input/cuencas_camels_cl/cuencas_camels_cl_v1.shp")

camels_ids <- create_well_area_join(
  wells = well_vect,
  polygons = basin_camels,
  grid_resolution = 0.01,
  use_buffer = FALSE,
  match_ids = metadata_df$cr2sub_id
)

camels_ids_char <- as.character(camels_ids)


metadata_df$cr2sub_in_basin_camels <- camels_ids


camels_pr_vals <- terra::extract(pr_mean_rast, basin_camels, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(camels_pr_vals) <- as.data.frame(basin_camels)[, "gauge_id"] |> as.character()

camels_pet_vals <- terra::extract(pet_mean_rast, basin_camels, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(camels_pet_vals) <- as.data.frame(basin_camels)[, "gauge_id"]

camels_snow_vals <- terra::extract(snow_mean_rast, basin_camels, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(camels_snow_vals) <- as.data.frame(basin_camels)[, "gauge_id"]

camels_elev_vals <- terra::extract(dem_raster, basin_camels, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(camels_elev_vals) <- as.data.frame(basin_camels)[, "gauge_id"]

camels_slope_vals <- terra::extract(slope_raster, basin_camels, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(camels_slope_vals) <- as.data.frame(basin_camels)[, "gauge_id"]

metadata_df$cr2sub_camels_pr_yr <- camels_pr_vals[camels_ids_char]
metadata_df$cr2sub_camels_pet_yr <- camels_pet_vals[camels_ids_char]
metadata_df$cr2sub_camels_snowf <- camels_snow_vals[camels_ids_char] / camels_pr_vals[camels_ids_char]
metadata_df$cr2sub_camels_aridity <- camels_pet_vals[camels_ids_char] / camels_pr_vals[camels_ids_char]
metadata_df$cr2sub_camels_elev <- camels_elev_vals[camels_ids_char]
metadata_df$cr2sub_camels_slp <- camels_slope_vals[camels_ids_char]


# bna related

basin_bna <- vect("input/cuencas_BNA/cuencas_BNA.shp") |> project(y = "epsg:4326")


bna_ids <- create_well_area_join(
  wells = well_vect,
  polygons = basin_bna,
  grid_resolution = 0.01,
  use_buffer = TRUE,
  match_ids = metadata_df$cr2sub_id
)

bna_ids_char <- as.character(bna_ids)


metadata_df$cr2sub_in_basin_bna <- bna_ids


bna_pr_vals <- terra::extract(pr_mean_rast, basin_bna, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(bna_pr_vals) <- as.data.frame(basin_bna)[, "COD_CUEN"] |> as.character()

bna_pet_vals <- terra::extract(pet_mean_rast, basin_bna, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(bna_pet_vals) <- as.data.frame(basin_bna)[, "COD_CUEN"]

bna_snow_vals <- terra::extract(snow_mean_rast, basin_bna, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(bna_snow_vals) <- as.data.frame(basin_bna)[, "COD_CUEN"]

bna_elev_vals <- terra::extract(dem_raster, basin_bna, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(bna_elev_vals) <- as.data.frame(basin_bna)[, "COD_CUEN"]

bna_slope_vals <- terra::extract(slope_raster, basin_bna, fun = mean, na.rm = TRUE)[, 2] |> as.numeric()
names(bna_slope_vals) <- as.data.frame(basin_bna)[, "COD_CUEN"]

metadata_df$cr2sub_bna_pr_yr <- bna_pr_vals[bna_ids_char]
metadata_df$cr2sub_bna_pet_yr <- bna_pet_vals[bna_ids_char]
metadata_df$cr2sub_bna_snowf <- bna_snow_vals[bna_ids_char] / bna_pr_vals[bna_ids_char]
metadata_df$cr2sub_bna_aridity <- bna_pet_vals[bna_ids_char] / bna_pr_vals[bna_ids_char]
metadata_df$cr2sub_bna_elev <- bna_elev_vals[bna_ids_char]
metadata_df$cr2sub_bna_slp <- bna_slope_vals[bna_ids_char]



# mean, sd y cv


clsmap <- read.csv("input/CLSoilMaps/cr2sub_CLSoilMap_data.csv", row.names = 1)

clsmap <- clsmap[as.character(metadata_df$cr2sub_id), ]

soil_cols <- c(
  "cr2sub_clsoilmap_awc_0_100cm" = "AWC_0_100",
  "cr2sub_clsoilmap_awc_100_200cm" = "AWC_100_200",
  "cr2sub_clsoilmap_bulkd_0_100cm" = "Bulkd_0_100",
  "cr2sub_clsoilmap_bulkd_100_200" = "Bulkd_100_200",
  "cr2sub_clsoilmap_clay_0_100cm" = "Clay_0_100",
  "cr2sub_clsoilmap_clay_100_200cm" = "Clay_100_200",
  "cr2sub_clsoilmap_ksat_0_100cm" = "ksat_0_100",
  "cr2sub_clsoilmap_ksat_100_200cm" = "ksat_100_200",
  "cr2sub_clsoilmap_sand_0_100cm" = "Sand_0_100",
  "cr2sub_clsoilmap_sand_100_200cm" = "Sand_100_200"
)

clsmap <- clsmap[, soil_cols]
colnames(clsmap) <- names(soil_cols)

metadata_df <- cbind.data.frame(metadata_df, clsmap)


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
