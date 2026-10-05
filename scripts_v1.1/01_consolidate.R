# -----------------------------------------------------------------------------
# Script: 01_consolidate.R
# Purpose: Consolidate groundwater level XLS files from DGA, extracting metadata
#   and time series to build monthly aggregates.
# Author: Rodrigo Marinao Rivas, Camila Alvarez Garretón
# Date: [2025-07-08]
# Updated : [2025-09-16]
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section: Dependencies
# -----------------------------------------------------------------------------

library(readxl)
library(zoo)
library(terra)

# Load custom helper functions for directory creation, coordinate conversion
source("scripts_v1.1/functions/functions_consolidate_cr2sub.R")

# -----------------------------------------------------------------------------
# Section: Configuration
# -----------------------------------------------------------------------------

tag <- "cr2sub"
version <- "v1.1"
out_folder <- file.path(tag)
project_tmp_dir <- file.path(".tmp")
system_tmp_dir <- file.path(tempdir(), "cr2sub")

create_dir_if_not_exists(out_folder)
create_dir_if_not_exists(project_tmp_dir)
create_dir_if_not_exists(system_tmp_dir)

# -----------------------------------------------------------------------------
# Section: Discover Input File
# -----------------------------------------------------------------------------

input_file <- file.path("input/DGA_GWL_obs_descarga", "AM006T0034616_1.xlsx")

if (!file.exists(input_file)) {
  stop(sprintf("No se encontró el archivo %s", input_file))
}

# -----------------------------------------------------------------------------
# Section: Extract Time Series And Metadata
# -----------------------------------------------------------------------------

raw <- suppressMessages(as.data.frame(read_excel(input_file, sheet = 1)))

colnames(raw) <- c(
  "dga_well_name", "dga_well_code", "dga_well_province", "dga_well_commune",
  "dga_well_basin", "dga_well_subbasin", "dga_well_subsubbasin",
  "dga_well_utm_north", "dga_well_utm_east", "dga_well_lat_raw",
  "dga_well_lon_raw", "dga_well_elev", "fecha_raw", "nivel"
)

raw$fecha <- as.Date(trimws(raw$fecha_raw), format = "%d-%m-%Y")
raw$nivel <- as.numeric(raw$nivel)

valid_rows <- !is.na(raw$dga_well_code) & !is.na(raw$fecha) & !is.na(raw$nivel)
raw <- raw[valid_rows, ]

# Average duplicate (well, date) measurements before building each series
raw_series <- aggregate(
  nivel ~ dga_well_code + fecha,
  data = raw,
  FUN = mean
)

well_codes <- unique(raw_series$dga_well_code)

for (code_idx in seq_along(well_codes)) {
  dga_well_code <- well_codes[code_idx]
  well_rows <- raw_series[raw_series$dga_well_code == dga_well_code, ]

  pozo_dir <- file.path(system_tmp_dir, dga_well_code)
  create_dir_if_not_exists(pozo_dir)

  df_out <- data.frame(
    fecha = well_rows$fecha,
    nivel = well_rows$nivel
  )
  colnames(df_out)[2] <- dga_well_code

  file_out <- file.path(
    pozo_dir,
    sprintf(
      "Nivel_%s_%s_%s.csv",
      dga_well_code,
      min(df_out$fecha),
      max(df_out$fecha)
    )
  )

  write.csv(df_out, file_out, row.names = FALSE)

  progress <- 100 * code_idx / length(well_codes)
  cat(sprintf("Avance: %.1f%% (%d/%d)\n", progress, code_idx, length(well_codes)))
  flush.console()
}

# -----------------------------------------------------------------------------
# Section: Consolidate Metadata
# -----------------------------------------------------------------------------

first_rows <- raw[!duplicated(raw$dga_well_code), ]

full_details <- data.frame(
  dga_well_code = first_rows$dga_well_code,
  dga_well_name = first_rows$dga_well_name,
  dga_well_basin = first_rows$dga_well_basin,
  dga_well_subbasin = first_rows$dga_well_subbasin,
  dga_well_elev = as.numeric(first_rows$dga_well_elev),
  dga_well_lat = unname(sapply(first_rows$dga_well_lat_raw, deg2dec_compact)),
  dga_well_lon = unname(sapply(first_rows$dga_well_lon_raw, deg2dec_compact)),
  dga_well_utm_north = as.numeric(first_rows$dga_well_utm_north),
  dga_well_utm_east = as.numeric(first_rows$dga_well_utm_east),
  stringsAsFactors = FALSE
)
rownames(full_details) <- full_details$dga_well_code

details <- full_details
details$cr2sub_id <- as.numeric(sub("-[0-9Kk]$", "", details$dga_well_code))

# -----------------------------------------------------------------------------
# Section: Spatial Validation
# -----------------------------------------------------------------------------

details_vect <- vect(details,
  geom = c("dga_well_lon", "dga_well_lat"),
  crs = "epsg:4326"
)

chile_boundaries <-
  vect("input/other_data/Chile_boundaries_simplified_from_SIIT.gpkg")

points_in_chile <- mask(details_vect, chile_boundaries)

flag_details <- details$cr2sub_id %in% points_in_chile$cr2sub_id
details_filtered <- details[flag_details, ]

write.csv(details_filtered, file.path(
  project_tmp_dir,
  paste0(tag, "_", version, "_metadata.csv")
),
row.names = FALSE
)

# -----------------------------------------------------------------------------
# Section: Merge Time Series By Well
# -----------------------------------------------------------------------------

wells_dirs <- list.dirs(system_tmp_dir, recursive = FALSE, full.names = TRUE)
wells_dirs <- wells_dirs[wells_dirs != system_tmp_dir]
wells_out <- basename(wells_dirs)
nombres_pozos <- wells_out
codigos_pozos <- rep(NA, length(wells_out))

for (pozo_idx in seq_along(wells_out)) {
  folder_out <- wells_out[pozo_idx]
  files_out <- list.files(file.path(system_tmp_dir, folder_out),
    full.names = TRUE
  )

  for (file_idx in seq_along(files_out)) {
    file_out <- files_out[file_idx]
    ts <- read.csv.zoo(file_out, format = "%Y-%m-%d")

    if (file_idx == 1) {
      codigos_pozos[pozo_idx] <-
        read.csv(file_out, header = FALSE, nrows = 1)[1, 2]
      ts_full <- ts
    } else {
      ts_full_tmp <- cbind.zoo(ts_full, ts)
      ts_full <- zoo(
        apply(coredata(ts_full_tmp), 1, mean, na.rm = TRUE),
        time(ts_full_tmp)
      )
    }
  }

  if (pozo_idx == 1) {
    ts_all <- zoo(matrix(coredata(ts_full), ncol = 1), time(ts_full))
  } else {
    ts_all <- cbind.zoo(ts_all, ts_full)
  }
}

codigos_pozos <- as.numeric(sub("-[0-9Kk]$", "", codigos_pozos))
colnames(ts_all) <- codigos_pozos

coredata(ts_all)[ts_all <= 0] <- NA
full_series <- -ts_all[, as.character(details_filtered$cr2sub_id)]

# -----------------------------------------------------------------------------
# Section: Monthly Aggregation And Export
# -----------------------------------------------------------------------------

index_mon <- format(time(full_series), "%Y-%m")
mon_series_cd <- aggregate(full_series, index_mon, FUN = mean, na.rm = TRUE)
mon_series0 <- zoo(mon_series_cd, as.Date(paste0(time(mon_series_cd), "-01")))

start_month <- as.Date(paste0(format(start(mon_series0), "%Y-%m"), "-01"))
end_month <- as.Date(paste0(format(end(mon_series0), "%Y-%m"), "-01"))
all_months <- seq(from = start_month, to = end_month, by = "month")
mon_series <- merge(zoo(, all_months), mon_series0)

df_mon <- data.frame(
  date = time(mon_series),
  coredata(mon_series), check.names = FALSE
)

write.csv(
  df_mon,
  file.path(out_folder, paste0(tag, "_", version, "_gwl_mon.csv")),
  quote = FALSE,
  row.names = FALSE
)
