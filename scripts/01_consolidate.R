# -----------------------------------------------------------------------------
# Script: 01_consolidate_cr2sub_v1.R
# Purpose: Consolidate groundwater level XLS files from DGA, extracting metadata
#   and time series to build monthly aggregates, backfill missing months, and
#   export curated outputs for CR2SUB.
# Author: Rodrigo Marinao Rivas, Camila Alvarez Garretón
# Date: [2025-07-08]
# Updated : [2025-09-16]
# -----------------------------------------------------------------------------

library(readxl)
library(zoo)
library(terra)

# Load custom helper functions for directory creation, coordinate conversion,.
source("scripts/functions/functions_consolidate_cr2sub.R")


tmp_dir <- gsub("//", "/", tempdir())

# Specify subfolder name for this data update run

# Construct paths for temporary and output CSV directories
tmp_folder <- file.path(tmp_dir, "tmp")


tag_name <- "cr2sub"
version <- "v1"

out_folder <- file.path(tag_name)

# -----------------------------------------------------------------------------
# Section: Read and process XLS files
# Loop through all Excel files and extract water level time series and metadata
# -----------------------------------------------------------------------------

# Ensure output directory exists
create_dir_if_not_exists(tmp_folder)
# Ensure output directory exists
create_dir_if_not_exists(out_folder)

# List all XLS/XLSX files across input directories matching input/dga_xls_*
input_dirs <- Sys.glob(file.path("input/DGA_GWL_observations", "dga_xls_*"))
files <- unlist(lapply(input_dirs, function(d) {
  list.files(d, full.names = TRUE, pattern = "\\.(xls|xlsx)$", ignore.case = TRUE)
}))

# Stop if no files were found
if (length(files) == 0) {
  stop(sprintf("No se encontraron archivos en %s", "input/dga_xls_*"))
}

# Loop over each input XLS file
files_details <- vector(mode = "list", length = length(files))


for (file_idx in seq_along(files)) {
  file <- files[file_idx]
  # Get sheet names from current XLS file
  sheets <- excel_sheets(path = file)
  # Loop over each sheet within the XLS file
  sheets_details <- vector(mode = "list", length = length(sheets))

  for (sheet_idx in seq_along(sheets)) {
    sheet <- sheets[sheet_idx]
    # Attempt to read and process each sheet, catch errors per sheet
    tryCatch(
      {
        # Read sheet into data frame (suppress tibble name repair message)
        xls <- suppressMessages(as.data.frame(read_excel(file, sheet = sheet)))
        # Extract metadata, dates, and data series from sheet
        processed_data <- process_excel_sheet(xls, sheet)
        # Guardar metadatos
        sheets_details[[sheet_idx]] <- processed_data$metadata
        # Create output data.frame with dates and water level values
        dga_well_code <- substr(xls[6, 3], 1, 10)

        df_out <- data.frame(
          fecha = processed_data$dates,
          nivel = processed_data$data
        )
        colnames(df_out)[2] <- dga_well_code
        # Define directory for this well code under temporary folder
        pozo_dir <- file.path(tmp_folder, dga_well_code)
        # Ensure output directory exists
        create_dir_if_not_exists(pozo_dir)
        # Write CSV for this well, named by code and date range
        file_out <- file.path(
          pozo_dir,
          sprintf(
            "Nivel_%s_%s_%s.csv",
            dga_well_code,
            min(processed_data$dates),
            max(processed_data$dates)
          )
        )

        write.csv(df_out, file_out, row.names = FALSE)
      },
      error = function(e) {
        # Warn and continue if sheet processing fails
        warning(sprintf(
          "Error procesando archivo %s, hoja %s: %s",
          file, sheet, e$message
        ))
      }
    )
  }

  files_details[[file_idx]] <- sheets_details

  # Progress update after processing each file
  progress <- 100 * file_idx / length(files)
  cat(sprintf("Avance: %.1f%% (%d/%d)\n", progress, file_idx, length(files)))
  flush.console()
}


# -----------------------------------------------------------------------------
# Section: Combine and clean metadata
# Extract metadata from all processed sheets and remove duplicates
# -----------------------------------------------------------------------------

# Define columns to extract for station metadata
metadata_columns <- c(
  "dga_well_code",
  "dga_well_name",
  "dga_well_basin",
  "dga_well_subbasin",
  "dga_well_elev",
  "dga_well_lat",
  "dga_well_lon",
  "dga_well_utm_north",
  "dga_well_utm_east"
)

# Combine metadata from all files and sheets into one data.frame
files_details_list <- lapply(files_details, function(details) {
  # Assemble raw metadata matrix into data.frame with defined columns
  df_details <- as.data.frame(matrix(NA,
    ncol = length(metadata_columns),
    nrow = length(details)
  ))
  colnames(df_details) <- metadata_columns

  for (j in seq_along(details)) {
    df_details[j, ] <- details[[j]]
  }
  df_details
})

# Assemble raw metadata matrix into data.frame with defined columns
df_details_raw <- do.call("rbind.data.frame", files_details_list)
df_details_prov <- df_details_raw[!is.na(df_details_raw[, 1]), ]

# Remove duplicate station entries based on key metadata fields
index_keep <- !duplicated(apply(
  df_details_prov[, c(
    "dga_well_code",
    "dga_well_name"
  )],
  FUN = paste, MARGIN = 1, collapse = ""
))
full_details <- df_details_prov[index_keep, ]

rownames(full_details) <- full_details[, "dga_well_code"]

details <- full_details

# Convert altitude to numeric and latitude/longitude to decimal degrees
details$cr2sub_id <- as.numeric(sub("-[0-9Kk]$", "", details$dga_well_code))
details$dga_well_elev <- as.numeric(full_details$dga_well_elev)
details$dga_well_lat <- unname(sapply(full_details$dga_well_lat, deg2dec))
details$dga_well_lon <- unname(sapply(full_details$dga_well_lon, deg2dec))
details$dga_well_utm_north <- as.numeric(full_details$dga_well_utm_north)
details$dga_well_utm_east <- as.numeric(full_details$dga_well_utm_east)

# -----------------------------------------------------------------------------
# Section: Spatial validation – keep only stations within Chilean boundaries
# -----------------------------------------------------------------------------

# Convert metadata table to spatial vector using longitude and latitude columns
# CRS is WGS84 (EPSG:4326), suitable for global geographic coordinates
details_vect <- vect(details,
  geom = c("dga_well_lon", "dga_well_lat"),
  crs = "epsg:4326"
)

# Load simplified polygon of Chilean boundaries from geopackage
chile_boundaries <- vect("input/other_data/Chile_boundaries_simplified_from_SIIT.gpkg")

# Spatially mask stations to keep only those within Chilean territory
# Any station falling outside the polygon will be excluded
points_in_chile <- mask(details_vect, chile_boundaries)

# Create a logical flag to identify valid stations inside Chile
flag_details <- details$cr2sub_id %in% points_in_chile$cr2sub_id

# Filter original metadata to keep only valid Chilean stations
details_filtered <- details[flag_details, ]


# Export cleaned metadata table to CSV
write.csv(details_filtered, file.path(
  out_folder,
  paste0(tag_name, "_", version, "_metadata.csv")
),
row.names = FALSE
)

# -----------------------------------------------------------------------------
# Section: Merge per-well CSVs into a complete raw time series matrix
# -----------------------------------------------------------------------------

# Read each well's CSV to build full time series per well
wells_out <- list.dirs(tmp_folder, recursive = FALSE, full.names = FALSE)
nombres_pozos <- wells_out
codigos_pozos <- rep(NA, length(wells_out))

# Procesar cada pozo
for (pozo_idx in seq_along(wells_out)) {
  folder_out <- wells_out[pozo_idx]
  files_out <- list.files(file.path(tmp_folder, folder_out), full.names = TRUE)

  # For each additional segment, merge by timestamp, averaging across duplicates
  for (file_idx in seq_along(files_out)) {
    file_out <- files_out[file_idx]
    ts <- read.csv.zoo(file_out, format = "%Y-%m-%d")

    if (file_idx == 1) {
      codigos_pozos[pozo_idx] <- read.csv(file_out,
        header = FALSE, nrows = 1
      )[1, 2]
      ts_full <- ts
    } else {
      ts_full00 <- cbind.zoo(ts_full, ts)
      ts_full <- zoo(
        apply(coredata(ts_full00), 1, mean, na.rm = TRUE),
        time(ts_full00)
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

# Replace non-positive values with NA and invert levels
colnames(ts_all) <- codigos_pozos
coredata(ts_all)[ts_all <= 0] <- NA

full_series <- -ts_all[, as.character(details_filtered$cr2sub_id)]

# -----------------------------------------------------------------------------
# Section: Monthly Aggregation
# -----------------------------------------------------------------------------
# Aggregate raw series to monthly means and fill missing months
index_mon <- format(time(full_series), "%Y-%m")


mon_series_cd <- aggregate(full_series, index_mon,
  FUN = mean,
  na.rm = TRUE
)

mon_series0 <- zoo(mon_series_cd, as.Date(paste0(time(mon_series_cd), "-01")))

# Ensure continuous monthly index, filling missing months with NAs
start_month <- as.Date(paste0(format(start(mon_series0), "%Y-%m"), "-01"))
end_month <- as.Date(paste0(format(end(mon_series0), "%Y-%m"), "-01"))
all_months <- seq(from = start_month, to = end_month, by = "month")
mon_series <- merge(zoo(, all_months), mon_series0)

# Save completed monthly series to CSV
df_mon <- data.frame(date = time(mon_series), coredata(mon_series), check.names = FALSE)
write.csv(df_mon, file.path(out_folder, paste0(tag_name, "_", version, "_mon.csv")), quote = FALSE, row.names = FALSE)
