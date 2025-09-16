# -----------------------------------------------------------------------------
# Script: 01_read_xls.R
# Purpose: Read groundwater level XLS files from DGA, extract metadata and time
#   series, aggregate monthly values, fill missing months, and export outputs.
# Author: Rodrigo Marinao Rivas, Camila Alvarez Garretón
# Date: [2025-07-08]
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section: User input – select which data folder to process
#
# IMPORTANT: The user should modify `data_folder` to point to the folder
# containing the new update (e.g. "update_2025-06-10").
#
# The option `base_2025-03-31` corresponds to the full historical base and
# can be used to audit the existing series, but takes "some" time.
# Subsequent runs should focus only on updates to speed up processing.
# -----------------------------------------------------------------------------
# data_folder <- "base_2025-03-31"
data_folder <- "update_2025-06-10"


# -----------------------------------------------------------------------------
# Section: Setup and configuration
# -----------------------------------------------------------------------------
# Automatically detect script location to set working directory
# If running interactively (e.g., in RStudio), setwd() moves one folder up from
this_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                rstudioapi::isAvailable()) {
  rstudioapi::getSourceEditorContext()$path
} else {
  sys.frames()[[1]]$ofile
}

wd <- normalizePath(file.path(dirname(this_file), ".."))
setwd(wd)

library(readxl)
library(zoo)
library(terra)

# Load custom helper functions for directory creation, coordinate conversion,.
source("scripts/00_functions.R")


tmp_dir <- gsub("//", "/", tempdir())

# Specify subfolder name for this data update run

# Construct paths for input XLS, temporary and output CSV directories
input_dir <- file.path("IN", data_folder)
tmp_folder <- file.path(tmp_dir, "tmp", data_folder)
out_folder <- file.path("OUT", data_folder)

# -----------------------------------------------------------------------------
# Section: Read and process XLS files
# Loop through all Excel files and extract water level time series and metadata
# -----------------------------------------------------------------------------

# Ensure output directory exists
create_dir_if_not_exists(tmp_folder)
# Ensure output directory exists
create_dir_if_not_exists(out_folder)

# List all XLS files in input directory
files <- list.files(input_dir, full.names = TRUE)
# Warn if no files found and exit loop
if (length(files) == 0) {
  warning(sprintf("No se encontraron archivos en %s", input_dir))
  next
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
    tryCatch({
      # Read sheet into data frame
      xls <- as.data.frame(read_excel(file, sheet = sheet))
      # Extract metadata, dates, and data series from sheet
      processed_data <- process_excel_sheet(xls, sheet)
      # Guardar metadatos
      sheets_details[[sheet_idx]] <- processed_data$metadata
      # Create output data.frame with dates and water level values
      well_id <- substr(xls[6, 3], 1, 10)

      df_out <- data.frame(
        fecha = processed_data$dates,
        nivel = processed_data$data
      )
      colnames(df_out)[2] <- well_id
      # Define directory for this well code under temporary folder
      pozo_dir <- file.path(tmp_folder, well_id)
      # Ensure output directory exists
      create_dir_if_not_exists(pozo_dir)
      # Write CSV for this well, named by code and date range
      file_out <- file.path(pozo_dir,
                            sprintf("Nivel_%s_%s_%s.csv",
                                    well_id,
                                    min(processed_data$dates),
                                    max(processed_data$dates)))

      write.csv(df_out, file_out, row.names = FALSE)

    }, error = function(e) {
      # Warn and continue if sheet processing fails
      warning(sprintf("Error procesando archivo %s, hoja %s: %s",
                      file, sheet, e$message))
    })
  }

  files_details[[file_idx]] <- sheets_details
}

# -----------------------------------------------------------------------------
# Section: Combine and clean metadata
# Extract metadata from all processed sheets and remove duplicates
# -----------------------------------------------------------------------------

# Define columns to extract for station metadata
metadata_columns <- c(
  "well_id", 
  "well_name", 
  "well_basin", 
  "well_subbasin", 
  "well_elev",
  "well_lat", 
  "well_lon", 
  "well_north", 
  "well_east"
)

# Combine metadata from all files and sheets into one data.frame
files_details_list <- lapply(files_details, function(details) {
  # Assemble raw metadata matrix into data.frame with defined columns
  df_details <- as.data.frame(matrix(NA,
                                     ncol = length(metadata_columns),
                                     nrow = length(details)))
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
index_keep <- !duplicated(apply(df_details_prov[, c("well_id",
                                                    "well_name", 
                                                    "well_north", 
                                                    "well_east")],
                                FUN = paste, MARGIN = 1, collapse = ""))
full_details <- df_details_prov[index_keep, ]

rownames(full_details) <- full_details[, "well_id"]

details <- full_details

# Convert altitude to numeric and latitude/longitude to decimal degrees
details$well_elev <- as.numeric(full_details$well_elev)
details$well_lat <- unname(sapply(full_details$well_lat, deg2dec))
details$well_lon <- unname(sapply(full_details$well_lon, deg2dec))
details$well_north <- as.numeric(full_details$well_north)
details$well_east <- as.numeric(full_details$well_east)

# -----------------------------------------------------------------------------
# Section: Spatial validation – keep only stations within Chilean boundaries
# -----------------------------------------------------------------------------

# Convert metadata table to spatial vector using longitude and latitude columns
# CRS is WGS84 (EPSG:4326), suitable for global geographic coordinates
details_vect <- vect(details, geom = c("well_lon", "well_lat"),
                     crs = "epsg:4326")

# Load simplified polygon of Chilean boundaries from geopackage
chile_boundaries <- vect("IN/Chile_boundaries_simplified_from_SIIT.gpkg")

# Spatially mask stations to keep only those within Chilean territory
# Any station falling outside the polygon will be excluded
points_in_chile <- mask(details_vect, chile_boundaries)

# Create a logical flag to identify valid stations inside Chile
flag_details <- details$well_id %in% points_in_chile$well_id

# Filter original metadata to keep only valid Chilean stations
details_filtered <- details[flag_details, ]


# Export cleaned metadata table to CSV
write.csv(details_filtered, file.path(out_folder, "cr2sub_metadata.csv"),
          row.names = FALSE)

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
                                          header = FALSE, nrows = 1)[1, 2]
      ts_full <- ts
    } else {
      ts_full00 <- cbind.zoo(ts_full, ts)
      ts_full <- zoo(apply(coredata(ts_full00), 1, mean, na.rm = TRUE),
                     time(ts_full00))
    }
  }

  if (pozo_idx == 1) {
    ts_all <- zoo(matrix(coredata(ts_full), ncol = 1), time(ts_full))
  } else {
    ts_all <- cbind.zoo(ts_all, ts_full)
  }
}

# Replace non-positive values with NA and invert levels
colnames(ts_all) <- codigos_pozos
coredata(ts_all)[ts_all <= 0] <- NA

full_series <- -ts_all[, details_filtered$well_id]

# Save raw combined time series to CSV
df <- data.frame(date = time(full_series), coredata(full_series), check.names = FALSE)
write.csv(df, file.path(out_folder, "cr2sub_raw.csv"), row.names = FALSE)



# -----------------------------------------------------------------------------
# Section: Monthly Aggregation
# -----------------------------------------------------------------------------
# Aggregate raw series to monthly means and fill missing months
index_mon <- format(time(full_series), "%Y-%m")


mon_series_cd <- aggregate(full_series, index_mon, FUN = mean,
                           na.rm = TRUE)

mon_series0 <- zoo(mon_series_cd, as.Date(paste0(time(mon_series_cd), "-01")))

# Ensure continuous monthly index, filling missing months with NAs
start_month <- as.Date(paste0(format(start(mon_series0), "%Y-%m"), "-01"))
end_month   <- as.Date(paste0(format(end(mon_series0),   "%Y-%m"), "-01"))
all_months  <- seq(from = start_month, to = end_month, by = "month")
mon_series  <- merge(zoo(, all_months), mon_series0)

# Save completed monthly series to CSV
df_mon <- data.frame(date = time(mon_series), coredata(mon_series), check.names = FALSE)
write.csv(df_mon, file.path(out_folder, "cr2sub_mon.csv"), row.names = FALSE)

