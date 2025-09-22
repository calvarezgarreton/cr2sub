# -----------------------------------------------------------------------------
# Script: 01_consolidate.R
# Purpose: Consolidate groundwater level XLS files from DGA, extracting metadata
#   and time series to build monthly aggregates, backfill missing months, and
#   export curated outputs for CR2SUB.
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
source("scripts/functions/functions_consolidate_cr2sub.R")

# -----------------------------------------------------------------------------
# Section: Configuration
# -----------------------------------------------------------------------------

tag <- "cr2sub"
version <- "v1.1"
tmp_dir <- gsub("//", "/", tempdir())
tmp_folder <- file.path(tmp_dir, "tmp")
out_folder <- file.path(tag)

create_dir_if_not_exists(tmp_folder)
create_dir_if_not_exists(out_folder)

# -----------------------------------------------------------------------------
# Section: Discover Input Files
# -----------------------------------------------------------------------------

input_dirs <- Sys.glob(file.path("input/DGA_GWL_observations", "dga_xls_*"))

files <- unlist(lapply(input_dirs, function(d) {
  list.files(d,
    full.names = TRUE, pattern = "\\.(xls|xlsx)$",
    ignore.case = TRUE
  )
}))

if (length(files) == 0) {
  stop(sprintf("No se encontraron archivos en %s", "input/dga_xls_*"))
}

# -----------------------------------------------------------------------------
# Section: Extract Time Series And Metadata By File
# -----------------------------------------------------------------------------

files_details <- vector(mode = "list", length = length(files))

for (file_idx in seq_along(files)) {
  file <- files[file_idx]
  sheets <- excel_sheets(path = file)
  sheets_details <- vector(mode = "list", length = length(sheets))

  for (sheet_idx in seq_along(sheets)) {
    sheet <- sheets[sheet_idx]

    tryCatch(
      {
        xls <- suppressMessages(as.data.frame(read_excel(file, sheet = sheet)))
        processed_data <- process_excel_sheet(xls, sheet)

        sheets_details[[sheet_idx]] <- processed_data$metadata

        dga_well_code <- substr(xls[6, 3], 1, 10)
        pozo_dir <- file.path(tmp_folder, dga_well_code)

        create_dir_if_not_exists(pozo_dir)

        df_out <- data.frame(
          fecha = processed_data$dates,
          nivel = processed_data$data
        )
        colnames(df_out)[2] <- dga_well_code

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
        warning(sprintf(
          "Error procesando archivo %s, hoja %s: %s",
          file, sheet, e$message
        ))
      }
    )
  }

  files_details[[file_idx]] <- sheets_details

  progress <- 100 * file_idx / length(files)
  cat(sprintf("Avance: %.1f%% (%d/%d)\n", progress, file_idx, length(files)))
  flush.console()
}

# -----------------------------------------------------------------------------
# Section: Consolidate Metadata
# -----------------------------------------------------------------------------

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

files_details_list <- lapply(files_details, function(details) {
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

df_details_raw <- do.call("rbind.data.frame", files_details_list)
df_details_prov <- df_details_raw[!is.na(df_details_raw[, 1]), ]

index_keep <- !duplicated(apply(
  df_details_prov[, c("dga_well_code", "dga_well_name")],
  FUN = paste,
  MARGIN = 1,
  collapse = ""
))

full_details <- df_details_prov[index_keep, ]
rownames(full_details) <- full_details[, "dga_well_code"]

details <- full_details

details$cr2sub_id <- as.numeric(sub("-[0-9Kk]$", "", details$dga_well_code))
details$dga_well_elev <- as.numeric(full_details$dga_well_elev)
details$dga_well_lat <- unname(sapply(full_details$dga_well_lat, deg2dec))
details$dga_well_lon <- unname(sapply(full_details$dga_well_lon, deg2dec))
details$dga_well_utm_north <- as.numeric(full_details$dga_well_utm_north)
details$dga_well_utm_east <- as.numeric(full_details$dga_well_utm_east)

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
  tmp_folder,
  paste0(tag, "_", version, "_metadata.csv")
),
row.names = FALSE
)

# -----------------------------------------------------------------------------
# Section: Merge Time Series By Well
# -----------------------------------------------------------------------------

wells_out <- list.dirs(tmp_folder, recursive = FALSE, full.names = FALSE)
nombres_pozos <- wells_out
codigos_pozos <- rep(NA, length(wells_out))

for (pozo_idx in seq_along(wells_out)) {
  folder_out <- wells_out[pozo_idx]
  files_out <- list.files(file.path(tmp_folder, folder_out), full.names = TRUE)

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
