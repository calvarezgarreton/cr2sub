# -----------------------------------------------------------------------------
# Script: 01b_extract_clsoilmaps.R
# Purpose: Create or update the version-specific CLSoilMaps CSV configured in
#   scripts/main_pipeline.py. A missing CSV is generated for all wells;
#   existing CSV rows are preserved and only new wells are extracted.
# Requires: .env with CR2SUB_HEAVY_DATA (see .env.example) and the rasters
#   downloaded with input/CLSoilMaps/download_clsoilmaps.R
# -----------------------------------------------------------------------------

library(terra)
library(sf)

source("scripts/functions/functions_pipeline_config.R")
source("scripts/functions/functions_env.R")
source("scripts/functions/functions_clsoilmaps.R")

load_dotenv()

metadata_file <- pipeline_config("metadata_file")
soil_file <- pipeline_config("clsoil_data_file")
buffer_m <- pipeline_config_integer("clsoil_buffer_m")
wells_file <- pipeline_config("wells_buffer_file")
new_wells_file <- pipeline_config("clsoil_new_wells_file")

if (!file.exists(metadata_file)) {
  stop(
    "No se encontró la metadata en ", dirname(metadata_file),
    "; ejecuta antes 01_consolidate.R."
  )
}

metadata <- read.csv(metadata_file)
# cr2sub_id has 8 digits; restore the leading zero lost when stored as numeric
metadata$cr2sub_id <- sprintf("%08d", as.integer(metadata$cr2sub_id))
soil <- if (file.exists(soil_file)) {
  read.csv(soil_file, colClasses = c(cr2sub_id = "character"))
} else {
  NULL
}

existing_ids <- if (is.null(soil)) character(0) else soil$cr2sub_id
new_wells <- metadata[!metadata$cr2sub_id %in% existing_ids, ]
if (nrow(new_wells) == 0) {
  message("[cr2sub] CLSoilMaps: no hay pozos nuevos; se usa ", soil_file)
  quit(save = "no", status = 0)
}

action <- if (is.null(soil)) "generando" else "actualizando"
message("[cr2sub] CLSoilMaps: ", action, " datos para ",
        nrow(new_wells), " pozos.")
has_coordinates <- complete.cases(
  new_wells[, c("dga_well_lon", "dga_well_lat")]
)
valid_new_wells <- new_wells[has_coordinates, ]
new_soil <- if (is.null(soil)) NULL else soil[FALSE, ]

if (nrow(valid_new_wells) > 0) {
  dirs <- get_clsoilmaps_dirs()
  if (is.null(dirs)) {
    stop(
      "Faltan los rasters de CLSoilMaps para extraer los pozos nuevos. ",
      "Descárgalos con input/CLSoilMaps/download_clsoilmaps.R."
    )
  }

  wells <- st_as_sf(
    data.frame(well_id = valid_new_wells$cr2sub_id,
               lon = valid_new_wells$dga_well_lon,
               lat = valid_new_wells$dga_well_lat),
    coords = c("lon", "lat"), crs = 4326
  )
  wells <- st_transform(st_buffer(st_transform(wells, 32719), buffer_m), 4326)
  st_write(wells, wells_file, delete_dsn = TRUE, quiet = TRUE)

  new_soil <- extract_clsoilmaps_attr(
    wells.shp = wells_file,
    phys.dir = dirs$phys,
    hydro.dir = dirs$hydro,
    out = new_wells_file
  )
  names(new_soil)[names(new_soil) == "well_id"] <- "cr2sub_id"
  new_soil$cr2sub_id <- as.character(new_soil$cr2sub_id)
  new_soil <- as.data.frame(new_soil)
  new_soil[is.na(new_soil)] <- NA
}

if (is.null(soil)) {
  if (nrow(valid_new_wells) == 0) {
    stop("No hay pozos con coordenadas válidas para crear el CSV de CLSoilMaps.")
  }
  soil <- new_soil
  new_soil <- soil[FALSE, ]
}

missing_coordinates <- new_wells$cr2sub_id[!has_coordinates]
if (length(missing_coordinates) > 0) {
  missing_soil <- data.frame(
    cr2sub_id = missing_coordinates,
    matrix(
      NA_real_,
      nrow = length(missing_coordinates),
      ncol = ncol(soil) - 1,
      dimnames = list(NULL, names(soil)[-1])
    ),
    check.names = FALSE
  )
  new_soil <- rbind(new_soil, missing_soil)
}

soil <- rbind(soil, new_soil[, names(soil), drop = FALSE])
soil <- soil[order(soil$cr2sub_id), ]
write.csv(soil, soil_file, row.names = FALSE, na = "NA")
message("[cr2sub] CLSoilMaps: ", nrow(soil), " pozos en ", soil_file)
