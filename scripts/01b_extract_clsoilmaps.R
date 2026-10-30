# -----------------------------------------------------------------------------
# Script: 01b_extract_clsoilmaps.R
# Purpose: Add CLSoilMaps soil attributes for wells that are not yet in
#   input/CLSoilMaps/cr2sub_CLSoilMaps_data.csv. Existing wells are kept, so
#   the heavy CLSoilMaps rasters are only needed when new wells appear.
# Requires: .env with CR2SUB_HEAVY_DATA (see .env.example) and the rasters
#   downloaded with scripts/00_download_clsoilmaps.R, only if there are new wells.
# -----------------------------------------------------------------------------

library(terra)
library(sf)

source("scripts/functions/functions_env.R")
source("scripts/functions/functions_clsoilmaps.R")

load_dotenv()

tag <- "cr2sub"
version <- "v1.2"
metadata_file <- file.path(".tmp", paste0(tag, "_", version, "_metadata.csv"))
soil_file <- "input/CLSoilMaps/cr2sub_CLSoilMaps_data.csv"
buffer_m <- 500

if (!file.exists(metadata_file)) {
  stop("No se encontró la metadata en .tmp; ejecuta antes 01_consolidate.R.")
}

metadata <- read.csv(metadata_file, colClasses = c(cr2sub_id = "character"))
soil <- read.csv(soil_file, colClasses = c(cr2sub_id = "character"))

new_wells <- metadata[!metadata$cr2sub_id %in% soil$cr2sub_id, ]
if (nrow(new_wells) == 0) {
  message("[cr2sub] CLSoilMaps: no hay pozos nuevos; se usa ", soil_file)
  quit(save = "no", status = 0)
}

message("[cr2sub] CLSoilMaps: ", nrow(new_wells), " pozos nuevos.")
dirs <- get_clsoilmaps_dirs()
if (is.null(dirs)) {
  warning(
    "Faltan los rasters de CLSoilMaps (CR2SUB_HEAVY_DATA en .env). ",
    "Los pozos nuevos quedarán con NA. Descárgalos con scripts/00_download_clsoilmaps.R.",
    call. = FALSE
  )
  quit(save = "no", status = 0)
}

# 500 m buffer around each well (metric CRS, then back to lon/lat)
wells <- st_as_sf(
  data.frame(well_id = new_wells$cr2sub_id,
             lon = new_wells$dga_well_lon, lat = new_wells$dga_well_lat),
  coords = c("lon", "lat"), crs = 4326
)
wells <- st_transform(st_buffer(st_transform(wells, 32719), buffer_m), 4326)
wells_file <- file.path(".tmp", "wells_buffer_500m.geojson")
st_write(wells, wells_file, delete_dsn = TRUE, quiet = TRUE)

result <- extract_clsoilmaps_attr(
  wells.shp = wells_file,
  phys.dir = dirs$phys,
  hydro.dir = dirs$hydro,
  out = file.path(".tmp", "clsoilmaps_new_wells.csv")
)
names(result)[names(result) == "well_id"] <- "cr2sub_id"
result$cr2sub_id <- as.character(result$cr2sub_id)

soil <- rbind(soil, result[, names(soil)])
soil <- soil[order(soil$cr2sub_id), ]
write.csv(soil, soil_file, row.names = FALSE)
message("[cr2sub] CLSoilMaps: ", nrow(soil), " pozos en ", soil_file)
