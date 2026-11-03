# -----------------------------------------------------------------------------
# Script: 01b_extract_clsoilmaps.R
# Purpose: Extract CLSoilMaps soil attributes for ALL wells from the original
#   rasters (500 m buffer around each well) and write
#   input/CLSoilMaps/cr2sub_CLSoilMaps_data.csv. Nothing is inherited from
#   previous versions of that file.
# Requires: .env with CR2SUB_HEAVY_DATA (see .env.example) and the rasters
#   downloaded with scripts/00_download_clsoilmaps.R.
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

dirs <- get_clsoilmaps_dirs()
if (is.null(dirs)) {
  stop(
    "Faltan los rasters de CLSoilMaps (CR2SUB_HEAVY_DATA en .env). ",
    "Descárgalos con scripts/00_download_clsoilmaps.R."
  )
}

metadata <- read.csv(metadata_file)
# cr2sub_id has 8 digits; restore the leading zero lost when stored as numeric
metadata$cr2sub_id <- sprintf("%08d", as.integer(metadata$cr2sub_id))
message("[cr2sub] CLSoilMaps: extrayendo ", nrow(metadata), " pozos desde los rasters.")

# 500 m buffer around each well (metric CRS, then back to lon/lat)
wells <- st_as_sf(
  data.frame(well_id = metadata$cr2sub_id,
             lon = metadata$dga_well_lon, lat = metadata$dga_well_lat),
  coords = c("lon", "lat"), crs = 4326
)
wells <- st_transform(st_buffer(st_transform(wells, 32719), buffer_m), 4326)
wells_file <- file.path(".tmp", "wells_buffer_500m.geojson")
st_write(wells, wells_file, delete_dsn = TRUE, quiet = TRUE)

result <- extract_clsoilmaps_attr(
  wells.shp = wells_file,
  phys.dir = dirs$phys,
  hydro.dir = dirs$hydro,
  out = file.path(".tmp", "clsoilmaps_all_wells.csv")
)
names(result)[names(result) == "well_id"] <- "cr2sub_id"
result$cr2sub_id <- as.character(result$cr2sub_id)
result <- as.data.frame(result)
result[is.na(result)] <- NA  # NaN (all cells NA in the buffer) -> NA
result <- result[order(result$cr2sub_id), ]

write.csv(result, soil_file, row.names = FALSE)
message("[cr2sub] CLSoilMaps: ", nrow(result), " pozos en ", soil_file,
        " (", sum(!complete.cases(result)), " con NA)")
