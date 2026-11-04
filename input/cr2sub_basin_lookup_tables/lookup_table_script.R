# -----------------------------------------------------------------------------
# Script: lookup_table_script.R
# Purpose: Build the lookup tables that assign each CR2SUB well to the
#   CAMELS-CL v2 and BNA basins that contain it. Run from the repository
#   root, after 01_consolidate.R (needs .tmp/<tag>_<version>_metadata.csv):
#     Rscript input/cr2sub_basin_lookup_tables/lookup_table_script.R
#   Rerun only when wells are added or removed.
# Output: input/cr2sub_basin_lookup_tables/<tag>_<version>_join_table_with_*.csv
#   (columns: cr2sub_id, basin id, basin area in km2)
# -----------------------------------------------------------------------------

library(terra)

tag <- "cr2sub"
version <- "v1.2"
out_dir <- "input/cr2sub_basin_lookup_tables"
metadata_file <- file.path(".tmp", paste0(tag, "_", version, "_metadata.csv"))

if (!file.exists(metadata_file)) {
  stop("No se encontró la metadata en .tmp; ejecuta antes 01_consolidate.R.")
}

# Basin layers: file, name of the id column, and how the id is derived
basin_types <- list(
  camels_cl_v2_basins = list(
    file = "input/camels_cl_basins/camels_cl_v2_boundaries.gpkg",
    id_col = "gauge_id",
    get_id = function(v) as.character(as.data.frame(v)[, "gauge_id"])
  ),
  bna_basins = list(
    file = "input/BNA_basins/Cuencas_BNA.shp",
    id_col = "bna_id",
    get_id = function(v) as.character(as.numeric(as.data.frame(v)[, "COD_CUEN"]))
  )
)

# Wells (cr2sub_id is 8 digits; restore leading zeros lost when stored as number)
metadata <- read.csv(metadata_file)
wells_vect <- vect(
  data.frame(
    well_id = sprintf("%08d", as.integer(metadata$cr2sub_id)),
    lon = metadata$dga_well_lon,
    lat = metadata$dga_well_lat
  ),
  geom = c("lon", "lat"), crs = "EPSG:4326"
)

# Basins that intersect each well (point or buffer), with their area (km2). Vector
# intersection (low memory).
match_wells <- function(wells, basins) {
  pairs <- terra::relate(wells, basins, relation = "intersects", pairs = TRUE)
  data.frame(
    well_id = wells$well_id[pairs[, 1]],
    id = basins$id[pairs[, 2]],
    area = round(basins$area[pairs[, 2]], 4)
  )
}

for (type in names(basin_types)) {
  cfg <- basin_types[[type]]
  basins <- vect(cfg$file) |> project("EPSG:4326")
  basins$id <- cfg$get_id(basins)
  basins$area <- expanse(basins, unit = "km")

  # Basins that contain the well; wells that fall in none are retried within 500 m
  inner_join <- match_wells(wells_vect, basins)

  unmatched <- wells_vect[!wells_vect$well_id %in% inner_join$well_id, ]
  if (nrow(unmatched) > 0) {
    unmatched_buffer <- buffer(unmatched, width = 500)
    unmatched_buffer$well_id <- unmatched$well_id
    inner_join <- rbind(inner_join, match_wells(unmatched_buffer, basins))
  }
  inner_join <- unique(inner_join[order(inner_join$well_id, inner_join$area), ])
  colnames(inner_join) <- c("cr2sub_id", cfg$id_col, "area")

  write.csv(
    inner_join,
    file.path(out_dir, paste0(tag, "_", version, "_join_table_with_", type, ".csv")),
    row.names = FALSE
  )
  message("[cr2sub] ", type, ": ", length(unique(inner_join$cr2sub_id)), " pozos asignados")
}
