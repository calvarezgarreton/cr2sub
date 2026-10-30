# Re-extract CLSoilMaps attributes for ALL wells with the pipeline's extractor and
# save them next to the values stored in input/CLSoilMaps (computed earlier by D. Dinamarca).
# Run from the repository root (needs .env with CR2SUB_HEAVY_DATA and the rasters):
#   Rscript review_dga_xlsx_migration/scripts/validate_clsoilmaps_extraction.R
suppressMessages({library(terra); library(sf)})
source("scripts/functions/functions_env.R")
source("scripts/functions/functions_clsoilmaps.R")
load_dotenv()

attrs <- read.csv("cr2sub/cr2sub_v1.2_attributes.csv")
wells <- st_as_sf(
  data.frame(well_id = sprintf("%08d", as.integer(attrs$cr2sub_id)),
             lon = attrs$dga_well_lon, lat = attrs$dga_well_lat),
  coords = c("lon", "lat"), crs = 4326
)
wells <- st_transform(st_buffer(st_transform(wells, 32719), 500), 4326)
tmp <- tempfile(fileext = ".geojson")
st_write(wells, tmp, quiet = TRUE)

dirs <- get_clsoilmaps_dirs()
if (is.null(dirs)) stop("Rasters de CLSoilMaps no disponibles (CR2SUB_HEAVY_DATA).")
out <- "review_dga_xlsx_migration/comparison_v1.1_vs_v1.2/clsoilmaps_reextraction.csv"
res <- extract_clsoilmaps_attr(tmp, dirs$phys, dirs$hydro, out, verbose = FALSE)
names(res)[1] <- "cr2sub_id"
write.csv(res, out, row.names = FALSE)
message("Listo: ", nrow(res), " pozos -> ", out)
