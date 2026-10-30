# -----------------------------------------------------------------------------
# Script: 00_download_clsoilmaps.R
# Purpose: Download and unzip the CLSoilMaps rasters used by 01b_extract_clsoilmaps.R
#   into $CR2SUB_HEAVY_DATA/CLSoilMaps. Source: https://zenodo.org/records/7464210
#   (Dinamarca et al., 2023). About 16.5 GB zipped; run only if there are new wells.
# -----------------------------------------------------------------------------

source("scripts/functions/functions_env.R")
load_dotenv()

heavy <- Sys.getenv("CR2SUB_HEAVY_DATA")
if (!nzchar(heavy)) {
  stop("Define CR2SUB_HEAVY_DATA en .env (ver .env.example).")
}

target <- file.path(heavy, "CLSoilMaps")
dir.create(target, recursive = TRUE, showWarnings = FALSE)
options(timeout = 24 * 3600)

base_url <- "https://zenodo.org/records/7464210/files"
for (name in c("SoilMaps_MEAN", "ROSETTA_MEAN")) {
  zip_file <- file.path(target, paste0(name, ".zip"))
  if (dir.exists(file.path(target, name))) {
    message(name, ": ya existe, se omite.")
    next
  }
  if (!file.exists(zip_file)) {
    message("Descargando ", name, "...")
    download.file(
      paste0(base_url, "/", name, ".zip?download=1"), zip_file, mode = "wb"
    )
  }
  unzip(zip_file, exdir = target)
}
message("Listo: ", target)
