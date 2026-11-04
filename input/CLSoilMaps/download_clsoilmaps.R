# -----------------------------------------------------------------------------
# Script: 00_download_clsoilmaps.R
# Purpose: Download and unzip the CLSoilMaps rasters used by 01b_extract_clsoilmaps.R
#   into $CR2SUB_HEAVY_DATA/CLSoilMaps. Source: https://zenodo.org/records/7464210
#   (Dinamarca et al., 2023). About 16.5 GB zipped; run only if there are new wells.
# -----------------------------------------------------------------------------

source("../../scripts/functions/functions_env.R")
load_dotenv("../../.env")

heavy <- Sys.getenv("CR2SUB_HEAVY_DATA")
if (!nzchar(heavy)) {
  stop("Define CR2SUB_HEAVY_DATA en .env (ver .env.example).")
}

target <- file.path(heavy, "CLSoilMaps")
dir.create(target, recursive = TRUE, showWarnings = FALSE)
options(timeout = 24 * 3600)

curl <- Sys.which("curl")
if (!nzchar(curl)) {
  stop("No se encontró curl; instálalo antes de descargar CLSoilMaps.")
}

base_url <- "https://zenodo.org/records/7464210/files"
for (name in c("SoilMaps_MEAN", "ROSETTA_MEAN")) {
  zip_file <- file.path(target, paste0(name, ".zip"))
  archive_url <- paste0(base_url, "/", name, ".zip?download=1")

  test_archive <- function() {
    system2(
      "unzip",
      args = c("-tq", shQuote(zip_file)),
      stdout = FALSE,
      stderr = FALSE
    )
  }

  archive_status <- if (file.exists(zip_file)) test_archive() else 1L
  if (archive_status != 0) {
    message("Descargando o reanudando ", name, "...")
    status <- system2(
      curl,
      args = c(
        "--fail",
        "-L",
        "-C", "-",
        "--retry", "5",
        "--retry-all-errors",
        "--progress-bar",
        "-o", zip_file,
        archive_url
      )
    )
    if (status != 0) {
      stop(
        "Falló la descarga de ", name,
        ". El archivo parcial se conserva para reanudarla al volver a ejecutar."
      )
    }
  }

  if (test_archive() != 0) {
    stop("El ZIP de ", name, " está incompleto o corrupto: ", zip_file)
  }

  archive_files <- utils::unzip(zip_file, list = TRUE)$Name
  archive_files <- archive_files[!endsWith(archive_files, "/")]
  missing_files <- archive_files[!file.exists(file.path(target, archive_files))]
  if (length(missing_files) > 0) {
    message("Extrayendo ", length(missing_files), " archivos de ", name, "...")
    status <- system2(
      "unzip",
      args = c("-o", "-q", "-d", shQuote(target), shQuote(zip_file))
    )
    if (status != 0) {
      stop("Falló la extracción de ", name, " (unzip: ", status, ").")
    }
  } else {
    message(name, ": ZIP y archivos extraídos ya están completos.")
  }

  if (any(!file.exists(file.path(target, archive_files)))) {
    stop("La extracción de ", name, " quedó incompleta.")
  }
}
message("Listo: ", target)
