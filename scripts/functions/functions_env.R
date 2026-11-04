# -----------------------------------------------------------------------------
# Section: Local environment (.env)
# -----------------------------------------------------------------------------
# Function: load_dotenv
# Purpose: Reads KEY=VALUE lines from a .env file into environment variables.
#          Existing environment variables are not overwritten.
# Arguments:
#   path: Path to the .env file.
# Returns:
#   Invisibly, TRUE if the file was read, FALSE if it does not exist.
load_dotenv <- function(path = ".env") {
  if (!file.exists(path)) {
    return(invisible(FALSE))
  }
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#") & grepl("=", lines)]
  keys <- trimws(sub("=.*$", "", lines))
  values <- trimws(sub("^[^=]*=", "", lines))
  values <- gsub("^[\"']|[\"']$", "", values)
  for (i in seq_along(keys)) {
    if (!nzchar(Sys.getenv(keys[i]))) {
      do.call(Sys.setenv, setNames(list(values[i]), keys[i]))
    }
  }
  invisible(TRUE)
}

# Function: get_clsoilmaps_dirs
# Purpose: Resolves the CLSoilMaps raster folders from CR2SUB_HEAVY_DATA.
# Returns:
#   Named list with `phys` and `hydro`, or NULL if the folders are missing.
get_clsoilmaps_dirs <- function() {
  heavy <- Sys.getenv("CR2SUB_HEAVY_DATA")
  if (!nzchar(heavy)) {
    return(NULL)
  }
  dirs <- list(
    phys = file.path(heavy, "CLSoilMaps", "SoilMaps_MEAN"),
    hydro = file.path(heavy, "CLSoilMaps", "ROSETTA_MEAN")
  )
  if (!all(dir.exists(unlist(dirs)))) {
    return(NULL)
  }
  dirs
}
