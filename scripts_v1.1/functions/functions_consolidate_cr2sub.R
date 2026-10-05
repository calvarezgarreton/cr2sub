# -----------------------------------------------------------------------------
# Section: Directory Management
# -----------------------------------------------------------------------------
# Function: create_dir_if_not_exists
# Purpose: Creates a directory at the given path if it does not already exist,
#          including any necessary parent directories.
# Arguments:
#   dir_path: Character string specifying the path of the directory to create.
# Returns:
#   NULL. Creates directory as side effect.
create_dir_if_not_exists <- function(dir_path) {
  # Check if the directory exists
  if (!dir.exists(dir_path)) {
    # Create the directory recursively
    dir.create(dir_path, recursive = TRUE)
  }
}

# -----------------------------------------------------------------------------
# Section: Coordinate Conversion
# -----------------------------------------------------------------------------
# Function: deg2dec
# Purpose: Converts coordinates from degrees-minutes-seconds format to decimal
#          degrees.
# Arguments:
#   x: Character string in the format "D° M' S\"" (degrees, minutes, seconds).
# Returns:
#   Numeric value representing decimal degrees, rounded to 6 decimal places.
deg2dec <- function(x) {
  # Remove degree and minute symbols
  y <- gsub("'", "", gsub("°", "", x))
  # Split cleaned string into components: degrees, minutes, seconds
  z <- as.numeric(strsplit(y, " ")[[1]])
  # Calculate decimal degrees and round to 6 decimals
  out <- -1 * round(z[1] + z[2] / 60 + z[3] / 3600, 6)
  out
}

# Function: deg2dec_compact
# Purpose: Converts coordinates from the DGA portal's compact degrees-minutes-
#          seconds format ("DDDMMSS", zero-padded, e.g. "0201212") to decimal
#          degrees.
# Arguments:
#   x: Character string with 7 digits: 3 for degrees, 2 for minutes, 2 for
#      seconds (trailing whitespace allowed).
# Returns:
#   Numeric value representing decimal degrees, rounded to 6 decimal places.
deg2dec_compact <- function(x) {
  y <- trimws(x)
  deg <- as.numeric(substr(y, 1, 3))
  min <- as.numeric(substr(y, 4, 5))
  sec <- as.numeric(substr(y, 6, 7))
  out <- -1 * round(deg + min / 60 + sec / 3600, 6)
  out
}

