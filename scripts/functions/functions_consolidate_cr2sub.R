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


# -----------------------------------------------------------------------------
# Section: Excel Sheet Processing
# -----------------------------------------------------------------------------
# Function: process_excel_sheet
# Purpose: Processes a sheet from an Excel object and extracts metadata, dates,
#          and data series.
# Arguments:
#   xls: Data frame of the entire Excel sheet.
#   sheet_name: Name of the sheet (for reference).
# Returns:
#   A list with elements:
#     metadata: Character vector of extracted metadata.
#     dates: Date vector of measurement dates.
#     data: Numeric vector of measured values.
process_excel_sheet <- function(xls, sheet_name) {
  # Extract metadata from specified rows and columns
  metadata <- c(xls[c(6,5,7:8),3], xls[6:8, 10], xls[6:7, 15])

  # Identify start of data by locating 'Fecha'
  row_ini <- which(xls[, 1] %in% "Fecha") + 1
  # Identify end of data by locating 'Indicadores'
  row_fin <- which(xls[, 1] %in% "Indicadores") - 1

  # Process data
  header <- xls[row_ini - 1, ]
  xls_clip <- xls[row_ini:row_fin, ]
  # Assign column names using header row
  colnames(xls_clip) <- header

  # Convert Excel serial numbers to Date objects (origin 1899-12-30)
  dates <- as.Date(as.numeric(t(as.matrix(xls_clip[, header %in% "Fecha"]))),
                   origin = as.Date("1899-12-30"))
  # Extract and convert level data to numeric
  data <- as.numeric(t(as.matrix(xls_clip[, header %in% "Nivel (m)"])))

  # Remove entries with missing dates
  valid_data <- !is.na(dates)
  dates <- dates[valid_data]
  data <- data[valid_data]

  # Assemble output list with metadata, dates, and data vectors
  out <- list(
    metadata = metadata,
    dates = dates,
    data = data
  )

  out
}
