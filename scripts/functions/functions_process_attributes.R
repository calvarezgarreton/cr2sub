
create_well_area_join <- function(
  wells,
  polygons,
  buffer_width = 500,
  use_buffer = TRUE,
  grid_extent = c(xmin = -76, xmax = -65, ymin = -57, ymax = -16),
  grid_resolution = 0.01,
  id_priority = c("gauge_id", "COD_CUEN"),
  match_ids = NULL
) {
  ensure_spatvector <- function(obj, label) {
    if (inherits(obj, "SpatVector")) {
      return(obj)
    }
    if (is.character(obj) && length(obj) == 1) {
      return(terra::vect(obj))
    }
    stop(label, " must be a SpatVector or a single file path.")
  }

  wells_vect <- ensure_spatvector(wells, "`wells`")
  polygons_vect <- ensure_spatvector(polygons, "`polygons`")

  polygons_vect$area <- terra::expanse(polygons_vect)

  ensure_epsg4326 <- function(vec, label) {
    crs_text <- terra::crs(vec)
    if (!nzchar(crs_text)) {
      stop(label, " must have a defined CRS.")
    }
    if (!isTRUE(terra::same.crs(vec, "EPSG:4326"))) {
      return(terra::project(vec, "EPSG:4326"))
    }
    vec
  }

  wells_vect <- ensure_epsg4326(wells_vect, "`wells`")
  polygons_vect <- ensure_epsg4326(polygons_vect, "`polygons`")

  if (!"cr2sub_id" %in% names(wells_vect)) {
    stop("`wells` must include a `cr2sub_id` column.")
  }

  id_candidates <- intersect(id_priority, names(polygons_vect))
  if (!length(id_candidates)) {
    stop("`polygons` must include one of the columns: ",
      paste(id_priority, collapse = ", ")
    )
  }
  id_col <- id_candidates[1]

  extract_by_area <- function(polys, wells_vect, grid, id_col) {
    total_polys <- nrow(polys)
    polygons_list <- lapply(seq_len(total_polys), function(i) polys[i, ])
    ids <- sapply(
      polygons_list,
      function(x) as.data.frame(x)[, id_col]
    )

    polygons_area <- lapply(
      polygons_list,
      terra::rasterize,
      y = grid,
      field = "area"
    )
    polygons_stack <- do.call(c, polygons_area)
    names(polygons_stack) <- ids

    extract_values <- terra::extract(polygons_stack, wells_vect, ID = FALSE)
    extract_matrix <- t(as.matrix(extract_values))
    rownames(extract_matrix) <- names(polygons_stack)
    colnames(extract_matrix) <- wells_vect$cr2sub_id

    id_vec <- rep(rownames(extract_matrix), times = ncol(extract_matrix))
    well_vec <- rep(colnames(extract_matrix), each = nrow(extract_matrix))

    full_join <- data.frame(
      cr2sub_id = well_vec,
      id = id_vec,
      area = round(as.numeric(extract_matrix), 4)
    )
    colnames(full_join) <- c("cr2sub_id", id_col, "area")
    full_join
  }

  grid <- terra::rast(
    xmin = grid_extent["xmin"],
    xmax = grid_extent["xmax"],
    ymin = grid_extent["ymin"],
    ymax = grid_extent["ymax"],
    res = grid_resolution,
    crs = "EPSG:4326"
  ) |>
    terra::crop(terra::ext(polygons_vect), snap = "out")

  full_no_buffer <- extract_by_area(
    polygons_vect,
    wells_vect,
    grid,
    id_col
  )

  no_match_ids <- unique(full_no_buffer$cr2sub_id[is.na(full_no_buffer$area)])
  wells_unmatched <- wells_vect[wells_vect$cr2sub_id %in% no_match_ids, ]

  inner_buffered <- full_no_buffer[FALSE, ]
  if (use_buffer && nrow(wells_unmatched)) {
    polygons_buffered <- terra::buffer(polygons_vect, width = buffer_width)
    full_buffered <- extract_by_area(
      polygons_buffered,
      wells_unmatched,
      grid,
      id_col
    )
    inner_buffered <- full_buffered[!is.na(full_buffered$area), ]
  }

  inner_no_buffer <- full_no_buffer[!is.na(full_no_buffer$area), ]
  inner_join <- rbind(inner_no_buffer, inner_buffered)
  inner_join <- inner_join[order(inner_join$cr2sub_id), ]
  inner_join <- inner_join[!is.na(inner_join$area), ]

  if (nrow(inner_join)) {
    inner_join <- inner_join[
      order(inner_join$cr2sub_id, -inner_join$area),
    ]
    inner_join <- inner_join[!duplicated(inner_join$cr2sub_id), ]
  }

  if (is.null(match_ids)) {
    match_ids <- wells_vect$cr2sub_id
  }

  match_ids_chr <- as.character(match_ids)
  lookup <- stats::setNames(inner_join[[id_col]], inner_join$cr2sub_id)
  match_vector <- lookup[match_ids_chr]
  match_vector <- utils::type.convert(as.character(match_vector), as.is = TRUE)
  match_vector
}
