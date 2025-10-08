rm(list = ls())

# 1. Configurar directorio de trabajo
wd <- "/Users/rmarinao/Library/CloudStorage/Dropbox/RM_fdcyt_cag"
setwd(wd)

# 2. Cargar librerías necesarias
library(terra)

# 3. Definir subdirectorio de salida y tipo de unidad territorial
type <- "Communes"

# 4. Asociar tipo con nombre de columna identificadora
basin_name_id <- c(
  "camels_cl_v2021_basins" = "gauge_",
  "bna_basins" = "bna_id",
  "Communes" = "id"
)[type]

# 5. Leer vectores de pozos y polígonos SIN buffer
wells_vect <- vect(
  "OUT/attribute_data_frames/gwl/gwl_dga_attributes_epsg4326.gpkg"
)

somebasins_raw <- vect(
  paste0("IN/polygons/", type, "_boundaries_epsg4326.gpkg")
)

# Función auxiliar para extraer join por área
extract_by_area <- function(somebasins, wells_vect, grid, basin_name_id) {
  to <- nrow(somebasins)
  somebasins_list <- lapply(1:to, function(i) somebasins[i, ])
  id <- sapply(somebasins_list, function(x) {
    as.data.frame(x)[, basin_name_id]
  })

  somebasins_area <- lapply(
    somebasins_list, rasterize, y = grid, field = "area"
  )
  somebasins_stack <- do.call(c, somebasins_area)
  names(somebasins_stack) <- id

  extract_matrix <- t(
    terra::extract(somebasins_stack, wells_vect, ID = FALSE)
  )
  colnames(extract_matrix) <- wells_vect$well_id

  id_vec <- array(
    rownames(extract_matrix),
    dim = dim(extract_matrix)
  ) |> as.character()

  wellcode <- t(array(
    colnames(extract_matrix),
    dim = rev(dim(extract_matrix))
  )) |> as.character()

  full_join <- data.frame(
    well_id = wellcode,
    id = id_vec,
    area = round(as.numeric(extract_matrix), 4)
  )
  colnames(full_join) <- c("well_id", basin_name_id, "area")
  full_join
}

# 6. Crear grilla base
grid <- terra::rast(
  xmin = -76, xmax = -65, ymin = -57, ymax = -16,
  res = 0.01, crs = "EPSG:4326"
) |> crop(ext(somebasins_raw), snap = "out")

# 7. Extraer join SIN buffer
full_no_buffer <- extract_by_area(somebasins_raw, wells_vect, grid,
                                  basin_name_id)

# 8. Identificar pozos sin asignación
no_match_ids <- unique(full_no_buffer$well_id[is.na(full_no_buffer$area)])
wells_unmatched <- wells_vect[wells_vect$well_id %in% no_match_ids, ]

# 9. Extraer join CON buffer solo para pozos sin asignación
somebasins_buffered <- buffer(somebasins_raw, width = 500)
full_buffered <- extract_by_area(somebasins_buffered, wells_unmatched,
                                 grid, basin_name_id)

# 10. Unir ambos resultados
inner_no_buffer <- full_no_buffer[!is.na(full_no_buffer$area), ]
inner_buffered <- full_buffered[!is.na(full_buffered$area), ]
inner_join <- rbind(inner_no_buffer, inner_buffered)
inner_join <- inner_join[order(inner_join$well_id), ]

# 11. Exportar resultado
if (!dir.exists(subdir)) dir.create(subdir)

write.csv(
  inner_join,
  paste0("cr2sub_v1.1_join_table_with_", type, ".csv"),
  row.names = FALSE
)