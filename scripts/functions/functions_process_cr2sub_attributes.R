
which_fun <- function(x, fun) {
  value <- sapply(list(x), FUN = fun)
  return(which(x == value))
}


select_basin_from_join_table <- function(well, join_table, fun) {
  col_candidates <- intersect(c("gauge_id", "bna_id"), colnames(join_table))
  if (!length(col_candidates)) {
    stop("`join_table` must include either `gauge_id` or `COD_CUEN`.")
  }
  col_sel <- col_candidates[1]

  well_filter <- join_table[, "cr2sub_id"] == well
  well_basins <- join_table[well_filter, , drop = FALSE]

  index_sel <- which_fun(well_basins[, "area"], fun = fun) # criterio según fun
  id <- well_basins[index_sel, col_sel] |> as.character()
  id[1]
}
