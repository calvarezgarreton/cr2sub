pipeline_config <- function(name) {
  key <- paste0("CR2SUB_", toupper(name))
  value <- Sys.getenv(key, unset = "")
  if (!nzchar(value)) {
    stop(sprintf("Pipeline configuration %s is missing.", key), call. = FALSE)
  }
  value
}

pipeline_config_integer <- function(name) {
  value <- suppressWarnings(as.integer(pipeline_config(name)))
  if (is.na(value)) {
    stop(
      sprintf("Pipeline configuration CR2SUB_%s must be an integer.", toupper(name)),
      call. = FALSE
    )
  }
  value
}
