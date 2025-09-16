rm(list = ls())

wd <- "/Users/rmarinao/Library/CloudStorage/Dropbox/RM_fdcyt_cag"
setwd(wd)

library(terra)
library(hydroTSM)

markham_index_zoo <- function(z) {
  if (!zoo::is.zoo(z)) stop("'z' debe ser 'zoo'")
  S <- numeric(ncol(z)); names(S) <- colnames(z)
  months <- as.integer(format(index(z), "%m"))
  for (j in seq_len(ncol(z))) {
    x <- coredata(z[, j])
    mclim <- tapply(x, months, mean, na.rm = TRUE)
    full  <- rep(NA_real_, 12); names(full) <- 1:12
    full[names(mclim)] <- mclim; mclim <- full
    mu <- mean(mclim, na.rm = TRUE)
    S[j] <- 0.5 * sum(abs(mclim - mu), na.rm = TRUE) /
      sum(mclim, na.rm = TRUE)
  }
  return(S)
}


# -----------------------------------------------------------
# decadal_slope_zoo_simple
# Pendiente lineal por década para cada columna de un zoo mensual
#   z      : zoo mensual (índice Date)
#   start  : (opcional) "YYYY-MM-DD" o Date; inicio de la ventana
#   end    : (opcional) idem; fin de la ventana
#   dpy    : días por año (default 365.25)
# -----------------------------------------------------------
decadal_slope_zoo <- function(z, start = NULL, end = NULL, dpy = 365.25) {
  if (!requireNamespace("zoo", quietly = TRUE) || !zoo::is.zoo(z))
    stop("'z' debe ser un objeto 'zoo' con fechas.")
  
  if (!is.null(start) || !is.null(end))
    z <- window(z, start = start, end = end)
  
  if (nrow(z) < 2)
    stop("Menos de dos observaciones tras aplicar la ventana.")
  
  t_days <- as.numeric(index(z))        # días desde 1970-01-01
  factor  <- dpy * 10                   # días por década
  
  n_col  <- ncol(z)
  slope  <- numeric(n_col)
  names(slope) <- colnames(z)
  
  for (j in seq_len(n_col)) {
    y  <- coredata(z[, j])
    ok <- !is.na(y)
    if (sum(ok) < 2) { slope[j] <- NA_real_; next }
    
    slope[j] <- coef(lm(y[ok] ~ t_days[ok]))[2] * factor
  }
  return(slope)
}



subdir <- "OUT/attribute_data_frames/gwl" 
source("scripts/fun/fun_build_attribute_data_frames.R") # some functions


gwl <- vect("IN/gwl/gwl_m_dga_mon_info_1957_2024_epsg4326.gpkg")


gwl_ts0  <- read.csv.zoo("IN/gwl/gwl_m_dga_mon_ts_1957_2024.csv",
                       format = "%Y-%m-%d", check.names = FALSE)


gwl_ts <- gwl_ts0[,as.character(gwl$well_id)]

gwl_avg <- apply(coredata(gwl_ts), 2, mean, na.rm = TRUE)
gwl_sd <- apply(coredata(gwl_ts), 2, sd, na.rm = TRUE)
gwl_p10 <- apply(coredata(gwl_ts), 2, quantile, probs = 0.1, na.rm = TRUE)
gwl_p25 <- apply(coredata(gwl_ts), 2, quantile, probs = 0.25, na.rm = TRUE)
gwl_p50 <- apply(coredata(gwl_ts), 2, quantile, probs = 0.5, na.rm = TRUE)
gwl_p75 <- apply(coredata(gwl_ts), 2, quantile, probs = 0.75, na.rm = TRUE)
gwl_p90 <- apply(coredata(gwl_ts), 2, quantile, probs = 0.9, na.rm = TRUE)

# gwl_Sea <- markham_index_zoo(gwl_ts)
gwl_trend <- decadal_slope_zoo(gwl_ts, start = "1980-01-01")


any(names(gwl_p90) != names(gwl_trend))

gwl$gwl_avg <- round(gwl_avg,5)
gwl$gwl_sd  <- round(gwl_sd ,5)
gwl$gwl_cv  <- gwl$gwl_sd/gwl$gwl_avg
gwl$gwl_p10 <- round(gwl_p10,5)
gwl$gwl_p25 <- round(gwl_p25,5)
gwl$gwl_p50 <- round(gwl_p50,5)
gwl$gwl_p75 <- round(gwl_p75,5)
gwl$gwl_p90 <- round(gwl_p90,5)
# gwl$gwl_Sea <- gwl_Sea
gwl$gwl_trend <- gwl_trend

### topo

elev_raster <- rast("IN/topo/elevation_srtm_v4.1_Chile_2000_1200m_epsg4326.nc")
slp_raster <- rast("IN/topo/slope_srtm_v4.1_Chile_2000_1200m_epsg4326.nc")



elev <- terra::extract(elev_raster, y = gwl)[,2]

slp <- terra::extract(slp_raster, y = gwl)[,2]

gwl$gauge_elev <- round(elev,5)
gwl$gauge_slp <- round(slp,5)

## meteo

pr <- rast("IN/meteo/cr2met_v2p5/ann/pr_mm_cr2met_2p5best_ann_1960_2024_0p05deg_epsg4326.nc") |>
        app(fun = mean, na.rm = TRUE)
et0 <- rast("IN/meteo/cr2met_v2p5/ann/et0_mm_cr2met_2p5best_ann_1960_2024_0p05deg_epsg4326.nc") |>
        app(fun = mean, na.rm = TRUE)
snow <- rast("IN/meteo/cr2met_v2p5/ann/snow_mm_cr2met_2p5best_ann_1960_2024_0p05deg_epsg4326.nc") |>
           app(fun = mean, na.rm = TRUE)


snowf_rast <- snow/pr
aridity_rast <- et0/pr

terra::plot(snowf_rast)
terra::plot(aridity_rast)

pr_spring <- rast("IN/meteo/cr2met_v2p5/season/pr_mm_cr2met_2p5best_south_spring_1960_2024_0p05deg_epsg4326.nc") |>
  app(fun = mean, na.rm = TRUE)

pr_summer <- rast("IN/meteo/cr2met_v2p5/season/pr_mm_cr2met_2p5best_south_summer_1960_2025_0p05deg_epsg4326.nc") |>
  app(fun = mean, na.rm = TRUE)

pr_autumn <- rast("IN/meteo/cr2met_v2p5/season/pr_mm_cr2met_2p5best_south_autumn_1960_2024_0p05deg_epsg4326.nc") |>
  app(fun = mean, na.rm = TRUE)

pr_winter <- rast("IN/meteo/cr2met_v2p5/season/pr_mm_cr2met_2p5best_south_winter_1960_2024_0p05deg_epsg4326.nc") |>
  app(fun = mean, na.rm = TRUE)

pr_yr_rast <- pr
pr_om_rast <- pr_spring + pr_summer
pr_as_rast <- pr_autumn + pr_winter





gwl_pr_yr <- terra::extract(pr_yr_rast, y = gwl)[,2]
gwl_pr_om <- terra::extract(pr_om_rast, y = gwl)[,2]
gwl_pr_as <- terra::extract(pr_as_rast, y = gwl)[,2]

gwl_aridity <- terra::extract(aridity_rast, y = gwl)[,2]
gwl_snowf   <- terra::extract(snowf_rast, y = gwl)[,2]


gwl$gauge_pr_yr <- round(gwl_pr_yr,5)
gwl$gauge_pr_om <- round(gwl_pr_om,5)
gwl$gauge_pr_as <- round(gwl_pr_as,5)

gwl$gauge_aridity <- round(gwl_aridity,5)
gwl$gauge_snowf   <- round(gwl_snowf,5)

gwl_df <- as.data.frame(gwl)

if(!dir.exists(subdir)) dir.create(subdir)

writeVector(gwl,  paste0(subdir, "/gwl_dga_attributes_epsg4326.gpkg"), overwrite = TRUE)
write.csv(gwl_df, paste0(subdir, "/gwl_dga_attributes.csv"), row.names = FALSE)

