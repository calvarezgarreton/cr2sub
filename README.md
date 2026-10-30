# CR2SUB: monthly groundwater level database for Chile

CR2SUB compiles groundwater level (GWL) observations from wells of the Chilean Water Bureau (Dirección General de Aguas, DGA), downloaded from the [DGA portal](https://snia.mop.gob.cl/BNAConsultas/reportes), into a homogenized monthly database (1957-2025), plus metadata and topographic, hydroclimatic and soil attributes for each well.

**Developer**: Center for Climate and Resilience Research (CR2, ANID/FONDAP/1523A0002), project ANID/FONDECYT/11240924. Responsible: Camila Alvarez-Garreton. Collaborators: Rodrigo Marinao Rivas, Juan P. Boisier, Diego Dinamarca.

**License**: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

**Citation** (v1.1): Alvarez-Garreton, C., Marinao, R., Boisier, J. P., & Dinamarca, D. (2025). CR2SUB: monthly groundwater level database for Chile (v1.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17591525

## Versions

Each version lives in its own branch with the same folder structure (see [docs/versionado_por_ramas.md](docs/versionado_por_ramas.md)); use `git switch <branch>`.

| Branch | Status |
|---|---|
| `main` | Common documentation, no data |
| `cr2sub_v1.1` | **Obsolete.** Published release (Zenodo); superseded by v1.2 |
| `cr2sub_v1.2` | **Current.** DGA transparency data update (October 2026): single DGA export as input, CAMELS-CL v2 basins |

## Output (`cr2sub/`)

- `cr2sub_v1.2_gwl_mon.csv`: monthly GWL per well (one column per well). Multiple records in a month are averaged.
- `cr2sub_v1.2_gwl_mon_clean.csv`: same series after outlier removal.
- `cr2sub_v1.2_attributes.csv`: one row per well (`cr2sub_id` = DGA well code without the last check digit):

| Columns | Content |
|---|---|
| `dga_*` | Metadata as reported by DGA (code, name, basin, sub-basin, coordinates, elevation) |
| `cr2sub_lat/lon`, `cr2sub_utm_*_h19` | Well coordinates (degrees; UTM 19S, m) |
| `cr2sub_elev`, `cr2sub_slp` | Elevation (m a.s.l.) and slope (degrees) at the well, from FABDEM resampled to 300 m |
| `cr2sub_mean/sd/cv_gwl`, `cr2sub_*_gwl_clean` | GWL mean, standard deviation and coefficient of variation, raw and after outlier removal |
| `cr2sub_in_basin_camels`, `cr2sub_camels_*` | Smallest CAMELS-CL v2 basin containing the well and its mean elevation, slope, precipitation, PET, aridity (PET/P) and snow fraction (CR2MET v2.5, 1980-2010) |
| `cr2sub_in_basin_bna`, `cr2sub_bna_*` | Same for the smallest BNA basin |
| `cr2sub_clsoilmap_*` | Available water content (mm), saturated conductivity (cm/day), bulk density (g/cm³), clay and sand (%) for 0-100 and 100-200 cm, from CLSoilMaps (500 m buffer around the well). NA outside CLSoilMaps coverage |

## Input (`input/`)

| Folder | Content |
|---|---|
| `DGA_GWL_obs_descarga/` | DGA export `AM006T0034616_1.xlsx` (raw observations and well metadata) |
| `camels_cl_basins/` | CAMELS-CL v2 basin boundaries (`.gpkg`, 860 basins), from [CAMELS-CL](https://www.camels.cr2.cl/) |
| `BNA_basins/` | BNA basin boundaries, from the [DGA repository](https://dga.mop.gob.cl/mapoteca-digital/) |
| `cr2sub_basin_lookup_tables/` | Well-to-basin tables and the script that builds them |
| `dem/` | FABDEM elevation and slope, 300 m ([source](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn)) |
| `cr2met/` | CR2MET v2.5 annual precipitation, PET and snowfall ([source](https://ftp.cr2.cl/browse/cr2met/v2.5)) |
| `CLSoilMaps/` | Soil attributes already extracted for each well |
| `other_data/` | Simplified national boundary (BCN-SIIT), used to validate well coordinates |

The CLSoilMaps rasters are large (~16.5 GB) and are not in git. They are only needed when new wells are added: copy `.env.example` to `.env`, set `CR2SUB_HEAVY_DATA` to a local folder and run `Rscript scripts/00_download_clsoilmaps.R` ([source](https://zenodo.org/records/7464210)).

## Running the pipeline

Requires R (packages `readxl`, `zoo`, `terra`, `sf`, `dplyr`, `readr`, `tibble`) and Python with `scripts/functions/python_requirements.txt`. From the repository root:

```bash
python3 scripts/main_pipeline.py
```

| Step | Script | Does |
|---|---|---|
| 01 | `01_consolidate.R` | Reads the DGA export, builds well metadata and monthly series |
| 01b | `01b_extract_clsoilmaps.R` | Adds soil attributes for wells not yet in `input/CLSoilMaps/` |
| 02 | `02_remove_outliers.ipynb` | Outlier removal |
| 03 | `03_process_attributes.R` | Builds `cr2sub_v1.2_attributes.csv` |

If wells are added or removed, rebuild the basin tables before step 03 with `Rscript input/cr2sub_basin_lookup_tables/lookup_table_script.R` (after step 01).
