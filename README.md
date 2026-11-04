# CR2SUB: monthly groundwater level database for Chile

CR2SUB compiles groundwater level (GWL) observations from wells of the Chilean Water Bureau (Dirección General de Aguas, DGA), downloaded from the [DGA portal](https://snia.mop.gob.cl/BNAConsultas/reportes), into a homogenized monthly database (1957-August 2026), plus metadata and topographic, hydroclimatic and soil attributes for each well.

**Developer**: Center for Climate and Resilience Research (CR2, ANID/FONDAP/1523A0002), project ANID/FONDECYT/11240924, project ANID IDeA TA TA25I10030. Responsible: Camila Alvarez-Garreton. Collaborators: Rodrigo Marinao Rivas, Juan P. Boisier, Diego Dinamarca.

**License**: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

**Citation**: TBD.

## Versions

Each version lives in its own branch with the same folder structure (see [docs/versionado_por_ramas.md](docs/versionado_por_ramas.md)); use `git switch <branch>`.

| Branch | Status |
|---|---|
| `main` | Common documentation, no data |
| `cr2sub_v1.1` | **Obsolete.** Published release (Zenodo); superseded by v1.2 |
| `cr2sub_v1.2` | **Current.** DGA transparency data update (October 2026): single DGA export as input, CAMELS-CL v2 basins. Data requested under the Chilean Transparency Law (see [Data source](#data-source-v12)) |

## What's New in v1.2

- **New DGA source**: single export obtained through the Transparency Law (see [Data source](#data-source-v12)). Series now extend to August 2026 (v1.1 ended in December 2024).
- **1,133 wells** (v1.1: 1,137). Four wells of the Taltal basin are not in the new export.
- **CAMELS-CL v2** basins (860) replace CAMELS-CL v1 (516); basin assignment now uses the basin that contains the well.
- **CLSoilMaps extraction** is part of the pipeline and requires the heavy rasters (~16.5 GB).
- **Self-contained pipeline**: local paths go in `.env`, and the review of v1.1 vs v1.2 is in `review_dga_xlsx_migration/`.

## Data source (v1.2)

The v1.2 groundwater level data were requested from the DGA under the Chilean Transparency Law (Ley de Transparencia, Ley 20.285) through the [Portal de Transparencia](https://www.portaltransparencia.cl/PortalPdT/ingreso-sai-v2). Request N° AM006T0034616 was received on 2026-09-21 (static water levels in wells, January 1950 to August 2026, for 1,137 stations) and answered by the DGA on 2026-10-01. The response is the single export `input/DGA_GWL_obs_descarga/AM006T0034616_1.xlsx`. The DGA reply letter is in `input/DGA_GWL_obs_descarga/transparency_request_AM006T0034616_dga_reply.pdf`.

## Output (`cr2sub/`)

- `cr2sub_v1.2_gwl_mon.csv`: monthly GWL per well (one column per well). Multiple records in a month are averaged.
- `cr2sub_v1.2_gwl_mon_clean.csv`: same series after outlier removal.
- `cr2sub_v1.2_attributes.csv`: one row per well (`cr2sub_id` = DGA well code without the last check digit):

| Columns | Content |
|---|---|
| `dga_*` | Metadata as reported by DGA (code, name, basin, sub-basin, coordinates, elevation) |
| `cr2sub_well_lat/lon`, `cr2sub_well_utm_*_h19` | Well coordinates (degrees; UTM 19S, m); the original DGA coordinates remain in `dga_well_*` |
| `cr2sub_well_elev`, `cr2sub_well_slp` | Elevation (m a.s.l.) and slope (degrees) at the well, from FABDEM resampled to 300 m |
| `cr2sub_well_gwl_mon_{stat}_{start}_{end}`, `cr2sub_well_gwl_mon_clean_{stat}_{start}_{end}` | Mean, standard deviation (`sd`) and coefficient of variation (`cv`) of monthly GWL observations for the configured period, raw and after outlier removal; `NA` when a well has no observations in that period |
| `cr2sub_well_in_basin_camels`, `cr2sub_camels_*` | Smallest CAMELS-CL v2 basin containing the well, plus basin-mean elevation and slope, period-mean annual precipitation (`pr`), mean temperature (`tmean`) and reference evapotranspiration (`et0`), and period-mean annual aridity and snow fraction |
| `cr2sub_well_in_basin_bna`, `cr2sub_bna_*` | Same basin and climate attributes for the smallest BNA basin |
| `cr2sub_clsoilmap_*` | Available water content (mm), saturated conductivity (cm/day), bulk density (g/cm³), clay and sand (%) for 0-100 and 100-200 cm, from CLSoilMaps (500 m buffer around the well). NA outside CLSoilMaps coverage |

The period-based columns use the configured start and end years in their names (for example, `cr2sub_camels_pr_yr_mean_1980_2010`, `cr2sub_camels_aridity_1980_2010`, and `cr2sub_bna_sf_1980_2010`). Basin aridity and snow fraction are calculated annually as the ratio of basin-mean annual values (`et0/pr` and `sf/pr`), then averaged over the configured period.

## Input (`input/`)

| Folder | Content |
|---|---|
| `DGA_GWL_obs_descarga/` | DGA export `AM006T0034616_1.xlsx` (raw observations and well metadata), obtained through the Transparency Law |
| `camels_cl_basins/` | CAMELS-CL v2 basin boundaries (`.gpkg`, 860 basins), from [CAMELS-CL](https://www.camels.cr2.cl/) |
| `BNA_basins/` | BNA basin boundaries, from the [DGA repository](https://dga.mop.gob.cl/mapoteca-digital/) |
| `cr2sub_basin_lookup_tables/` | Well-to-basin tables and the script that builds them |
| `dem/` | FABDEM elevation and slope, 300 m ([source](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn)) |
| `cr2met/` | CR2MET download and daily-processing scripts, plus the annual NetCDF rasters used by the pipeline; daily source files stay under `CR2MET_DATA_DIR` ([source](https://ftp.cr2.cl/browse/cr2met/v2.5)) |
| `CLSoilMaps/` | Version-specific soil attributes CSV, named with the CR2SUB tag and version; generated for all wells when missing, then updated for new wells |
| `other_data/` | Simplified national boundary (BCN-SIIT), used to validate well coordinates |

The CLSoilMaps rasters are large (~16.5 GB) and are not in git. They are needed when the version-specific CSV is first generated or when new wells are added: copy `.env.example` to `.env`, set `CR2SUB_HEAVY_DATA` to a local folder, then run `cd input/CLSoilMaps && Rscript download_clsoilmaps.R` from the repository root ([source](https://zenodo.org/records/7464210)).

CR2MET daily NetCDF files are kept outside the repository in the folder configured by `CR2MET_DATA_DIR` in `.env`. Install the Python requirements with `python3 -m pip install -r scripts/functions/python_requirements.txt`, then, from the repository root, download the daily files with `python3 input/cr2met/download_cr2met_ftp.py` and generate annual rasters with `python3 input/cr2met/process_cr2met_daily.py`. The processor reads daily inputs from `CR2MET_DATA_DIR` and writes only annual files to `input/cr2met/`; these are consumed by the main pipeline and can be versioned in Git.

## Running the pipeline

Configure the tag, version, statistical period, input/output paths, CLSoilMaps buffer and diagnostic plot settings in `scripts/main_pipeline.py`. For a new CR2MET release, update its paths there. `.env` defines `CR2SUB_HEAVY_DATA` for the large CLSoilMaps rasters and `CR2MET_DATA_DIR` for CR2MET daily downloads. The processor writes annual aggregates to `input/cr2met/`, where the pipeline reads them.

Requires R (packages `readxl`, `zoo`, `terra`, `sf`, `dplyr`, `readr`, `tibble`) and Python with the packages in `scripts/functions/python_requirements.txt`. Configure `.env` as described above. From the repository root:

```bash
python3 scripts/main_pipeline.py
```

| Step | Script | Does |
|---|---|---|
| 01 | `01_consolidate.R` | Reads the DGA export, builds well metadata and monthly series |
| 01b | `01b_extract_clsoilmaps.R` | Creates the version-specific CLSoilMaps CSV on the first run, then adds attributes for wells not yet included |
| 02 | `02_remove_outliers.ipynb` | Outlier removal |
| 03 | `03_process_attributes.R` | Builds `cr2sub_v1.2_attributes.csv` |

The CLSoilMaps CSV filename is configured in `scripts/main_pipeline.py` and includes the tag and CR2SUB version (for example, `cr2sub_CLSoilMaps_data_cr2sub_v1.2.csv`) so each version remains separate and can be compared with previous ones. Step 01b creates it for all wells if missing, or preserves existing rows and adds new wells. This requires the rasters downloaded with `cd input/CLSoilMaps && Rscript download_clsoilmaps.R`.

If wells are added or removed, rebuild the basin tables before step 03 with `Rscript input/cr2sub_basin_lookup_tables/lookup_table_script.R` (after step 01).
