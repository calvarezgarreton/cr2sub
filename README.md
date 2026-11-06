# CR2SUB: monthly groundwater level database for Chile

CR2SUB compiles groundwater level (GWL) observations from wells of the Chilean Water Bureau (Dirección General de Aguas, DGA), downloaded from the [DGA portal](https://snia.mop.gob.cl/BNAConsultas/reportes), into a homogenized monthly database (1957-August 2026), plus metadata and topographic, hydroclimatic and soil attributes for each well.

**Developer**: Center for Climate and Resilience Research (CR2, ANID/FONDAP/1523A0002), project ANID/FONDECYT/11240924, project ANID IDeA TA TA25I10030. Responsible: Camila Alvarez-Garreton. Collaborators: Rodrigo Marinao Rivas, Juan P. Boisier, Diego Dinamarca.

**License**: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

**Citation**: Please cite the [Zenodo record](https://doi.org/10.5281/zenodo.15851099) when using this dataset. Publications describing the dataset will be added here when available.

## Versions and releases

`main` is the canonical branch for current development and contains the latest
integrated state. A release is a fixed snapshot identified by an immutable,
annotated Git tag (for example, `v1.2.0`); version-specific branches are not
used for releases. The old `cr2sub_v1.1` and `cr2sub_v1.2` branches are kept
as archived, read-only references and are not active development branches.

| Reference | Status |
|---|---|
| `main` | Canonical branch for current development; includes the integrated v1.2 release |
| `v1.2.0` | Immutable release tag for the v1.2 dataset and pipeline |
| `cr2sub_v1.2` | Archived, read-only historical branch; not for new development |
| `cr2sub_v1.1` | Archived, read-only branch for the earlier published release; not for new development |

## What's new in v1.2

- **New DGA source**: single export obtained through the Transparency Law (see [Data source](#data-source-v12)). Series now extend to August 2026 (v1.1 ended in December 2024).
- **1,133 wells** (v1.1: 1,137). Four wells of the Taltal basin are not in the new export.
- **CAMELS-CL v2** basins (860) replace CAMELS-CL v1 (516); basin assignment now uses the basin that contains the well.
- **CLSoilMaps extraction** is part of the pipeline and requires the heavy rasters (~16.5 GB).
- **Self-contained pipeline**: local paths go in `.env`, and the review of v1.1 vs v1.2 is in `review_v1.2_vs_v1.1/`.

## Data source (v1.2)

The v1.2 groundwater level data were requested from the DGA under the Chilean Transparency Law (Ley de Transparencia, Ley 20.285) through the [Portal de Transparencia](https://www.portaltransparencia.cl/PortalPdT/ingreso-sai-v2). Request N° AM006T0034616 was received on 2026-09-21 (static water levels in wells, January 1950 to August 2026, for 1,137 stations) and answered by the DGA on 2026-10-01. The response is the single export `input/DGA_GWL_obs_descarga/AM006T0034616_1.xlsx`. The DGA reply letter is in `input/DGA_GWL_obs_descarga/transparency_request_AM006T0034616_dga_reply.pdf`.

## Output (`cr2sub/`)

- `cr2sub_v1.2_gwl_mon.csv`: monthly GWL per well (one column per well). Multiple records in a month are averaged.
- `cr2sub_v1.2_gwl_mon_clean.csv`: same series after outlier removal.
- `cr2sub_v1.2_attributes.csv`: one row per well (`cr2sub_id` = DGA well code without the last check digit):

| Columns | Content |
|---|---|
| `dga_*` | Metadata as reported by DGA (code, name, basin, sub-basin, coordinates, elevation) |
| `cr2sub_well_lat`, `cr2sub_well_lon`, `cr2sub_well_utm_north_h19`, `cr2sub_well_utm_east_h19` | Well coordinates (latitude and longitude in degrees; UTM 19S northing and easting in m); the original DGA coordinates remain in `dga_well_*` |
| `cr2sub_well_elev`, `cr2sub_well_slp` | Elevation (m a.s.l.) and slope (degrees) at the well, from FABDEM resampled to 300 m |
| `cr2sub_well_gwl_mon_{stat}_1980_2010`, `cr2sub_well_gwl_mon_clean_{stat}_1980_2010` | Raw and outlier-cleaned monthly GWL statistics: mean (`mean`), standard deviation (`sd`) and coefficient of variation (`cv`). Values are `NA` if no observations are available in the period. |
| `cr2sub_well_in_basin_camels`, `cr2sub_camels_*` | Smallest CAMELS-CL v2 basin containing the well, with basin-mean elevation and slope and climate attributes (`pr`, `tmean`, `et0`, `aridity`, `sf`). |
| `cr2sub_well_in_basin_bna`, `cr2sub_bna_*` | Same attributes for the smallest BNA basin containing the well. |
| `cr2sub_clsoilmap_*` | Available water content (mm), saturated conductivity (cm/day), bulk density (g/cm³), clay and sand (%) for 0-100 and 100-200 cm, from CLSoilMaps (500 m buffer around the well). NA outside CLSoilMaps coverage |

For v1.2, all period-based statistics use 1980-2010 inclusive, and the years appear in column names (e.g. `cr2sub_well_gwl_mon_mean_1980_2010`). Basin aridity and snow fraction are annual ratios of basin-mean values (`et0/pr` and `sf/pr`), averaged over this period. Basin climate attributes use CR2MET v2.5 and must be recalculated for a future release when CR2MET v3 becomes available.

## Input (`input/`)

| Folder | Content |
|---|---|
| `DGA_GWL_obs_descarga/` | DGA export `AM006T0034616_1.xlsx` (raw observations and well metadata), obtained through the Transparency Law |
| `camels_cl_basins/` | CAMELS-CL v2 basin boundaries (`.gpkg`, 860 basins), from [CAMELS-CL](https://www.camels.cr2.cl/) |
| `BNA_basins/` | BNA basin boundaries, from the [DGA repository](https://dga.mop.gob.cl/mapoteca-digital/) |
| `cr2sub_basin_lookup_tables/` | Well-to-basin tables and the script that builds them |
| `dem/` | FABDEM elevation and slope, 300 m ([source](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn)) |
| `cr2met/` | CR2MET v2.5 download and daily-processing scripts, plus the annual NetCDF rasters used by the pipeline ([source](https://ftp.cr2.cl/browse/cr2met/v2.5)) |
| `CLSoilMaps/` | Version-specific soil attributes CSV; heavy rasters from [Zenodo](https://zenodo.org/records/7464210) are used to generate it |
| `other_data/` | Simplified national boundary (BCN-SIIT), used to validate well coordinates |

## Large data files

Large raw datasets and machine-specific paths stay outside Git and are configured locally in `.env`. Copy `.env.example` to `.env` before using the data below.

- **CLSoilMaps:** The rasters (~16.5 GB) are needed to create a version-specific CSV for the first time or add wells. Set `CR2SUB_HEAVY_DATA` to a local folder, then run `cd input/CLSoilMaps && Rscript download_clsoilmaps.R` from the repository root.
- **CR2MET:** Daily NetCDF files stay outside the repository under `CR2MET_DATA_DIR`. Install dependencies with `python3 -m pip install -r scripts/functions/python_requirements.txt`, then run `python3 input/cr2met/download_cr2met_ftp.py` and `python3 input/cr2met/process_cr2met_daily.py` from the repository root. The processor writes annual rasters to `input/cr2met/`, which are used by the pipeline and can be versioned in Git.

## Running the pipeline

Set the tag, version, statistical period, input/output paths, CLSoilMaps buffer and diagnostic plot settings in `scripts/main_pipeline.py`. Update the input paths there when adopting a new CR2MET release.

Requires R (packages `readxl`, `zoo`, `terra`, `sf`, `dplyr`, `readr`, `tibble`) and Python packages listed in `scripts/functions/python_requirements.txt`. From the repository root, run:

```bash
python3 scripts/main_pipeline.py
```

| Step | Script | Does |
|---|---|---|
| 01 | `01_consolidate.R` | Reads the DGA export, builds well metadata and monthly series |
| 01b | `01b_extract_clsoilmaps.R` | Creates the version-specific CSV (named with the CR2SUB tag and version, e.g. `cr2sub_CLSoilMaps_data_cr2sub_v1.2.csv`) for all wells on the first run; later preserves existing rows and adds new wells |
| 02 | `02_remove_outliers.ipynb` | Outlier removal |
| 03 | `03_process_attributes.R` | Builds `cr2sub_v1.2_attributes.csv` |

If wells are added or removed, rebuild the basin tables before step 03 with `Rscript input/cr2sub_basin_lookup_tables/lookup_table_script.R` (after step 01).


## Historical references

Tags preserve exact release snapshots. Older branches may be kept for comparison or reproducibility. Make future fixes from `main` and publish them with new tags; never modify existing release tags.