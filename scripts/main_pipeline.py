#!/usr/bin/env python3
"""Run the CR2SUB processing pipeline (01 → 01b → 02 → 03)."""
import json
import os
import subprocess
import sys
from pathlib import Path


PIPELINE_CONFIG = {
    "TAG": "cr2sub",
    "VERSION": "v1.2",
    "INI_MEAN_VALUES": "1980",
    "END_MEAN_VALUES": "2010",
    "CLSOIL_BUFFER_M": "500",
    "SAVE_DIAGNOSTIC_PLOTS": "0",

    "PATHS": {
        "TMP_DIR": ".tmp",
        "OUTPUT_DIR": "{tag}",
        "DGA_INPUT_FILE": "input/DGA_GWL_obs_descarga/AM006T0034616_1.xlsx",
        "CHILE_BOUNDARIES_FILE": "input/other_data/Chile_boundaries_simplified_from_SIIT.gpkg",
        "METADATA_FILE": "{tmp_dir}/{tag}_{version}_metadata.csv",
        "GWL_TS_FILE": "{output_dir}/{tag}_{version}_gwl_mon.csv",
        "GWL_TS_CLEAN_FILE": "{output_dir}/{tag}_{version}_gwl_mon_clean.csv",
        "ATTRIBUTES_FILE": "{output_dir}/{tag}_{version}_attributes.csv",
        "CLSOIL_DATA_FILE": "input/CLSoilMaps/cr2sub_CLSoilMaps_data_{tag}_{version}.csv",
        "WELLS_BUFFER_FILE": "{tmp_dir}/wells_buffer_500m.geojson",
        "CLSOIL_NEW_WELLS_FILE": "{tmp_dir}/clsoilmaps_new_wells.csv",
        "CR2MET_PR_FILE": (
            "input/cr2met/CR2MET_pr_v2.5_v2.5_R1_yr_1960_2025.nc"
        ),
        "CR2MET_ET0_FILE": (
            "input/cr2met/CR2MET_et0_v2.5_v2.5_R1_yr_1960_2025.nc"
        ),
        "CR2MET_SF_FILE": (
            "input/cr2met/CR2MET_sf_v2.5_v2.5_R1_yr_1960_2025.nc"
        ),
        "CR2MET_TMEAN_FILE": (
            "input/cr2met/CR2MET_tmean_v2.5_v2.5_R1_yr_1960_2025.nc"
        ),
        "DEM_FILE": "input/dem/dem_masl_fabdemv1.2_2015_300m_epsg4326.tif",
        "SLOPE_FILE": "input/dem/slope_deg_fabdemv1.2_2015_300m_epsg4326.tif",
        "CAMELS_BASINS_FILE": "input/camels_cl_basins/camels_cl_v2_boundaries.gpkg",
        "BNA_BASINS_FILE": "input/BNA_basins/Cuencas_BNA.shp",
        "CAMELS_JOIN_FILE": (
            "input/cr2sub_basin_lookup_tables/"
            "cr2sub_{version}_join_table_with_camels_cl_v2_basins.csv"),
        "BNA_JOIN_FILE": (
            "input/cr2sub_basin_lookup_tables/"
            "cr2sub_{version}_join_table_with_bna_basins.csv"),
        "PLOTS_DIR": "figures/outlier_removal",
    },
}

REQUIRED_INPUT_PATHS = (
    "DGA_INPUT_FILE",
    "CHILE_BOUNDARIES_FILE",
    "DEM_FILE",
    "SLOPE_FILE",
    "CR2MET_PR_FILE",
    "CR2MET_ET0_FILE",
    "CR2MET_SF_FILE",
    "CR2MET_TMEAN_FILE",
    "CAMELS_BASINS_FILE",
    "BNA_BASINS_FILE",
    "CAMELS_JOIN_FILE",
    "BNA_JOIN_FILE",
)

OUTPUT_PATHS = (
    "METADATA_FILE",
    "GWL_TS_FILE",
    "GWL_TS_CLEAN_FILE",
    "ATTRIBUTES_FILE",
    "WELLS_BUFFER_FILE",
    "CLSOIL_NEW_WELLS_FILE",
)


def load_local_environment(project_root: Path) -> None:
    """Load .env values without overriding variables already in the environment."""
    env_file = project_root / ".env"
    if not env_file.is_file():
        return
    for line in env_file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip("\"'"))


def build_pipeline_environment(project_root: Path) -> dict[str, str]:
    """Resolve the user-editable pipeline settings into child-process variables."""
    load_local_environment(project_root)
    tag = PIPELINE_CONFIG["TAG"]
    version = PIPELINE_CONFIG["VERSION"]
    paths = PIPELINE_CONFIG["PATHS"]
    path_context = {
        "tag": tag,
        "version": version,
        "tmp_dir": paths["TMP_DIR"],
        "output_dir": paths["OUTPUT_DIR"].format(tag=tag, version=version),
    }

    config = {
        key: str(value)
        for key, value in PIPELINE_CONFIG.items()
        if key != "PATHS"
    }
    start_year = int(config["INI_MEAN_VALUES"])
    end_year = int(config["END_MEAN_VALUES"])
    if start_year > end_year:
        raise ValueError("INI_MEAN_VALUES must not exceed END_MEAN_VALUES.")

    config.update({
        key: str((project_root / value.format(**path_context)).resolve())
        for key, value in paths.items()
    })

    missing = [
        config[key]
        for key in REQUIRED_INPUT_PATHS
        if not Path(config[key]).is_file()
    ]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Pipeline input file(s) not found:\n{formatted}")

    return {
        f"CR2SUB_{key}": value
        for key, value in config.items()
    }


def run_r_script(script_path: Path, workdir: Path, env: dict[str, str]) -> None:
    """Execute an R script with the given working directory."""
    print(f"[cr2sub] Running {script_path.name}...")
    subprocess.run([
        "Rscript",
        str(script_path)
    ], check=True, cwd=str(workdir), env=env)
    print(f"[cr2sub] Finished {script_path.name}.\n")


def execute_notebook(
    notebook_path: Path, workdir: Path, env: dict[str, str]
) -> None:
    """Execute a lightweight .ipynb notebook by running its code cells sequentially."""
    print(f"[cr2sub] Running {notebook_path.name}...")
    with notebook_path.open("r", encoding="utf-8") as nb_file:
        notebook = json.load(nb_file)

    current_cell = None
    runtime_env = {
        "__name__": "__main__",
        "__file__": str(notebook_path),
        **env,
    }

    original_cwd = Path.cwd()
    original_sys_path = list(sys.path)
    try:
        os.chdir(workdir)
        if str(workdir) not in sys.path:
            sys.path.insert(0, str(workdir))

        for idx, cell in enumerate(notebook.get("cells", [])):
            if cell.get("cell_type") != "code":
                continue

            source = "".join(cell.get("source", []))
            if not source.strip():
                continue

            current_cell = idx
            print(f"[cr2sub] Executing notebook cell {idx}...")
            code_obj = compile(source, f"{notebook_path.name}#cell{idx}", "exec")
            exec(code_obj, runtime_env)
    except Exception as exc:
        cell_msg = f"cell {current_cell}" if current_cell is not None else "initialisation"
        raise RuntimeError(
            f"Execution failed while running {notebook_path.name} at {cell_msg}."
        ) from exc
    finally:
        os.chdir(original_cwd)
        sys.path[:] = original_sys_path

    print(f"[cr2sub] Finished {notebook_path.name}.\n")


def main() -> None:
    script_path = Path(__file__).resolve()
    scripts_dir = script_path.parent
    project_root = scripts_dir.parent
    pipeline_env = build_pipeline_environment(project_root)
    for key in OUTPUT_PATHS:
        Path(pipeline_env[f"CR2SUB_{key}"]).parent.mkdir(
            parents=True, exist_ok=True
        )

    env = os.environ.copy()
    env.update(pipeline_env)

    run_r_script(project_root / "scripts/01_consolidate.R", project_root, env)
    run_r_script(project_root / "scripts/01b_extract_clsoilmaps.R", project_root, env)
    soil_file = Path(pipeline_env["CR2SUB_CLSOIL_DATA_FILE"])
    if not soil_file.is_file():
        raise FileNotFoundError(
            f"CLSoilMaps data file was not generated: {soil_file}"
        )
    execute_notebook(scripts_dir / "02_remove_outliers.ipynb", scripts_dir, env)
    run_r_script(project_root / "scripts/03_process_attributes.R", project_root, env)


if __name__ == "__main__":
    main()
