#!/usr/bin/env python3
"""Run the CR2SUB processing pipeline (01 → 02 → 03)."""
import json
import os
import subprocess
import sys
from pathlib import Path


def run_r_script(script_path: Path, workdir: Path) -> None:
    """Execute an R script with the given working directory."""
    print(f"[cr2sub] Running {script_path.name}...")
    subprocess.run([
        "Rscript",
        str(script_path)
    ], check=True, cwd=str(workdir))
    print(f"[cr2sub] Finished {script_path.name}.\n")


def execute_notebook(notebook_path: Path, workdir: Path) -> None:
    """Execute a lightweight .ipynb notebook by running its code cells sequentially."""
    print(f"[cr2sub] Running {notebook_path.name}...")
    with notebook_path.open("r", encoding="utf-8") as nb_file:
        notebook = json.load(nb_file)

    current_cell = None
    runtime_env = {
        "__name__": "__main__",
        "__file__": str(notebook_path)
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

    run_r_script(project_root / "scripts/01_consolidate.R", project_root)
    execute_notebook(scripts_dir / "02_remove_outliers.ipynb", scripts_dir)
    run_r_script(project_root / "scripts/03_process_attributes.R", project_root)


if __name__ == "__main__":
    main()


