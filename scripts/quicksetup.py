#!/usr/bin/env python3
"""Quick setup helper for the CR2SUB processing workflow.

This module offers two entry points:
  * :func:`run_quicksetup` – prepare the minimal runtime environment by writing
    the `.env` file and installing Python/R dependencies.
  * :func:`run_full_setup` – run the quick setup, ensure the runtime environment
    variables are available, verify required executables, and execute the
    pipeline steps that are passed in.

Run the standalone quick setup from the repository root:
    python scripts/quicksetup.py
"""

from __future__ import annotations

import os
import os.path
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import importlib
import shutil

REPO_ROOT = Path(__file__).resolve().parent.parent
ENV_PATH = REPO_ROOT / ".env"

PYTHON_PACKAGES = [
    "numpy",
    "pandas",
    "matplotlib",
    "seaborn",
    "jupyter",
    "nbconvert",
]

R_PACKAGES = [
    "readxl",
    "zoo",
    "terra",
]

ENV_DEFAULTS = {
    "CR2SUB_TAG": "cr2sub",
    "CR2SUB_VERSION": "v1",
    "CR2SUB_INPUT_DIR": "input",
    "CR2SUB_OUTPUT_DIR": "cr2sub",
    "CR2SUB_TMP_DIR": "tmp",
}


@dataclass
class CommandResult:
    description: str
    returncode: int


def load_env_file(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not path.exists():
        return values

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def write_env_file(path: Path, data: Dict[str, str]) -> None:
    header = [
        "# Environment configuration for the CR2SUB processing workflow",
        "# Update any value that differs from your local setup.",
    ]
    body = [f"{key}={data[key]}" for key in sorted(data)]
    content = "\n".join(header + body) + "\n"
    path.write_text(content)


def prepare_env_file() -> Tuple[Dict[str, str], list[str]]:
    env_values = load_env_file(ENV_PATH)
    new_keys: list[str] = []
    for key, value in ENV_DEFAULTS.items():
        if key not in env_values:
            env_values[key] = value
            new_keys.append(key)
    write_env_file(ENV_PATH, env_values)
    return env_values, new_keys


def run_command(cmd: list[str], description: str) -> CommandResult:
    print(f"[setup] {description}")
    completed = subprocess.run(cmd, cwd=REPO_ROOT)
    return CommandResult(description=description, returncode=completed.returncode)


def install_python_packages() -> CommandResult:
    pip_cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "pip"]
    upgrade_res = run_command(pip_cmd, "Upgrading pip")
    if upgrade_res.returncode != 0:
        return upgrade_res

    cmd = [sys.executable, "-m", "pip", "install", *PYTHON_PACKAGES]
    return run_command(cmd, "Installing Python packages")


def install_r_packages() -> CommandResult | None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        print("[warn] Rscript not found. Skipping R package installation.")
        print("       Install R and rerun this setup to install packages automatically.")
        return None

    pkg_vector = ", ".join(f'"{pkg}"' for pkg in R_PACKAGES)
    expr = (
        "pkgs <- c(" + pkg_vector + "); "
        "to_install <- setdiff(pkgs, rownames(installed.packages())); "
        "if (length(to_install)) install.packages(to_install, repos=\"https://cloud.r-project.org\")"
    )
    cmd = [rscript, "-e", expr]
    return run_command(cmd, "Installing R packages")


def ensure_env_loaded(env_values: Dict[str, str]) -> None:
    empty = [key for key, value in env_values.items() if not value]
    if empty:
        raise ValueError(
            "Found empty environment defaults for: " + ", ".join(empty)
        )

    applied = []
    for key, value in env_values.items():
        if os.environ.get(key) is None:
            os.environ[key] = value
            applied.append(key)

    if applied:
        print(f"[info] Exported defaults into session for: {', '.join(applied)}")
    else:
        print("[info] Runtime environment already defines required variables")


def verify_python_environment() -> None:
    conda_prefix = os.environ.get("CONDA_PREFIX")
    venv_path = os.environ.get("VIRTUAL_ENV")
    using_virtualenv = sys.prefix != sys.base_prefix or venv_path is not None
    if conda_prefix:
        name = os.environ.get("CONDA_DEFAULT_ENV") or os.path.basename(conda_prefix)
        print(f"[info] Active Conda environment: {name}")
    elif using_virtualenv:
        name = os.path.basename(venv_path or sys.prefix)
        print(f"[info] Active virtualenv: {name}")
    else:
        print("[warn] No dedicated Python environment detected; using system interpreter")


def verify_binaries() -> None:
    required = {"Rscript": "Rscript"}
    missing = [name for name, binary in required.items() if shutil.which(binary) is None]
    if missing:
        raise FileNotFoundError(
            "Missing required executables: " + ", ".join(missing)
        )

    try:
        importlib.import_module("nbconvert")
    except ModuleNotFoundError as exc:
        raise FileNotFoundError(
            "Missing Python package 'nbconvert'. Run the quick setup to install it."
        ) from exc


def run_step(description: str, cmd: List[str]) -> int:
    print(f"[pipeline] {description}")
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode == 0:
        print("[ok] Step completed")
    else:
        print(f"[error] Step failed with code {result.returncode}")
    return result.returncode


def run_pipeline(steps: List[Tuple[str, List[str]]]) -> int:
    for description, cmd in steps:
        code = run_step(description, cmd)
        if code != 0:
            return code
    return 0


def run_quicksetup() -> int:
    print("Preparing CR2SUB environment...")

    env_values, new_keys = prepare_env_file()

    if new_keys:
        print(f"[ok] .env created/updated with defaults for: {', '.join(new_keys)}")
    else:
        print("[ok] .env already contained all required keys")

    py_res = install_python_packages()
    if py_res.returncode != 0:
        print("[error] Python dependency installation failed")
        return py_res.returncode
    print("[ok] Python dependencies installed")

    r_res = install_r_packages()
    if r_res is not None and r_res.returncode != 0:
        print("[error] R dependency installation failed")
        return r_res.returncode
    if r_res is None:
        print("[warn] Skipped R packages (Rscript unavailable)")
    else:
        print("[ok] R dependencies installed")

    print("Environment bootstrap complete.")
    return 0


def run_full_setup(steps: List[Tuple[str, List[str]]]) -> int:
    setup_code = run_quicksetup()
    if setup_code != 0:
        print("[error] Quick setup failed; aborting pipeline execution.")
        return setup_code

    env_values, _ = prepare_env_file()
    try:
        ensure_env_loaded(env_values)
        verify_python_environment()
        verify_binaries()
    except (ValueError, FileNotFoundError) as exc:
        print(f"[error] {exc}")
        return 1

    return run_pipeline(steps)


if __name__ == "__main__":
    sys.exit(run_quicksetup())
