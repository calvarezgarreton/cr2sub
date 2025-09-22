#!/usr/bin/env bash
set -euo pipefail

# Setup local Python venv for the notebook and install deps.
# Uses existing .venv if present; otherwise creates one with python3.

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
VENV_DIR="$ROOT_DIR/.venv"

if [[ ! -d "$VENV_DIR" ]]; then
  echo "[setup] Creating venv at .venv"
  python3 -m venv "$VENV_DIR"
fi

PY="$VENV_DIR/bin/python"
PIP="$VENV_DIR/bin/pip"

echo "[setup] Ensuring pip is available"
"$PY" -m ensurepip --upgrade || true

echo "[setup] Upgrading pip/setuptools/wheel"
"$PY" -m pip install -U pip setuptools wheel

echo "[setup] Installing dependencies from requirements.txt"
"$PIP" install -r "$ROOT_DIR/scripts/functions/python_requirements.txt"

echo "[setup] Registering Jupyter kernel for this venv"
"$PY" -m ipykernel install --user --name cr2sub-py311 --display-name "Python (cr2sub)"

echo "[setup] Done. Select kernel: Python (cr2sub) in the notebook."

