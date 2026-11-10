#!/usr/bin/env bash
# Run with bash setup.sh on Linux, WSL or Codespaces.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

usage() {
  cat <<'HELP'
Usage: bash setup.sh [options]

Creates .venv, installs the CPU dependencies, fetches models and public data,
and checks imports plus a real prediction. No licensed software is installed.

  --env-dir PATH       Create/reuse a virtual environment at PATH (default .venv)
  --use-current-env    Install into the active venv/Conda environment instead
  --skip-models        Install packages without fetching/checking MHCflurry weights
  --skip-data          Do not fetch/check public datasets
  -h, --help           Show this help

PYTHON=/path/to/python3.11 bash setup.sh selects the base interpreter.
Supported Python: 3.10, 3.11 and 3.12. Conda alternative: see QUICKSTART.md.
HELP
}

ENV_DIR=".venv"
CURRENT_ENV=0
CUSTOM_ENV_DIR=0
SKIP_MODELS=0
SKIP_DATA=0
while (($#)); do
  case "$1" in
    --env-dir)
      if (($# < 2)) || [[ -z "$2" || "$2" == --* ]]; then
        echo "--env-dir needs a path." >&2; exit 2
      fi
      ENV_DIR="$2"; CUSTOM_ENV_DIR=1; shift 2 ;;
    --use-current-env) CURRENT_ENV=1; shift ;;
    --skip-models) SKIP_MODELS=1; shift ;;
    --skip-data) SKIP_DATA=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done
if ((CURRENT_ENV && CUSTOM_ENV_DIR)); then
  echo "Choose either --env-dir or --use-current-env." >&2; exit 2
fi

if ((CURRENT_ENV)); then
  PY="${PYTHON:-python}"
else
  PY="${PYTHON:-python3}"
fi
# Reuse the environment interpreter even if the base Python is no longer on PATH.
if ((!CURRENT_ENV)) && [[ -f "$ENV_DIR/pyvenv.cfg" ]]; then
  PY="$ENV_DIR/bin/python"
  [[ -x "$PY" ]] || { echo "Broken virtual environment: $ENV_DIR" >&2; exit 1; }
fi
command -v "$PY" >/dev/null || {
  echo "Python not found: $PY. Install Python 3.10-3.12 or use Conda (QUICKSTART.md)." >&2
  exit 1
}
"$PY" - <<'PY'
import sys
if not (3, 10) <= sys.version_info[:2] <= (3, 12):
    sys.exit("Use Python 3.10-3.12 (recommended: 3.11). See QUICKSTART.md for Conda.")
print("Python:", sys.version.split()[0], "at", sys.executable)
PY

if ((CURRENT_ENV)); then
  "$PY" - <<'PY'
import os, sys
from pathlib import Path
conda = os.environ.get("CONDA_PREFIX")
is_conda = conda and Path(conda).resolve() == Path(sys.prefix).resolve()
if sys.prefix == sys.base_prefix and not is_conda:
    sys.exit("Activate a venv/Conda environment before using --use-current-env.")
PY
else
  if [[ ! -f "$ENV_DIR/pyvenv.cfg" ]]; then
    if [[ -e "$ENV_DIR" ]]; then
      echo "Refusing to overwrite $ENV_DIR: it is not a Python virtual environment." >&2
      exit 1
    fi
    if ! "$PY" -m venv "$ENV_DIR"; then
      echo "Could not create $ENV_DIR. On Ubuntu install python3-venv (matching your Python)," >&2
      echo "or use the Conda instructions in QUICKSTART.md." >&2
      exit 1
    fi
  fi
  PY="$ENV_DIR/bin/python"
fi

if ! "$PY" -m pip --version >/dev/null 2>&1; then
  echo "The environment has no pip. Install the matching python3-venv package or use Conda." >&2
  echo "See QUICKSTART.md and select a fresh --env-dir for the next attempt." >&2
  exit 1
fi

echo "[1/4] Installing packages into $PY"
# Install only into the selected environment.
"$PY" -m pip install -r requirements-cpu.txt
"$PY" -m pip check

if ((SKIP_MODELS)); then
  echo "[2/4] MHCflurry weights skipped"
else
  echo "[2/4] MHCflurry presentation weights (reuses existing downloads)"
  # mhcflurry.downloads has no CLI entry point. Invoke the real downloads CLI
  # with this environment's Python, regardless of the parent shell's PATH.
  "$PY" -c 'from mhcflurry.downloads_command import run; run()' fetch models_class1_presentation
fi

if ((SKIP_DATA)); then
  echo "[3/4] Public datasets skipped"
else
  echo "[3/4] Public datasets"
  "$PY" utils/download_data.py
fi

echo "[4/4] Checking the environment"
CHECK_ARGS=(utils/check_environment.py)
if ((SKIP_MODELS)); then CHECK_ARGS+=(--packages-only); fi
if ((!SKIP_DATA)); then CHECK_ARGS+=(--check-data); fi
"$PY" "${CHECK_ARGS[@]}"

echo
echo "Setup complete."
if ((!CURRENT_ENV)); then
  printf 'Activate in your current terminal: source %q/bin/activate\n' "$ENV_DIR"
fi
if ((SKIP_MODELS || SKIP_DATA)); then
  echo "Models or datasets were skipped. Rerun setup without those flags before the demo."
else
  echo "Then run: python solutions/ex00_fhaC_demo.py"
fi
