#!/usr/bin/env bash
# Setup for module M03. Works in GitHub Codespaces, in a plain Linux shell and
# on macOS. Run it once, from inside modules/M03_immunoinformatics.
#
#     bash setup.sh
#
# It does three things and nothing else:
#   1. installs the Python packages this module needs
#   2. downloads the MHCflurry weights, about 200 MB
#   3. fetches the public datasets from UniProt, DeepImmuno and NetTCR
#
# It does NOT install NetMHCpan, NetMHCIIpan, SignalP or DeepTMHMM. Those need
# an academic licence and a manual download, and they cannot be installed from
# a script. The exercises that use them read from data/precomputed instead, so
# nothing is blocked.

set -euo pipefail

cd "$(dirname "$0")"

echo "=============================================================="
echo " M03 Immunoinformatics, setup"
echo "=============================================================="

PY="${PYTHON:-python3}"
echo
echo "[1/3] Python packages"
"$PY" -m pip install --quiet --upgrade pip
"$PY" -m pip install --quiet -r requirements.txt
echo "      done"

echo
echo "[2/3] MHCflurry weights, this is the slow step"
if "$PY" -c "from mhcflurry import Class1PresentationPredictor as P; P.load()" 2>/dev/null; then
  echo "      already present"
else
  "$PY" -m mhcflurry.downloads fetch models_class1_presentation
  echo "      done"
fi

echo
echo "[3/3] Public datasets"
"$PY" utils/download_data.py

echo
echo "=============================================================="
echo " Ready. Start here:"
echo
echo "   python3 solutions/ex00_fhaC_demo.py"
echo
echo " That runs the whole workflow on fhaC and trains nothing."
echo " See RUNBOOK.md for what every step means."
echo "=============================================================="
