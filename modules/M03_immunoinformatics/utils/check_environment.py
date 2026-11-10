#!/usr/bin/env python3
"""Check installed packages and run small CPU predictions without rewriting results."""
from __future__ import annotations

import argparse
import importlib
import importlib.metadata
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packages-only", action="store_true",
                        help="Skip pretrained model checks (no weights needed)")
    parser.add_argument("--check-data", action="store_true",
                        help="Also check the downloaded demo/training inputs")
    args = parser.parse_args()
    if not (3, 10) <= sys.version_info[:2] <= (3, 12):
        parser.error("Use Python 3.10-3.12; see QUICKSTART.md")
    print("Interpreter:", sys.executable)
    for package, module in [
        ("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"),
        ("scikit-learn", "sklearn"), ("matplotlib", "matplotlib"),
        ("biopython", "Bio"), ("torch", "torch"),
        ("mhcflurry", "mhcflurry"), ("mhcgnomes", "mhcgnomes"),
    ]:
        importlib.import_module(module)
        print("[ok]", package, importlib.metadata.version(package))

    import numpy as np
    import torch
    from encoding import load_blosum50, load_paratopes, load_pca12

    assert load_blosum50().shape == (20, 20), "Invalid BLOSUM50 matrix"
    if args.check_data:
        from download_data import FILES
        for name in ["fhaC.fasta", *FILES]:
            path = ROOT / "data" / name
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing {path}. Run python utils/download_data.py")
        assert load_pca12().shape == (21, 12), "Invalid PCA matrix"
        assert not load_paratopes().empty, "Empty HLA paratope table"
        print("[ok] Public inputs and encoding tables")
    if args.packages_only:
        print("Package checks passed. Pretrained models were not checked.")
        return 0

    sys.path.insert(0, str(ROOT / "solutions"))
    from ex12_immunogenicity import DeepImmunoCNN
    from ex14_prioritization import HLA_PANEL
    from mhcflurry import Class1PresentationPredictor

    model = DeepImmunoCNN()
    weights = ROOT / "data" / "precomputed" / "ex12_immunogenicity_cnn.pt"
    model.load_state_dict(torch.load(weights, map_location="cpu", weights_only=True))
    model.eval()
    with torch.no_grad():
        score = model(torch.zeros(1, 10, 12), torch.zeros(1, 46, 12))
    assert torch.isfinite(score).all(), "Non-finite DeepImmuno prediction"
    print("[ok] Tutorial CNN checkpoint and CPU inference")

    predictor = Class1PresentationPredictor.load()
    missing = set(HLA_PANEL) - set(predictor.supported_alleles)
    assert not missing, f"Missing tutorial HLA alleles: {sorted(missing)}"
    result = predictor.predict(peptides=["RLFDPAPLV"], alleles=["HLA-A*02:01"],
                               include_affinity_percentile=True, verbose=0)
    columns = ["affinity", "affinity_percentile", "processing_score",
               "presentation_score", "presentation_percentile"]
    assert len(result) == 1 and np.isfinite(result[columns].to_numpy()).all(), result
    print("[ok] MHCflurry presentation prediction:")
    print(result[["peptide", *columns]].to_string(index=False))
    print("Environment checks passed.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        print("Run bash setup.sh (or --use-current-env for Conda). See QUICKSTART.md.",
              file=sys.stderr)
        sys.exit(1)
