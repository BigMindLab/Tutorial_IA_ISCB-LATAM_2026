#!/usr/bin/env python3
"""Download the inputs for module M03.1.

This repository does NOT redistribute third party data. Everything is fetched
from its original source so that each project's licence is respected.

Sources
-------
fhaC.fasta                     The sequence sent to participants. Identical to
                               UniProtKB P35077, Bordetella pertussis Tohama I,
                               and shipped in the repository already.
P35077.json                    UniProtKB entry with its curated annotation
                               (signal peptide, transmembrane strands, domains).
                               Exercises 11 and 14 use it as a topology filter.
after_pca.txt                  DeepImmuno (Li et al., Brief Bioinform 2021).
                               21 x 12 matrix: 12 principal components of 553
                               AAindex properties, one row per amino acid plus X.
hla2paratopeTable_aligned.txt  DeepImmuno. 46 residue paratope pseudo sequence
                               for each HLA allele.
remove0123.csv                 DeepImmuno. 9,055 IEDB peptides with an
                               immunogenicity label. Slide 22 cites 9,056.
train_ab_90_alphabeta.csv      NetTCR-2.0 (Montemurro et al., Commun Biol 2021).
                               Paired CDR3 alpha and beta plus peptide, with
                               partitions for cross validation.

Usage:  python utils/download_data.py
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

DEEPIMMUNO = "https://raw.githubusercontent.com/frankligy/DeepImmuno/main/data/"
NETTCR = "https://raw.githubusercontent.com/mnielLab/NetTCR-2.0/main/data/"

FILES = {
    "P35077.json": "https://rest.uniprot.org/uniprotkb/P35077.json",
    "after_pca.txt": DEEPIMMUNO + "after_pca.txt",
    "hla2paratopeTable_aligned.txt": DEEPIMMUNO + "hla2paratopeTable_aligned.txt",
    "remove0123.csv": DEEPIMMUNO + "remove_low_negative/remove0123.csv",
    "train_ab_90_alphabeta.csv": NETTCR + "train_ab_90_alphabeta.csv",
}


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    failures = 0
    for name, url in FILES.items():
        target = DATA_DIR / name
        if target.exists() and target.stat().st_size > 0:
            print("[present ] %s (%d bytes)" % (name, target.stat().st_size))
            continue
        print("[fetching] %s" % name)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "iscb-tutorial"})
            with urllib.request.urlopen(request, timeout=180) as response:
                target.write_bytes(response.read())
            print("[ok      ] %s (%d bytes)" % (name, target.stat().st_size))
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print("[ERROR   ] %s: %s" % (name, exc))

    if not (DATA_DIR / "fhaC.fasta").exists():
        print("[ERROR   ] fhaC.fasta is missing from data/. It ships with the repository.")
        failures += 1

    print("\nData directory: %s" % DATA_DIR)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
