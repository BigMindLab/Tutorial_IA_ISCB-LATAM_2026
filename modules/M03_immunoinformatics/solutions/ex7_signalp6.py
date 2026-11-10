"""Exercise 7 - SignalP 6.0: does fhaC enter the secretory pathway? NOT TESTED HERE.

Needs a local install (DTU academic license, manual download) providing `signalp6`.
Organism: SignalP 6.0 only distinguishes "eukarya" from "other" (Gram-negative, Gram-positive,
archaea all use "other"). fhaC is from a Gram-negative bacterium, so the default is "other".

Usage:
    python ex7_signalp6.py fhaC.fasta [other|eukarya]
Output: results/signalp6_summary.csv with the predicted type, per-class probabilities and
cleavage-site position, parsed from SignalP's prediction_results.txt. Per-residue region
probabilities (n-, h-, c-region, slide 1 panel d) are in the file(s) SignalP writes to
results/signalp6/ (use --format all to also get plots).

Assumes the SignalP 6.0 prediction_results.txt layout: lines starting with '#' are comments;
data lines are tab-separated: ID, Prediction, OTHER, SP(Sec/SPI), LIPO(Sec/SPII), TAT(Tat/SPI),
TATLIPO(Tat/SPII), PILIN(Sec/SPIII), CS Position. Check against your own output.
Reminder from the slides: a signal peptide says where the protein goes, not whether it is antigenic.
"""
import csv
import os
import subprocess
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c


def main(fasta, organism="other"):
    out_dir = os.path.join(c.RESULTS, "signalp6")
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run(["signalp6", "--fastafile", fasta, "--organism", organism,
                    "--output_dir", out_dir, "--format", "txt", "--mode", "fast"], check=True)
    fields = ["id", "prediction", "OTHER", "SP_Sec_SPI", "LIPO_Sec_SPII", "TAT_Tat_SPI",
              "TATLIPO_Tat_SPII", "PILIN_Sec_SPIII", "cs_position"]
    rows = []
    with open(os.path.join(out_dir, "prediction_results.txt")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            parts += [""] * (len(fields) - len(parts))
            rows.append(dict(zip(fields, parts[:len(fields)])))
    c.write_csv(os.path.join(c.RESULTS, "signalp6_summary.csv"), rows, fields)
    for r in rows:
        print(r["id"], r["prediction"], r["cs_position"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta",
         sys.argv[2] if len(sys.argv) > 2 else "other")
