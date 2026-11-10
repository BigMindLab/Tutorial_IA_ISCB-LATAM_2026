"""Exercise 10 - NetTCR (TCR-peptide recognition): input preparation only.

fhaC alone gives no TCR. NetTCR needs CDR3 alpha and CDR3 beta sequences, which you must supply.
Put them in tcrs.csv with columns: name,CDR3a,CDR3b (no placeholders are invented here).

What this script does: takes the PRIME top candidates from results/prime_class_I.csv (run ex6
first), crosses them with your TCRs, and writes results/nettcr_input.csv with columns
peptide,CDR3a,CDR3b. Check these columns against the input requirements of the
NetTCR release you use before passing the table to that package.

Usage:
    python ex10_nettcr_prepare.py tcrs.csv [top_n]

Slide caveat, to discuss: the NetTCR authors state the model can at present only be applied to the
peptides in its training data. A score for an fhaC peptide, an unseen TCR or an unseen HLA context
is extrapolation. NetTCR-2.0 (slide) is not NetTCR-2.2 (2024), which has a different recipe.
"""
import csv
import os
import sys

import pandas as pd

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c


def main(tcr_file, top_n=20):
    prime = pd.read_csv(os.path.join(c.RESULTS, "prime_class_I.csv"))
    top = prime[prime.length == 9].sort_values("prime_rank").drop_duplicates("peptide").head(top_n)
    tcrs = list(csv.DictReader(open(tcr_file)))
    out = os.path.join(c.RESULTS, "nettcr_input.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["peptide", "CDR3a", "CDR3b"])
        for p in top.peptide:
            for t in tcrs:
                w.writerow([p, t["CDR3a"], t["CDR3b"]])
    print("wrote", out, "(%d peptides x %d TCRs)" % (len(top), len(tcrs)))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python ex10_nettcr_prepare.py tcrs.csv [top_n]")
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20)
