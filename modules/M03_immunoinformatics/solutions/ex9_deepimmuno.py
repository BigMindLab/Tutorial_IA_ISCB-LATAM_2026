"""Exercise 9 - DeepImmuno-CNN (class I immunogenicity). NOT TESTED HERE (needs TensorFlow + repo).

Setup:
    git clone https://github.com/frankligy/DeepImmuno.git
    pip install tensorflow pandas numpy      # versions per the DeepImmuno README
Prepares the input CSV from the fhaC 9- and 10-mers for the alleles in common.py, then prints the
command to run.

VERIFIED against deepimmuno-cnn.py in the current repo:
  * --mode multiple, --intdir <input file>, --outdir <output folder> are correct.
    The name --intdir is misleading: it takes a FILE path, not a directory.
  * The allele spelling HLA-A*0201 is correct, and A*02:01, A*24:02 and B*07:02
    are all present in data/hla2paratopeTable_aligned.txt (62 alleles total).
  * The input CSV must have NO HEADER ROW. DeepImmuno reads it with
    pd.read_csv(upload, sep=',', header=None) and assigns the column names
    itself. A header line would be scored as a peptide called "peptide".
    This script therefore writes the file without a header.
  * Run it from the DeepImmuno repository root. The script loads './data/...'
    and './models/...' by relative path and fails anywhere else.
  * Output lands in <outdir>/deepimmuno-cnn-result.txt, TAB separated, with the
    score in a column named immunogenicity.

Usage:
    python ex9_deepimmuno.py fhaC.fasta
Output: results/deepimmuno_input.csv with two unnamed columns, peptide and HLA.
DeepImmuno handles 9- and 10-mers only, so 8- and 11-mers are skipped here.
Slide caveat: trained on 9,056 IEDB peptides, mostly common alleles, low precision on SARS-CoV-2;
it ranks candidates and does not classify them.
"""
import csv
import os
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c


def deepimmuno_name(allele):
    return allele.replace(":", "")          # HLA-A*02:01 -> HLA-A*0201


def main(fasta):
    seq = list(c.read_fasta(fasta).values())[0]
    peps = sorted({p for _, k, p in c.tile_many(seq, [9, 10])})
    path = os.path.join(c.RESULTS, "deepimmuno_input.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        # No header row on purpose. DeepImmuno reads this file with header=None.
        for a in c.CLASS_I_ALLELES:
            for p in peps:
                w.writerow([p, deepimmuno_name(a)])
    print("wrote", path, "(%d rows, no header)" % (len(peps) * len(c.CLASS_I_ALLELES)))
    print("Next, from inside the DeepImmuno repository root:")
    print("  python deepimmuno-cnn.py --mode multiple --intdir %s --outdir %s"
          % (os.path.abspath(path), os.path.abspath(c.RESULTS)))
    print("Result: deepimmuno-cnn-result.txt, tab separated, score in column 'immunogenicity'.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta")
