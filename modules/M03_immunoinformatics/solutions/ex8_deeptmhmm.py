"""Exercise 8 - DeepTMHMM: is it a membrane protein, and which way round? NOT TESTED HERE.

DeepTMHMM runs through BioLib (cloud) or the DTU web server; it is not a plain local binary.
    pip install pybiolib
    biolib run DTU/DeepTMHMM --fasta fhaC.fasta
BioLib writes an output folder (default: biolib_results/) with predicted_topologies.3line
(header, sequence, one label per residue) and TMRs.gff3.

Usage:
    python ex8_deeptmhmm.py fhaC.fasta [path/to/predicted_topologies.3line]
If the 3line path is omitted, the script runs the biolib command above first.
Output: results/deeptmhmm_summary.txt with the residue count per label and the segments.

Per-residue labels (slide 2): S signal peptide, I inside (cytosol), O outside, M alpha-helical TM,
B beta-barrel TM, P periplasm. fhaC is an outer-membrane transporter, so expect B segments;
compare their number with what you know about the FhaC beta-barrel.
Reminder from the slides: neither SignalP nor DeepTMHMM predicts antigenicity. The 3line layout
(">header", sequence, labels) is assumed from the tool's documentation; check your own file.
"""
import itertools
import os
import subprocess
import sys
from collections import Counter

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c


def segments(labels):
    pos = 1
    for lab, grp in itertools.groupby(labels):
        n = len(list(grp))
        yield lab, pos, pos + n - 1
        pos += n


def main(fasta, three_line=None):
    if three_line is None:
        subprocess.run(["biolib", "run", "DTU/DeepTMHMM", "--fasta", fasta], check=True)
        three_line = os.path.join("biolib_results", "predicted_topologies.3line")
    lines = [l.strip() for l in open(three_line) if l.strip()]
    out = []
    for i in range(0, len(lines) - 2, 3):
        name, seq, labels = lines[i], lines[i + 1], lines[i + 2]
        out.append(name + "  length=%d" % len(seq))
        out.append("  residues per label: " + str(dict(Counter(labels))))
        for lab, a, b in segments(labels):
            if lab in "MBS":
                out.append("  %s  %d-%d" % (lab, a, b))
    text = "\n".join(out)
    with open(os.path.join(c.RESULTS, "deeptmhmm_summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta",
         sys.argv[2] if len(sys.argv) > 2 else None)
