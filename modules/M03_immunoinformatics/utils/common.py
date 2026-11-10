"""Shared helpers for the MHC prediction exercises (ASCII only, US English)."""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
# Results live at the repository root, one level above scripts/, so that every
# exercise writes to the same place regardless of where it is launched from.
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
os.makedirs(RESULTS, exist_ok=True)

# Default input for every exercise: the fhaC sequence sent to participants.
DEFAULT_FASTA = os.path.join(DATA, "fhaC.fasta")

# Editable allele panel. Class I alleles feed NetMHCpan, MHCflurry, MixMHCpred.
# Class II alleles feed NetMHCIIpan only.
CLASS_I_ALLELES = ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"]
CLASS_II_ALLELES = ["DRB1_0401", "DRB1_1501"]   # NetMHCIIpan naming
CLASS_I_LENGTHS = [8, 9, 10, 11]                # slide range is 8-14; trim for speed
CLASS_II_LENGTH = 15

# %rank thresholds taken from the slides
THRESHOLDS = {
    "I":  {"SB": 0.5, "WB": 2.0},
    "II": {"SB": 2.0, "WB": 10.0},
}


def read_fasta(path):
    """Return dict name -> sequence. Handles CRLF line endings."""
    seqs, name = {}, None
    with open(path, newline="") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            elif name is not None:
                seqs[name] += line.upper()
    return seqs


def tile(seq, length):
    """All overlapping k-mers of one length. Returns list of (start_1based, peptide)."""
    return [(i + 1, seq[i:i + length]) for i in range(len(seq) - length + 1)]


def tile_many(seq, lengths):
    out = []
    for k in lengths:
        out.extend((s, k, p) for s, p in tile(seq, k))
    return out


def classify(rank, mhc_class):
    t = THRESHOLDS[mhc_class]
    if rank < t["SB"]:
        return "SB"
    if rank < t["WB"]:
        return "WB"
    return "NB"


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    # Print a path relative to the repository, never the absolute one, so that
    # logs committed to the repository do not carry anybody's local directories.
    try:
        shown = os.path.relpath(path, os.path.join(HERE, ".."))
    except ValueError:
        shown = os.path.basename(path)
    print("wrote", shown.replace(os.sep, "/"), "(%d rows)" % len(rows))
