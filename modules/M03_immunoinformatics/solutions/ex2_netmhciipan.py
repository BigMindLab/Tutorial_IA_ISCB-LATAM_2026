"""Exercise 2 - NetMHCIIpan 4.x (class II, 15-mers). NOT TESTED HERE: needs a local install.

Get the Linux binary from https://services.healthtech.dtu.dk/ (academic license,
manual download), then put `netMHCIIpan` on your PATH or set NETMHCIIPAN_BIN.

Usage:
    python ex2_netmhciipan.py fhaC.fasta
Output: results/netmhciipan_class_II.csv with columns
    start, peptide, allele, score_EL, rank_EL, call
Class II slide thresholds: SB < 2 %rank, WB < 10 %rank.

The parser is header-driven: it finds the line starting with "Pos" and looks for the
columns named Peptide, Score_EL (or Score) and %Rank_EL (or %Rank), so small
differences between versions are tolerated. If it finds nothing, open the raw
stdout saved in results/netmhciipan_raw_<allele>.txt and adjust the names below.
"""
import os
import subprocess
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c

BIN = os.environ.get("NETMHCIIPAN_BIN", "netMHCIIpan")
RANK_NAMES = ("%Rank_EL", "%Rank")
SCORE_NAMES = ("Score_EL", "Score")


def run_one(pep_file, allele):
    # -inptype 1 = peptide list; each line is already a 15-mer
    cmd = [BIN, "-f", pep_file, "-inptype", "1", "-a", allele]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def pick(header, names):
    for n in names:
        if n in header:
            return header.index(n)
    return None


def parse(stdout, allele, start_of):
    rows, header = [], None
    for line in stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == "Pos" and "Peptide" in parts:
            header = parts
            continue
        if header is None or not parts or not parts[0].isdigit():
            continue
        i_pep = header.index("Peptide")
        i_rank, i_score = pick(header, RANK_NAMES), pick(header, SCORE_NAMES)
        if i_rank is None or len(parts) <= max(i_pep, i_rank):
            continue
        pep = parts[i_pep]
        rank = float(parts[i_rank])
        score = float(parts[i_score]) if i_score is not None else float("nan")
        rows.append({"start": start_of.get(pep, 0), "peptide": pep, "allele": allele,
                     "score_EL": score, "rank_EL": rank,
                     "call": c.classify(rank, "II")})
    return rows


def main(fasta):
    seq = list(c.read_fasta(fasta).values())[0]
    peps = c.tile(seq, c.CLASS_II_LENGTH)
    pep_file = os.path.join(c.RESULTS, "peptides_class_II.txt")
    with open(pep_file, "w") as fh:
        fh.write("\n".join(p for _, p in peps) + "\n")
    start_of = {p: s for s, p in peps}
    rows = []
    for allele in c.CLASS_II_ALLELES:
        out = run_one(pep_file, allele)
        with open(os.path.join(c.RESULTS, "netmhciipan_raw_%s.txt" % allele), "w") as fh:
            fh.write(out)
        rows += parse(out, allele, start_of)
    c.write_csv(os.path.join(c.RESULTS, "netmhciipan_class_II.csv"), rows,
                ["start", "peptide", "allele", "score_EL", "rank_EL", "call"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta")
