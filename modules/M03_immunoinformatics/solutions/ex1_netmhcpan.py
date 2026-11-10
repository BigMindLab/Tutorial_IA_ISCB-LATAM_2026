"""Exercise 1 - NetMHCpan 4.x (class I, EL + BA). NOT TESTED HERE: needs a local install.

Get the Linux binary from https://services.healthtech.dtu.dk/ (academic license,
manual download), then put `netMHCpan` on your PATH or set NETMHCPAN_BIN.

Usage:
    python ex1_netmhcpan.py fhaC.fasta
Output: results/netmhcpan_class_I.csv with columns
    start, length, peptide, allele, score_EL, rank_EL, score_BA, rank_BA, affinity_nM, call

Parsing assumes the NetMHCpan 4.1 stdout table printed with -BA:
    Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM) BindLevel
Check one line of your own output against this if your version differs.
"""
import os
import subprocess
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c

BIN = os.environ.get("NETMHCPAN_BIN", "netMHCpan")


def run_one(pep_file, allele):
    # NetMHCpan wants names without the asterisk, e.g. HLA-A02:01
    cmd = [BIN, "-p", pep_file, "-a", allele.replace("*", ""), "-BA"]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def parse(stdout, allele, start_of):
    rows = []
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) < 16 or not parts[0].isdigit() or not parts[1].startswith("HLA-"):
            continue
        pep = parts[2]
        try:
            score_el, rank_el = float(parts[11]), float(parts[12])
            score_ba, rank_ba, aff = float(parts[13]), float(parts[14]), float(parts[15])
        except ValueError:
            continue
        start, length = start_of[pep]
        rows.append({
            "start": start, "length": length, "peptide": pep, "allele": allele,
            "score_EL": score_el, "rank_EL": rank_el,
            "score_BA": score_ba, "rank_BA": rank_ba, "affinity_nM": aff,
            "call": c.classify(rank_el, "I"),
        })
    return rows


def main(fasta):
    seq = list(c.read_fasta(fasta).values())[0]
    peps = c.tile_many(seq, c.CLASS_I_LENGTHS)
    pep_file = os.path.join(c.RESULTS, "peptides_class_I.txt")
    with open(pep_file, "w") as fh:
        fh.write("\n".join(p for _, _, p in peps) + "\n")
    start_of = {p: (s, k) for s, k, p in peps}   # first occurrence wins
    rows = []
    for allele in c.CLASS_I_ALLELES:
        rows += parse(run_one(pep_file, allele), allele, start_of)
    fields = ["start", "length", "peptide", "allele", "score_EL", "rank_EL",
              "score_BA", "rank_BA", "affinity_nM", "call"]
    c.write_csv(os.path.join(c.RESULTS, "netmhcpan_class_I.csv"), rows, fields)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta")
