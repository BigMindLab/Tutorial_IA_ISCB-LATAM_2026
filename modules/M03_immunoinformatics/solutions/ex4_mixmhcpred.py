"""Exercise 4 - MixMHCpred 2.2 (class I, PWM motifs from MS ligands; not a neural network).

Setup (Linux/Mac):
    git clone --depth 1 --branch v2.2 https://github.com/GfellerLab/MixMHCpred.git MixMHCpred2.2
    cd MixMHCpred2.2
    g++ -O3 lib/MixMHCpred.cc -o lib/MixMHCpred.x
    # edit line 27 of the MixMHCpred script: lib_path="/full/path/to/MixMHCpred2.2/lib"
    chmod +x MixMHCpred lib/MixMHCpred.x
    export MIXMHCPRED=/full/path/to/MixMHCpred2.2/MixMHCpred
Note: the default branch of that repo is v3.0 (a neural network), NOT 2.2.

Usage:
    python ex4_mixmhcpred.py fhaC.fasta
Output: results/mixmhcpred_class_I.csv with columns
    start, length, peptide, allele, score, rank, call
MixMHCpred outputs a score and a %rank, not an affinity. Use %rank.
"""
import os
import subprocess
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c

BIN = os.environ.get("MIXMHCPRED", "./MixMHCpred2.2/MixMHCpred")


def mix_name(allele):
    """HLA-A*02:01 -> A0201"""
    return allele.replace("HLA-", "").replace("*", "").replace(":", "")


def main(fasta):
    seq = list(c.read_fasta(fasta).values())[0]
    peps = c.tile_many(seq, c.CLASS_I_LENGTHS)
    pep_file = os.path.join(c.RESULTS, "peptides_class_I.txt")   # no special chars in names
    with open(pep_file, "w") as fh:
        fh.write("\n".join(p for _, _, p in peps) + "\n")
    out_file = os.path.join(c.RESULTS, "mixmhcpred_raw.txt")
    alleles = ",".join(mix_name(a) for a in c.CLASS_I_ALLELES)
    subprocess.run([BIN, "-i", pep_file, "-o", out_file, "-a", alleles],
                   check=True, capture_output=True, text=True)

    start_of = {p: (s, k) for s, k, p in peps}
    # Output has comment lines starting with '#', then a header, then one row per peptide:
    # Peptide BestAllele %Rank_bestAllele Score_bestAllele  then, per allele, Score_X and %Rank_X.
    header, rows = None, []
    with open(out_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if parts[0] == "Peptide":
                header = parts
                continue
            rec = dict(zip(header, parts))
            pep = rec["Peptide"]
            start, length = start_of[pep]
            for a in c.CLASS_I_ALLELES:
                m = mix_name(a)
                rank, score = float(rec["%Rank_" + m]), float(rec["Score_" + m])
                rows.append({"start": start, "length": length, "peptide": pep,
                             "allele": a, "score": score, "rank": rank,
                             "call": c.classify(rank, "I")})
    c.write_csv(os.path.join(c.RESULTS, "mixmhcpred_class_I.csv"), rows,
                ["start", "length", "peptide", "allele", "score", "rank", "call"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta")
