"""Exercise 6 - PRIME 2.0 (class I immunogenicity ranking on top of MixMHCpred 2.2).

PRIME takes the presentation score from MixMHCpred and adds TCR-facing residue
preferences (28-input fully connected net, see the slide). It ranks candidates; it does not
classify them. The TCR is never modelled.

Setup (after the MixMHCpred 2.2 setup in the README):
    git clone --depth 1 --branch v2.0 https://github.com/GfellerLab/PRIME.git PRIME2.0
    cd PRIME2.0 && g++ -O3 lib/PRIME.cc -o lib/PRIME.x
    # edit line 30 of PRIME: lib_path="/full/path/to/PRIME2.0/lib"; chmod +x PRIME lib/PRIME.x
    export PRIME=/full/path/to/PRIME2.0/PRIME
    export MIXMHCPRED=/full/path/to/MixMHCpred2.2/MixMHCpred
Note: the default branch is PRIME 2.1, which needs MixMHCpred 3.0; the slide is PRIME 2.0.

Usage:
    python ex6_prime.py fhaC.fasta
Output: results/prime_class_I.csv with columns
    start, length, peptide, allele, prime_rank, prime_score, binding_rank, binder_call, prime_top1pct
binding_rank is the MixMHCpred %rank. binder_call uses the slide class I cutoffs (SB < 0.5,
WB < 2). prime_top1pct is an arbitrary teaching cutoff (PRIME %rank < 1); there is no official
cutoff in the PRIME README, which says to use the %rank for ranking.
"""
import os
import subprocess
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c

PRIME = os.environ.get("PRIME", "./PRIME2.0/PRIME")
MIX = os.environ.get("MIXMHCPRED", "./MixMHCpred2.2/MixMHCpred")


def short(allele):
    return allele.replace("HLA-", "").replace("*", "").replace(":", "")


def main(fasta):
    seq = list(c.read_fasta(fasta).values())[0]
    peps = c.tile_many(seq, c.CLASS_I_LENGTHS)
    pep_file = os.path.join(c.RESULTS, "peptides_class_I.txt")
    with open(pep_file, "w") as fh:
        fh.write("\n".join(p for _, _, p in peps) + "\n")
    out_file = os.path.join(c.RESULTS, "prime_raw.txt")
    alleles = ",".join(short(a) for a in c.CLASS_I_ALLELES)
    subprocess.run([PRIME, "-i", pep_file, "-o", out_file, "-a", alleles, "-mix", MIX],
                   check=True, capture_output=True, text=True)

    start_of = {p: (s, k) for s, k, p in peps}
    header, rows = None, []
    with open(out_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "Peptide":
                header = parts
                continue
            rec = dict(zip(header, parts))
            pep = rec["Peptide"]
            start, length = start_of[pep]
            for a in c.CLASS_I_ALLELES:
                m = short(a)
                pr = float(rec["%Rank_" + m])
                br = float(rec["%RankBinding_" + m])
                rows.append({"start": start, "length": length, "peptide": pep, "allele": a,
                             "prime_rank": pr, "prime_score": float(rec["Score_" + m]),
                             "binding_rank": br, "binder_call": c.classify(br, "I"),
                             "prime_top1pct": pr < 1.0})
    c.write_csv(os.path.join(c.RESULTS, "prime_class_I.csv"), rows,
                ["start", "length", "peptide", "allele", "prime_rank", "prime_score",
                 "binding_rank", "binder_call", "prime_top1pct"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta")
