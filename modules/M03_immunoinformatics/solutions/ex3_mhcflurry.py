"""Exercise 3 - MHCflurry 2.x (class I: affinity + antigen processing + presentation).

Setup:
    pip install mhcflurry
    mhcflurry-downloads fetch models_class1_presentation
Usage:
    python ex3_mhcflurry.py fhaC.fasta
Output: results/mhcflurry_class_I.csv with columns
    start, length, peptide, allele, affinity_nM, affinity_percentile,
    processing_score, presentation_score, presentation_percentile, call

call uses affinity_percentile with the class I slide thresholds (SB < 0.5, WB < 2).
Compare affinity_percentile (BA-style) with presentation_percentile (BA + processing)
to see the "binding is not presentation" point from the slides.
"""
import sys

import pandas as pd

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c


def main(fasta):
    from mhcflurry import Class1PresentationPredictor

    seq = list(c.read_fasta(fasta).values())[0]
    peps = c.tile_many(seq, c.CLASS_I_LENGTHS)
    info = {p: (s, k) for s, k, p in peps}
    unique = sorted(info)

    predictor = Class1PresentationPredictor.load()
    frames = []
    for allele in c.CLASS_I_ALLELES:
        # One allele at a time so each row is that allele, not the best of a genotype.
        df = predictor.predict(peptides=unique, alleles=[allele],
                               include_affinity_percentile=True, verbose=0)
        df = df.copy()
        df["allele"] = allele
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["start"] = out["peptide"].map(lambda p: info[p][0])
    out["length"] = out["peptide"].map(lambda p: info[p][1])
    out = out.rename(columns={"affinity": "affinity_nM"})
    out["call"] = out["affinity_percentile"].map(lambda r: c.classify(r, "I"))
    cols = ["start", "length", "peptide", "allele", "affinity_nM", "affinity_percentile",
            "processing_score", "presentation_score", "presentation_percentile", "call"]
    out = out[cols].sort_values(["allele", "start", "length"])
    out.to_csv("%s/mhcflurry_class_I.csv" % c.RESULTS, index=False)
    print("wrote results/mhcflurry_class_I.csv (%d rows)" % len(out))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fhaC.fasta")
