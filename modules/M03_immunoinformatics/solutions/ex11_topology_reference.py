"""Exercise 11 - Slides 18 and 19: scoring SignalP 6.0 and DeepTMHMM against curated annotation.

RUNS WITH NO LICENSED SOFTWARE. This is the exercise that works the moment you
clone the repository, and it becomes richer once SignalP and DeepTMHMM are
installed.

Why this exercise exists
------------------------
Slides 18 and 19 present two predictors and end with the same warning: neither
one predicts antigenicity. Running them is easy. Knowing whether to believe
them is the actual skill, and that needs something to compare against.

fhaC is a good case because UniProt P35077 carries curated annotation derived
from experimental work, including a crystal structure: a signal peptide, a
POTRA domain and sixteen transmembrane beta strands. That annotation is the
reference used here.

The reference is curated annotation, not ground truth. Where a predictor and
UniProt disagree, the interesting question is which one to trust and why, and
that is the discussion this exercise is built for.

Usage
-----
    python ex11_topology_reference.py

    # once you have run the other two, pass their output to compare
    python ex11_topology_reference.py --signalp results/signalp6/prediction_results.txt
    python ex11_topology_reference.py --deeptmhmm biolib_results/predicted_topologies.3line

Output: results/topology_reference.csv and a printed report.
"""
from __future__ import annotations

import argparse
import json
import os
from collections import Counter

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c

HERE = os.path.dirname(os.path.abspath(__file__))
UNIPROT_JSON = os.path.join(HERE, "..", "data", "P35077.json")

# DeepTMHMM per residue alphabet, from slide 19.
LABEL_HELP = {
    "S": "signal peptide",
    "I": "inside, cytosol",
    "O": "outside",
    "M": "alpha helical transmembrane",
    "B": "beta barrel transmembrane",
    "P": "periplasm",
}


def load_uniprot_regions(path):
    """Pull the topology features out of a UniProt entry in JSON form."""
    with open(path, encoding="utf-8") as fh:
        entry = json.load(fh)
    wanted = {"Signal", "Transmembrane", "Domain", "Chain"}
    regions = []
    for feat in entry.get("features", []):
        if feat.get("type") not in wanted:
            continue
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        regions.append({
            "type": feat["type"],
            "start": start,
            "end": end,
            "description": feat.get("description", ""),
        })
    return sorted(regions, key=lambda r: (r["start"], r["end"]))


def reference_labels(regions, length):
    """Build a per residue label string in the DeepTMHMM alphabet.

    This is an approximation on purpose, and the approximation is part of the
    lesson. UniProt annotates a signal peptide and the beta strands. It does
    not label every residue as inside or outside, so everything that is not a
    signal peptide or a strand is left as 'O'. A loop facing the periplasm is
    therefore marked 'O' here even though DeepTMHMM may well call it 'P'.
    Disagreement on those residues is expected and is not a predictor error.
    """
    labels = ["O"] * length
    for r in regions:
        if r["type"] == "Signal":
            for i in range(r["start"] - 1, min(r["end"], length)):
                labels[i] = "S"
    for r in regions:
        if r["type"] == "Transmembrane":
            for i in range(r["start"] - 1, min(r["end"], length)):
                labels[i] = "B"
    return "".join(labels)


def parse_signalp(path):
    """Read SignalP 6.0 prediction_results.txt and return the first data row."""
    fields = ["id", "prediction", "OTHER", "SP_Sec_SPI", "LIPO_Sec_SPII",
              "TAT_Tat_SPI", "TATLIPO_Tat_SPII", "PILIN_Sec_SPIII", "cs_position"]
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            parts += [""] * (len(fields) - len(parts))
            return dict(zip(fields, parts[:len(fields)]))
    return None


def cleavage_position(cs_text):
    """SignalP writes something like 'CS pos: 30-31. Pr: 0.87'. Pull out 30."""
    if not cs_text or "pos:" not in cs_text:
        return None
    try:
        chunk = cs_text.split("pos:")[1].split(".")[0].strip()
        return int(chunk.split("-")[0])
    except (IndexError, ValueError):
        return None


def parse_3line(path):
    """Read a DeepTMHMM predicted_topologies.3line file into (name, seq, labels)."""
    lines = [l.strip() for l in open(path, encoding="utf-8") if l.strip()]
    out = []
    for i in range(0, len(lines) - 2, 3):
        out.append((lines[i], lines[i + 1], lines[i + 2]))
    return out


def segments(labels):
    """Yield (label, start_1based, end_1based) for each run of equal labels."""
    import itertools
    pos = 1
    for lab, grp in itertools.groupby(labels):
        n = len(list(grp))
        yield lab, pos, pos + n - 1
        pos += n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", default=os.path.join(HERE, "..", "data", "fhaC.fasta"))
    ap.add_argument("--signalp", default=None,
                    help="SignalP 6.0 prediction_results.txt, optional")
    ap.add_argument("--deeptmhmm", default=None,
                    help="DeepTMHMM predicted_topologies.3line, optional")
    args = ap.parse_args()

    seq = list(c.read_fasta(args.fasta).values())[0]
    print("=" * 74)
    print("EXERCISE 11 - Topology reference for fhaC (slides 18 and 19)")
    print("=" * 74)
    print()
    print("  NO TRAINING HAPPENS IN THIS EXERCISE. Everything below is YOUR")
    print("  protein, fhaC, compared against its curated UniProt annotation.")
    print()
    print("  Note on counting: this exercise works RESIDUE BY RESIDUE, so it")
    print("  reports %d rows. Exercises 12 and 14 work PEPTIDE BY PEPTIDE and" % len(seq))
    print("  report %d, because %d residues hold %d overlapping 9-mers."
          % (len(seq) - 8, len(seq), len(seq) - 8))
    print("  Same protein, different unit of analysis.")
    print()
    print("sequence length: %d residues" % len(seq))

    regions = load_uniprot_regions(UNIPROT_JSON)
    signal = [r for r in regions if r["type"] == "Signal"]
    strands = [r for r in regions if r["type"] == "Transmembrane"]
    domains = [r for r in regions if r["type"] == "Domain"]

    print("\n[reference] UniProt P35077 curated annotation")
    for r in signal:
        print("  signal peptide      %d-%d" % (r["start"], r["end"]))
    for r in domains:
        print("  domain              %d-%d  %s" % (r["start"], r["end"], r["description"]))
    print("  transmembrane       %d beta strands, from %d to %d"
          % (len(strands), strands[0]["start"], strands[-1]["end"]) if strands else "")

    labels = reference_labels(regions, len(seq))
    rows = [{"position": i + 1, "residue": seq[i], "reference_label": labels[i]}
            for i in range(len(seq))]
    c.write_csv(os.path.join(c.RESULTS, "topology_reference.csv"), rows,
                ["position", "residue", "reference_label"])
    print("\n  per residue reference written, label counts: %s" % dict(Counter(labels)))

    # DeepTMHMM also labels a signal peptide, so read it once up front and use
    # it in both questions below.
    dtm_records = []
    if args.deeptmhmm and os.path.exists(args.deeptmhmm):
        dtm_records = parse_3line(args.deeptmhmm)

    # ---------------------------------------------------------------- SignalP
    print("\n[question 1] Does fhaC enter the secretory pathway? (slide 18)")
    ref_cs = signal[0]["end"] if signal else None

    if dtm_records:
        _, _, dtm_labels = dtm_records[0]
        s_segments = [(a, b) for lab, a, b in segments(dtm_labels) if lab == "S"]
        if s_segments:
            dtm_cs = s_segments[0][1]
            print("  DeepTMHMM also predicts a signal peptide, ending at %d." % dtm_cs)
            print("  UniProt annotates it ending at %s." % ref_cs)
            if ref_cs:
                print("  difference: %+d residues" % (dtm_cs - ref_cs))
            print("  Two independent predictions of the same boundary is already")
            print("  something to compare, before SignalP is even installed.")
            print()

    if args.signalp and os.path.exists(args.signalp):
        row = parse_signalp(args.signalp)
        pred_cs = cleavage_position(row.get("cs_position", "")) if row else None
        print("  SignalP prediction : %s" % (row.get("prediction") if row else "unreadable"))
        print("  SignalP cleavage   : %s" % (pred_cs if pred_cs else "not reported"))
        print("  UniProt cleavage   : %s" % ref_cs)
        if pred_cs and ref_cs:
            delta = pred_cs - ref_cs
            verdict = "exact match" if delta == 0 else "off by %+d residues" % delta
            print("  agreement          : %s" % verdict)
            print("\n  Discuss: a cleavage site off by one or two residues changes which")
            print("  peptides exist at the N terminus of the mature protein. Check whether")
            print("  any of your class I candidates straddle that boundary.")
    else:
        print("  SignalP output not supplied, so only the reference is shown.")
        print("  UniProt annotates a signal peptide at 1-%s." % ref_cs)
        print("  Run ex7_signalp6.py once SignalP 6.0 is installed, then pass")
        print("  --signalp results/signalp6/prediction_results.txt to compare.")
        print("\n  Predict before you run it: fhaC is from a Gram negative bacterium,")
        print("  so the organism argument is 'other', not 'eukarya'. Which of the five")
        print("  signal peptide types on slide 18 do you expect, and why?")

    # -------------------------------------------------------------- DeepTMHMM
    print("\n[question 2] Is it a membrane protein, and which way round? (slide 19)")
    print("  UniProt annotates %d transmembrane beta strands." % len(strands))
    if dtm_records:
        for name, pseq, plab in dtm_records:
            print("  %s  length=%d" % (name, len(pseq)))
            counts = dict(Counter(plab))
            print("  label counts: %s" % counts)
            pred_strands = [(a, b) for lab, a, b in segments(plab) if lab == "B"]
            print("  predicted beta strands: %d" % len(pred_strands))
            print("  UniProt beta strands  : %d" % len(strands))
            if len(plab) == len(labels):
                agree = sum(1 for x, y in zip(plab, labels) if x == y)
                print("  raw per residue agreement: %.1f%%" % (100.0 * agree / len(labels)))
                print("  That number looks poor, and it is misleading. The reference")
                print("  marks every non strand residue as 'O' because UniProt does not")
                print("  label inside and outside, so every 'P' DeepTMHMM predicts is")
                print("  counted as a mismatch even when it is right.")
                # The fair comparison is the question UniProt actually answers:
                # is this residue inside a transmembrane strand or not?
                strand_agree = sum(1 for x, y in zip(plab, labels)
                                   if (x == "B") == (y == "B"))
                print("  agreement on the strand or not strand question: %.1f%%"
                      % (100.0 * strand_agree / len(labels)))
                print("  That is the comparison the reference can actually support.")
            if len(pred_strands) % 2 != 0:
                print("  Note: an odd number of strands. Slide 19 says the state model")
                print("  forces an even number, so check how the tool resolved this.")
    else:
        print("  DeepTMHMM output not supplied.")
        print("  It runs through BioLib rather than as a local binary:")
        print("      pip install pybiolib")
        print("      biolib run DTU/DeepTMHMM --fasta data/fhaC.fasta")
        print("  BioLib now requires an account, and running locally is no longer")
        print("  supported, so expect to sign in first.")
        print("\n  Predict before you run it: fhaC is an outer membrane transporter.")
        print("  Do you expect 'M' labels, 'B' labels, or both? Slide 19 notes the")
        print("  model was built mainly for prokaryotic beta barrels.")

    # ------------------------------------------------------------ the bridge
    print("\n[why this matters downstream]")
    if ref_cs:
        print("  The first %d residues are cut off during secretion, so they are not" % ref_cs)
        print("  part of the mature protein. Any class I or class II candidate that")
        print("  falls entirely inside 1-%d should be dropped before ranking." % ref_cs)
        for k in (9, 15):
            n_removed = max(0, ref_cs - k + 1)
            print("    %d-mers fully inside the signal peptide: %d" % (k, n_removed))
    print("\n  That is the first filter of the prioritisation exercise, and it comes")
    print("  from slide 18, not from any binding predictor.")

    print("\n" + "=" * 74)
    print("QUESTIONS TO DISCUSS")
    print("=" * 74)
    print("  1. SignalP gives a signal peptide type and a cleavage site. Neither")
    print("     says anything about antigenicity. What does each one actually")
    print("     constrain in the rest of this workflow?")
    print("  2. The reference here is curated annotation, much of it from a solved")
    print("     structure. If a predictor disagrees with it, is the predictor wrong?")
    print("  3. DeepTMHMM labels every residue. SignalP labels only the N terminus.")
    print("     Why does that difference matter when you chain them together?")
    print("  4. fhaC is an extracellular bacterial protein. Which MHC class is the")
    print("     natural route for it, and what does that imply for exercises that")
    print("     only predict class I? (Look back at slide 20.)")


if __name__ == "__main__":
    main()
