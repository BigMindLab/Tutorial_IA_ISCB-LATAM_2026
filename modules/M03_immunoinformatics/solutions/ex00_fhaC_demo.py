"""Exercise 00 - START HERE. The whole fhaC story, with no training anywhere.

This is the script to show in the room.

It runs the complete workflow of slides 18 to 24 on fhaC and nothing else.
No model is trained while it runs. Every number it prints comes from YOUR
protein. It takes a few minutes, almost all of it MHCflurry scoring peptides.

What it chains
--------------
    exercise 11   slides 18 and 19   where does fhaC go
    exercise 14   slides 20 to 24    binding, immunogenicity, shortlist

Those two exercises are already training free, so this runner simply calls
them in order with clear headers. There is no duplicated logic here, which
means the demo can never drift away from the real exercises.

Where the model comes from
--------------------------
Exercise 14 loads pretrained weights from
data/precomputed/ex12_immunogenicity_cnn.pt, which ship with the repository.
Those weights were produced once by exercise 12 and committed, exactly the way
MHCflurry ships its own. Nothing is fitted at run time.

The complementary track
-----------------------
Exercises 12 and 13 are the ones that train. They exist to show what is inside
the architectures on slides 22 and 23, and to make two points that cannot be
made by prediction alone:

    exercise 12   presentation is not immunogenicity, as two AUC numbers
    exercise 13   NetTCR must not be used on fhaC, as a collapsing AUC

Run those when there is time for the methodology discussion. Run this one when
the goal is the biology of fhaC.

Usage
-----
    python ex00_fhaC_demo.py
"""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODULE = os.path.normpath(os.path.join(HERE, ".."))
DEEPTMHMM = os.path.join(MODULE, "data", "precomputed", "deeptmhmm_fhaC.3line")
MODEL = os.path.join(MODULE, "data", "precomputed", "ex12_immunogenicity_cnn.pt")

STEPS = [
    (
        "PART 1 of 2   Where does fhaC go?   (slides 18 and 19)",
        "Signal peptide and membrane topology, scored against the curated "
        "UniProt annotation of P35077.",
        [sys.executable, os.path.join(HERE, "ex11_topology_reference.py"),
         "--deeptmhmm", DEEPTMHMM],
    ),
    (
        "PART 2 of 2   From peptides to a shortlist   (slides 20 to 24)",
        "Binding and presentation with MHCflurry, immunogenicity from the "
        "pretrained model, then the pVACtools style filter chain.",
        [sys.executable, os.path.join(HERE, "ex14_prioritization.py")],
    ),
]


def banner(title: str, subtitle: str) -> None:
    print()
    print("#" * 78)
    print("#  " + title)
    print("#  " + subtitle)
    print("#" * 78)
    sys.stdout.flush()


def check_dependencies() -> bool:
    required = {
        "numpy": "numpy", "pandas": "pandas", "matplotlib": "matplotlib",
        "sklearn": "scikit-learn", "torch": "torch", "mhcflurry": "mhcflurry",
    }
    missing = [package for module, package in required.items()
               if importlib.util.find_spec(module) is None]
    if not missing:
        return True

    print("[ERROR] Missing packages: " + ", ".join(missing), file=sys.stderr)
    print("Python interpreter: " + sys.executable, file=sys.stderr)
    print("From the module directory, activate the tutorial environment:", file=sys.stderr)
    print("  source .venv/bin/activate", file=sys.stderr)
    print("  python solutions/ex00_fhaC_demo.py", file=sys.stderr)
    print("Or select its Python directly:", file=sys.stderr)
    print("  .venv/bin/python solutions/ex00_fhaC_demo.py", file=sys.stderr)
    print("For Conda, use: conda activate m03-immunoinformatics", file=sys.stderr)
    print("If packages are missing inside the active environment, run:", file=sys.stderr)
    print("  bash setup.sh --use-current-env", file=sys.stderr)
    return False


def main() -> int:
    if not check_dependencies():
        return 1

    print("=" * 78)
    print("fhaC DEMO - the whole workflow on one protein, no training")
    print("=" * 78)
    print()
    print("  Protein : fhaC, filamentous hemagglutinin transporter")
    print("            Bordetella pertussis, UniProtKB P35077, 584 residues")
    print("  Training: none. Pretrained weights are loaded from disk.")
    print("  Python  : " + sys.executable)
    print()
    print("  A note on counting, because it confuses people every time.")
    print("  Part 1 reports 584 rows, one per residue. Part 2 reports 576")
    print("  peptides, because 584 residues hold 576 overlapping 9-mers,")
    print("  and 2880 pairs, because each peptide is scored against 5 HLA")
    print("  alleles. Same protein throughout.")

    if not os.path.exists(MODEL):
        print()
        print("  [ERROR] Missing pretrained weights:")
        print("          data/precomputed/ex12_immunogenicity_cnn.pt")
        print("          Either restore that file or run ex12_immunogenicity.py once.")
        return 1

    if not os.path.exists(DEEPTMHMM):
        print()
        print("  [note] No DeepTMHMM output found, part 1 will show the UniProt")
        print("         reference only.")

    for title, subtitle, command in STEPS:
        banner(title, subtitle)
        result = subprocess.run(command, cwd=HERE)
        if result.returncode != 0:
            print("\n[ERROR] step failed: %s" % title)
            return result.returncode

    print()
    print("=" * 78)
    print("DONE. What you have now")
    print("=" * 78)
    print("  results/topology_reference.csv   per residue reference for fhaC")
    print("  results/ex14_shortlist.csv       the ranked candidate list")
    print("  results/ex14_full_table.csv      every peptide and allele pair")
    print("  results/ex14_prioritization.png  the funnel and the ranking")
    print()
    print("  None of this proves a T cell response. It is a ranked list of")
    print("  hypotheses, which is the honest output of the whole module.")
    print()
    print("  For the methodology behind the models, run ex12_immunogenicity.py")
    print("  and ex13_nettcr.py. Those two train, and that is their point.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
