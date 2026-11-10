"""Exercise 16 - Where three independent tools agree, and where class I meets class II.

RUNS WITH NO LICENSED SOFTWARE, once ex15 and ex3 have produced their CSVs.

Why this exercise exists
------------------------
Any single ranking is one model's opinion. Slide 21 puts three tools side by
side precisely because they are built differently: NetMHCpan is a neural
network trained on eluted ligands, MHCflurry adds an antigen processing term,
and MixMHCpred is position weight matrices from motif deconvolution with no
network at all.

When three methods that disagree about almost everything still call the same
peptide a strong binder, that is worth far more than any one of their ranks.

Then the second half does something the slides never ask for and that matters
for a vaccine. It checks whether the class I consensus peptides fall inside
the best class II regions. A stretch of protein that feeds both a CD8 and a
CD4 response is a better bet than either alone.

Run these first
---------------
    python3 ex15_netmhcpan_iedb.py     NetMHCpan class I and NetMHCIIpan
    python3 ex3_mhcflurry.py           MHCflurry class I

MixMHCpred output ships in data/precomputed.

Usage
-----
    python3 ex16_consensus.py
"""
from __future__ import annotations

import os
import sys

import pandas as pd

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "utils"))

import common as c  # noqa: E402

PRECOMPUTED = os.path.join(c.DATA, "precomputed")


def load(name: str, where: str):
    path = os.path.join(where, name)
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def main() -> int:
    buffer: list[str] = []

    def log(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        buffer.append(line)

    log("=" * 78)
    log("EXERCISE 16 - Consensus across tools, and class I meeting class II")
    log("=" * 78)

    sources = {
        "NetMHCpan": (load("netmhcpan_class_I.csv", c.RESULTS), "rank_EL"),
        "MHCflurry": (load("mhcflurry_class_I.csv", c.RESULTS)
                      if load("mhcflurry_class_I.csv", c.RESULTS) is not None
                      else load("mhcflurry_class_I.csv", PRECOMPUTED),
                      "presentation_percentile"),
        "MixMHCpred": (load("mixmhcpred_class_I.csv", c.RESULTS)
                       if load("mixmhcpred_class_I.csv", c.RESULTS) is not None
                       else load("mixmhcpred_class_I.csv", PRECOMPUTED),
                       "rank"),
    }

    log("\n[16.1] Strong binders per tool")
    sets = {}
    missing = []
    for name, (df, _) in sources.items():
        if df is None:
            missing.append(name)
            log("      %-12s MISSING, run its exercise first" % name)
            continue
        sb = df[df["call"] == "SB"]
        sets[name] = set(zip(sb.peptide, sb.allele))
        log("      %-12s %4d strong binders" % (name, len(sets[name])))

    if len(sets) < 2:
        log("\n      Need at least two tools. Missing: %s" % ", ".join(missing))
        return 1

    log("\n[16.2] Agreement")
    shared = set.intersection(*sets.values())
    log("      Strong binder in ALL %d tools: %d peptide and allele pairs"
        % (len(sets), len(shared)))
    log("")
    log("      Pairwise overlap:")
    names = list(sets)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            both = sets[a] & sets[b]
            log("        %-12s and %-12s  %3d" % (a, b, len(both)))

    rows = []
    for peptide, allele in shared:
        row = {"peptide": peptide, "allele": allele}
        for name, (df, rank_col) in sources.items():
            if df is None:
                continue
            hit = df[(df.peptide == peptide) & (df.allele == allele)]
            if hit.empty:
                continue
            row["start"] = int(hit.iloc[0]["start"])
            row["%s_rank" % name.lower()] = round(float(hit.iloc[0][rank_col]), 4)
        rows.append(row)

    consensus = pd.DataFrame(rows)
    rank_cols = [col for col in consensus.columns if col.endswith("_rank")]
    consensus["worst_rank"] = consensus[rank_cols].max(axis=1)
    consensus = consensus.sort_values("worst_rank").reset_index(drop=True)
    ordered = ["peptide", "start", "allele"] + rank_cols + ["worst_rank"]
    consensus = consensus[ordered]
    consensus.to_csv(os.path.join(c.RESULTS, "consensus_class_I.csv"), index=False)

    log("\n      Ranked by the WORST rank any tool gave them, which is the")
    log("      conservative way to read a consensus:")
    log(consensus.head(15).to_string(index=False))
    log("\n      -> results/consensus_class_I.csv")

    # ------------------------------------------------------------------ II
    log("\n[16.3] Does class I meet class II?")
    class_ii = load("netmhciipan_class_II.csv", c.RESULTS)
    if class_ii is None:
        log("      netmhciipan_class_II.csv not found, run ex15 first.")
    else:
        top_ii = class_ii[class_ii["call"] == "SB"].sort_values("rank_EL").head(12)
        log("      fhaC is extracellular, so class II is its natural route.")
        log("      A region that feeds BOTH a CD8 and a CD4 response is a better")
        log("      bet than either on its own. Checking the top class II hits:")
        log("")
        overlaps = []
        for _, r in top_ii.iterrows():
            start_ii = int(r.start)
            end_ii = start_ii + len(r.peptide) - 1
            inside = consensus[
                (consensus.start + consensus.peptide.str.len() - 1 >= start_ii)
                & (consensus.start <= end_ii)]
            tag = ""
            if not inside.empty:
                tag = "  overlaps class I: " + ", ".join(inside.peptide)
                for _, q in inside.iterrows():
                    overlaps.append({"region_start": start_ii, "region_end": end_ii,
                                     "class_II_peptide": r.peptide,
                                     "class_II_allele": r.allele,
                                     "class_II_rank": r.rank_EL,
                                     "class_I_peptide": q.peptide,
                                     "class_I_allele": q.allele})
            log("        %-16s %-10s %3d-%-3d rank=%-6s%s"
                % (r.peptide, r.allele, start_ii, end_ii, r.rank_EL, tag))

        if overlaps:
            df_ov = pd.DataFrame(overlaps)
            df_ov.to_csv(os.path.join(c.RESULTS, "class_I_II_overlap.csv"), index=False)
            regions = df_ov.groupby(["region_start", "region_end"]).size()
            log("")
            log("      %d dual regions found:" % len(regions))
            for (a, b), n in regions.items():
                peps = df_ov[(df_ov.region_start == a)].class_I_peptide.unique()
                log("        residues %d-%d, %d class I peptide(s): %s"
                    % (a, b, len(peps), ", ".join(peps)))
            log("")
            log("      -> results/class_I_II_overlap.csv")
            log("")
            log("      This is the most useful result in the module and it is not")
            log("      on any slide. Still a hypothesis, but one supported from")
            log("      two independent directions at once.")

    log("\n" + "=" * 78)
    log("QUESTIONS TO DISCUSS")
    log("=" * 78)
    log("  1. Agreement between three tools is reassuring. Could they agree and")
    log("     all be wrong together? What would they have in common that would")
    log("     cause that?")
    log("  2. We ranked by the WORST rank across tools. What changes if you rank")
    log("     by the average instead, and which is the honest default?")
    log("  3. A dual class I and class II region looks ideal. What experiment")
    log("     would you run first to test it?")
    log("  4. MixMHCpred is not a neural network. Does that make its agreement")
    log("     with the other two worth more, or less?")

    with open(os.path.join(c.RESULTS, "ex16_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(buffer))
    print("\n[log] results/ex16_log.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
