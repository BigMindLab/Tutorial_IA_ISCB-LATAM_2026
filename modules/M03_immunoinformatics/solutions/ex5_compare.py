"""Exercise 5 - Compare tools, call SB / WB / NB, and add immunogenicity and localization context.

Reads whatever result files exist in results/ (run ex1-ex4, ex6, ex7 first; missing ones are skipped):
    netmhcpan_class_I.csv    rank_EL                      -> binding call (SB/WB/NB)
    mhcflurry_class_I.csv    affinity_ and presentation_percentile -> binding call
    mixmhcpred_class_I.csv   rank                         -> binding call
    prime_class_I.csv        prime_rank                   -> rank only, NO SB/WB call
    signalp6_summary.csv     cleavage site                -> flags peptides in the signal peptide
    netmhciipan_class_II.csv class II, reported on its own

PRIME is kept out of the SB/WB calls and the consensus vote on purpose: its README gives no
official cutoff and says to use the %rank for ranking. The prime_top1pct column (%rank < 1) is an
arbitrary teaching cutoff. SignalP's flag is information, not a filter: peptides overlapping the
signal peptide are marked in_signal_peptide so you can discuss them, not dropped.

Outputs:
    results/class_I_merged.csv      one row per peptide x allele, one rank column per tool
    results/class_I_consensus.csv   peptides called SB by >= 2 binding tools (or 1, if only one ran)
    results/summary.txt             counts, rank correlations, agreement, PRIME and SignalP notes

Usage:
    python ex5_compare.py
"""
import os
import re

import pandas as pd

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "utils"))
import common as c

R = c.RESULTS
SOURCES = {   # binding tool label -> (file, rank column)
    "NetMHCpan_EL":       ("netmhcpan_class_I.csv", "rank_EL"),
    "MHCflurry_affinity": ("mhcflurry_class_I.csv", "affinity_percentile"),
    "MHCflurry_present":  ("mhcflurry_class_I.csv", "presentation_percentile"),
    "MixMHCpred":         ("mixmhcpred_class_I.csv", "rank"),
}


def load_binding():
    merged, used = None, []
    for label, (fname, col) in SOURCES.items():
        path = os.path.join(R, fname)
        if not os.path.exists(path):
            print("skipping %s (no %s)" % (label, fname))
            continue
        df = pd.read_csv(path)[["start", "length", "peptide", "allele", col]]
        df = df.drop_duplicates(["peptide", "allele"]).rename(columns={col: label})
        merged = df if merged is None else merged.merge(
            df.drop(columns=["start", "length"]), on=["peptide", "allele"], how="outer")
        used.append(label)
    return merged, used


def add_prime(merged):
    path = os.path.join(R, "prime_class_I.csv")
    if not os.path.exists(path):
        print("skipping PRIME (no prime_class_I.csv)")
        return merged, False
    p = pd.read_csv(path)[["peptide", "allele", "prime_rank"]].drop_duplicates(["peptide", "allele"])
    merged = merged.merge(p, on=["peptide", "allele"], how="left")
    merged["prime_top1pct"] = merged["prime_rank"] < 1.0
    return merged, True


def signalp_cleavage():
    """Return the position of the last residue of the signal peptide, or None.
    Assumes the cs_position text contains 'CS pos: A-B' (SignalP 6.0 layout); confirm on your output."""
    path = os.path.join(R, "signalp6_summary.csv")
    if not os.path.exists(path):
        print("skipping SignalP (no signalp6_summary.csv)")
        return None
    df = pd.read_csv(path)
    if df.empty or str(df.loc[0, "prediction"]).upper() == "OTHER":
        return 0   # no signal peptide predicted
    m = re.search(r"CS pos:\s*(\d+)", str(df.loc[0, "cs_position"]))
    return int(m.group(1)) if m else None


def main():
    merged, used = load_binding()
    if merged is None:
        print("no class I binding results found; run ex1, ex3 or ex4 first")
        return
    lines = []
    for t in used:
        merged["call_" + t] = merged[t].map(lambda r: c.classify(r, "I") if pd.notna(r) else "NA")
    merged, have_prime = add_prime(merged)
    cs = signalp_cleavage()
    if cs is not None:
        merged["in_signal_peptide"] = merged["start"] <= cs

    lines.append("Binding tools used: " + ", ".join(used))
    lines.append("PRIME included: %s; SignalP included: %s" % (have_prime, cs is not None))
    lines.append("Thresholds (class I %%rank): SB < %.1f, WB < %.1f" %
                 (c.THRESHOLDS["I"]["SB"], c.THRESHOLDS["I"]["WB"]))
    lines.append("")
    lines.append("Call counts per tool and allele:")
    for t in used:
        tab = merged.groupby(["allele", "call_" + t]).size().unstack(fill_value=0)
        lines.append("-- " + t)
        lines.append(tab.to_string())

    rank_cols = used + (["prime_rank"] if have_prime else [])
    if len(rank_cols) > 1:
        lines.append("")
        lines.append("Spearman correlation of %rank (lower = better):")
        lines.append(merged[rank_cols].corr(method="spearman").round(3).to_string())
    n_sb = sum((merged["call_" + t] == "SB").astype(int) for t in used)
    if len(used) > 1:
        lines.append("")
        lines.append("Number of binding tools calling SB per peptide x allele:")
        lines.append(n_sb.value_counts().sort_index().to_string())

    need = 2 if len(used) > 1 else 1
    cons = merged[n_sb >= need].copy()
    cons["n_tools_SB"] = n_sb[n_sb >= need]
    cons = cons.sort_values(["n_tools_SB", used[0]], ascending=[False, True])
    cons.to_csv(os.path.join(R, "class_I_consensus.csv"), index=False)
    merged.to_csv(os.path.join(R, "class_I_merged.csv"), index=False)
    lines.append("")
    lines.append("Consensus SB pairs (>= %d binding tools): %d" % (need, len(cons)))

    if have_prime:
        sb = cons
        lines.append("  of which PRIME %%rank < 1: %d; PRIME %%rank >= 1: %d" %
                     (int(sb["prime_top1pct"].sum()), int((~sb["prime_top1pct"]).sum())))
        lines.append("  (PRIME ranks immunogenic potential on top of presentation; it never classifies.)")
    if cs is not None:
        if cs == 0:
            lines.append("SignalP: no signal peptide predicted for this protein.")
        else:
            n_in = int(cons["in_signal_peptide"].sum())
            lines.append("SignalP: signal peptide ends at residue %d; %d consensus SB pairs start inside it." % (cs, n_in))
            lines.append("  Localization and antigenicity are separate questions (slide 1).")
    show = ["start", "length", "peptide", "allele"] + used
    if have_prime:
        show.append("prime_rank")
    if cs is not None:
        show.append("in_signal_peptide")
    lines.append(cons.head(15)[show].round(4).to_string(index=False))

    p2 = os.path.join(R, "netmhciipan_class_II.csv")
    if os.path.exists(p2):
        d2 = pd.read_csv(p2)
        lines.append("")
        lines.append("Class II (NetMHCIIpan, 15-mers; SB < 2, WB < 10):")
        lines.append(d2.groupby(["allele", "call"]).size().unstack(fill_value=0).to_string())

    text = "\n".join(lines)
    with open(os.path.join(R, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
