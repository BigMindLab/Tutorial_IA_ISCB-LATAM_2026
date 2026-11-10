"""Exercise 15 - NetMHCpan 4.1 and NetMHCIIpan 4.0 without a local licence.

RUNS WITH NO LICENSED SOFTWARE, but it does need internet.

Why this exercise exists
------------------------
Slides 20 and 21 are built on NetMHCpan and NetMHCIIpan. Both need an academic
licence from DTU and a manual install, which stops a workshop before it starts.

IEDB runs the same tools behind a public web API. This exercise sends the fhaC
sequence there and parses the result, so everyone gets real NetMHCpan numbers
on the case study even with nothing installed.

What this is and is not
-----------------------
It IS the real NetMHCpan and NetMHCIIpan, run on IEDB's cluster.
It is NOT a substitute for a local install when you need thousands of
sequences, offline work, or the binding affinity output. The API here returns
the eluted ligand (EL) output, which is the one to rank on anyway.

Be considerate. The API is a free shared service. This script sends one
request per allele, waits between calls, and should not be put in a loop.

Usage
-----
    python3 ex15_netmhcpan_iedb.py
    python3 ex15_netmhcpan_iedb.py --class-i-only

Outputs
-------
    results/netmhcpan_class_I.csv      same columns as the local wrapper
    results/netmhciipan_class_II.csv   same columns as the local wrapper
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request

_sys_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_sys_dir, "..", "utils"))

import common as c  # noqa: E402

CLASS_I_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
CLASS_II_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
PAUSE_SECONDS = 3


def post(url: str, fields: dict) -> str:
    data = urllib.parse.urlencode(fields).encode()
    # Send no Accept header. The IEDB endpoint answers 406 Not Acceptable
    # when one is present, even for text/plain, which is what it returns.
    request = urllib.request.Request(
        url, data=data, headers={"User-Agent": "iscb-latam-tutorial"})
    with urllib.request.urlopen(request, timeout=300) as response:
        return response.read().decode("utf-8", "replace")


def parse_table(text: str) -> tuple[list[str], list[list[str]]]:
    """IEDB answers with a tab separated table and one header line."""
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:] if "\t" in l]
    return header, [r for r in rows if len(r) == len(header)]


def column(header: list[str], *names: str):
    for name in names:
        if name in header:
            return header.index(name)
    return None


def run_class_i(sequence: str, log) -> list[dict]:
    rows_out: list[dict] = []
    lengths = ",".join(str(k) for k in c.CLASS_I_LENGTHS)
    for allele in c.CLASS_I_ALLELES:
        log("      NetMHCpan 4.1  %s  lengths %s" % (allele, lengths))
        text = post(CLASS_I_URL, {
            "method": "netmhcpan_el",
            "sequence_text": sequence,
            "allele": ",".join([allele] * len(c.CLASS_I_LENGTHS)),
            "length": lengths,
        })
        header, rows = parse_table(text)
        if not rows:
            log("      [warning] empty answer for %s: %s" % (allele, text[:160]))
            continue
        i_start = column(header, "start")
        i_len = column(header, "length")
        i_pep = column(header, "peptide")
        i_score = column(header, "score")
        i_rank = column(header, "percentile_rank", "rank")
        for r in rows:
            try:
                rank = float(r[i_rank])
            except (TypeError, ValueError):
                continue
            rows_out.append({
                "start": int(r[i_start]),
                "length": int(r[i_len]),
                "peptide": r[i_pep],
                "allele": allele,
                "score_EL": float(r[i_score]),
                "rank_EL": rank,
                # The API returns the eluted ligand output only, so the two
                # binding affinity columns of the local wrapper stay empty
                # rather than being filled with something they are not.
                "score_BA": "",
                "rank_BA": "",
                "affinity_nM": "",
                "call": c.classify(rank, "I"),
            })
        log("        %d peptides" % len([x for x in rows_out if x["allele"] == allele]))
        time.sleep(PAUSE_SECONDS)
    return rows_out


def run_class_ii(sequence: str, log) -> list[dict]:
    rows_out: list[dict] = []
    for allele in c.CLASS_II_ALLELES:
        # common.py keeps the NetMHCIIpan spelling, DRB1_0401. The API wants
        # DRB1*04:01, so convert rather than making people edit two lists.
        api_allele = allele
        if "_" in allele:
            family, digits = allele.split("_", 1)
            if len(digits) == 4:
                api_allele = "%s*%s:%s" % (family, digits[:2], digits[2:])
        log("      NetMHCIIpan 4.0  %s  sent as %s" % (allele, api_allele))
        text = post(CLASS_II_URL, {
            "method": "netmhciipan_el",
            "sequence_text": sequence,
            "allele": api_allele,
            "length": str(c.CLASS_II_LENGTH),
        })
        header, rows = parse_table(text)
        if not rows:
            log("      [warning] empty answer for %s: %s" % (allele, text[:160]))
            continue
        i_start = column(header, "start")
        i_pep = column(header, "peptide")
        i_score = column(header, "score")
        i_rank = column(header, "rank", "percentile_rank")
        for r in rows:
            try:
                rank = float(r[i_rank])
            except (TypeError, ValueError):
                continue
            rows_out.append({
                "start": int(r[i_start]),
                "peptide": r[i_pep],
                "allele": allele,
                "score_EL": float(r[i_score]),
                "rank_EL": rank,
                "call": c.classify(rank, "II"),
            })
        log("        %d peptides" % len([x for x in rows_out if x["allele"] == allele]))
        time.sleep(PAUSE_SECONDS)
    return rows_out


def summarise(rows, label, log):
    if not rows:
        log("      no rows for %s" % label)
        return
    counts = {}
    for r in rows:
        counts[r["call"]] = counts.get(r["call"], 0) + 1
    log("      %s: %d rows, %s" % (label, len(rows), counts))
    best = sorted(rows, key=lambda r: r["rank_EL"])[:5]
    log("      top 5 by %rank:")
    for r in best:
        log("        %-16s %-12s rank_EL=%-8.4g %s"
            % (r["peptide"], r["allele"], r["rank_EL"], r["call"]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", default=c.DEFAULT_FASTA)
    ap.add_argument("--class-i-only", action="store_true")
    ap.add_argument("--class-ii-only", action="store_true")
    args = ap.parse_args()

    buffer: list[str] = []

    def log(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        buffer.append(line)

    log("=" * 78)
    log("EXERCISE 15 - NetMHCpan and NetMHCIIpan through the IEDB API")
    log("=" * 78)
    log("")
    log("  These are the real DTU tools, run on IEDB's public cluster, so no")
    log("  local licence is needed. The API returns the eluted ligand output,")
    log("  which is the column to rank on anyway.")
    log("")
    log("  Thresholds, from slide 21, and they differ by class:")
    log("    class I   SB below %.1f   WB below %.1f" % (c.THRESHOLDS["I"]["SB"],
                                                        c.THRESHOLDS["I"]["WB"]))
    log("    class II  SB below %.1f   WB below %.1f" % (c.THRESHOLDS["II"]["SB"],
                                                         c.THRESHOLDS["II"]["WB"]))

    sequence = list(c.read_fasta(args.fasta).values())[0]
    log("\n      protein: %d residues from %s"
        % (len(sequence), os.path.basename(args.fasta)))

    if not args.class_ii_only:
        log("\n[15.1] Class I, presented to CD8 cytotoxic T cells")
        try:
            rows = run_class_i(sequence, log)
        except Exception as exc:  # noqa: BLE001
            log("      [ERROR] %s" % exc)
            log("      The IEDB API needs internet. Behind a proxy or offline,")
            log("      use the local wrapper ex1_netmhcpan.py instead.")
            rows = []
        if rows:
            fields = ["start", "length", "peptide", "allele", "score_EL", "rank_EL",
                      "score_BA", "rank_BA", "affinity_nM", "call"]
            c.write_csv(os.path.join(c.RESULTS, "netmhcpan_class_I.csv"), rows, fields)
            summarise(rows, "class I", log)

    if not args.class_i_only:
        log("\n[15.2] Class II, presented to CD4 helper T cells")
        log("      This is the one that matters most for fhaC. An extracellular")
        log("      bacterial protein is endocytosed and presented on class II,")
        log("      so class I results are a hypothesis about cross presentation.")
        try:
            rows2 = run_class_ii(sequence, log)
        except Exception as exc:  # noqa: BLE001
            log("      [ERROR] %s" % exc)
            rows2 = []
        if rows2:
            fields2 = ["start", "peptide", "allele", "score_EL", "rank_EL", "call"]
            c.write_csv(os.path.join(c.RESULTS, "netmhciipan_class_II.csv"), rows2, fields2)
            summarise(rows2, "class II", log)

    log("\n" + "=" * 78)
    log("QUESTIONS TO DISCUSS")
    log("=" * 78)
    log("  1. The API gives the eluted ligand output and not binding affinity.")
    log("     Slide 21 says BA is trained on measured affinities and EL on")
    log("     peptides eluted from cells. Which one would you rank on, and why?")
    log("  2. Compare these class I hits against MHCflurry and MixMHCpred. How")
    log("     many peptides do all three call a strong binder?")
    log("  3. Class II uses a threshold of 2 where class I uses 0.5. Why can a")
    log("     looser threshold still be the right call for class II?")
    log("  4. A class II 15-mer is scored through a 9-mer core. Look at the core")
    log("     column: do several 15-mers share one core? What does that mean for")
    log("     counting candidates?")

    with open(os.path.join(c.RESULTS, "ex15_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(buffer))
    print("\n[log] results/ex15_log.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
