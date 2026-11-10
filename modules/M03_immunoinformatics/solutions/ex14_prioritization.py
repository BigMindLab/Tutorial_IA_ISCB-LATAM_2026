"""Exercise 14 - Slide 24: how do we prioritise candidates?

Question Q5 of the module. Out of everything above, what goes to the bench?

RUNS WITH NO LICENSED SOFTWARE.

The slide makes two claims and this exercise turns both into code:

  (a) pVACtools is a PIPELINE, not a model. It is filters chained over
      predictions that already exist. Here that chain is built by hand over
      fhaC, using the topology from exercise 11, binding from MHCflurry and
      the immunogenicity model from exercise 12.

  (b) TIDE answers a DIFFERENT question. It works at the level of a tumour,
      over gene expression, and there is no peptide and no HLA anywhere in
      its diagram. Section 14.5 shows why it cannot be applied to a pathogen
      antigen such as fhaC.

Steps
-----
  14.1  Generate every 9-mer of fhaC and annotate it with the curated UniProt
        topology, which is the same reference exercise 11 builds.
  14.2  Binding and presentation with MHCflurry 2.0 over a panel of HLAs.
  14.3  Immunogenicity from the exercise 12 model, and the TCR column that
        honestly stays empty after exercise 13.
  14.4  pVACtools style filters, giving a ranked shortlist.
  14.5  The contrast with TIDE.

Run exercise 12 first, because this one loads the model it saves.

Usage
-----
    python ex14_prioritization.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from encoding import (  # noqa: E402
    DATA,
    RESULTS,
    encode_hla_pca,
    encode_peptide_pca,
    kmers,
    load_paratopes,
    load_pca12,
    read_fasta,
)
from ex12_immunogenicity import DeepImmunoCNN  # noqa: E402

SEED = 2026
np.random.seed(SEED)

# Allele panel: broad population coverage, and supported both by MHCflurry and
# by the DeepImmuno paratope table.
HLA_PANEL = {
    "HLA-A*02:01": "HLA-A*0201",
    "HLA-A*01:01": "HLA-A*0101",
    "HLA-A*24:02": "HLA-A*2402",
    "HLA-B*07:02": "HLA-B*0702",
    "HLA-B*44:02": "HLA-B*4402",
}

# Thresholds from slide 21: class I, SB below 0.5 %rank, WB below 2 %rank.
SB_THRESHOLD = 0.5
WB_THRESHOLD = 2.0


def uniprot_regions(path: Path) -> list[dict]:
    """Signal peptide, transmembrane strands and domains from a UniProt entry."""
    with open(path, encoding="utf-8") as fh:
        entry = json.load(fh)
    wanted = {"Signal", "Transmembrane", "Domain", "Topological domain", "Chain"}
    regions = []
    for feat in entry.get("features", []):
        if feat.get("type") not in wanted:
            continue
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        regions.append({"type": feat["type"], "start": start, "end": end,
                        "description": feat.get("description", "")})
    return regions


def annotate_region(start: int, end: int, regions: list[dict]) -> str:
    """Label a peptide by the UniProt region it falls inside."""
    labels = []
    for r in regions:
        if start >= r["start"] and end <= r["end"]:
            if r["type"] == "Signal":
                return "signal_peptide"
            if r["type"] in ("Transmembrane", "Topological domain"):
                labels.append(r["description"] or r["type"])
            elif r["type"] == "Domain":
                labels.append(r["description"])
    return ";".join(dict.fromkeys(labels)) if labels else "mature"


# --------------------------------------------------------------------------
def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    buffer: list[str] = []

    def log(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        buffer.append(line)

    log("=" * 78)
    log("EXERCISE 14 - Prioritising candidates (slide 24)")
    log("=" * 78)
    log("")
    log("  NO TRAINING HAPPENS IN THIS EXERCISE.")
    log("  Every number below comes from YOUR protein, fhaC. The only thing")
    log("  borrowed from elsewhere is the model weights saved by exercise 12,")
    log("  which were trained there on IEDB peptides.")

    # ---------------------------------------------------------------- 14.1
    log("\n[14.1] [fhaC] Candidates and topology")
    header, seq = read_fasta(DATA / "fhaC.fasta")
    log("      %s" % header[:90])
    log("      %d residues" % len(seq))

    regions = uniprot_regions(DATA / "P35077.json")
    log("\n      UniProt annotation used as a filter (%d regions):" % len(regions))
    for r in regions[:12]:
        log("        %-20s %4d-%-4d %s" % (r["type"], r["start"], r["end"],
                                           r["description"][:45]))
    if len(regions) > 12:
        log("        ... and %d more" % (len(regions) - 12))

    cand = kmers(seq, 9)
    cand["region"] = [annotate_region(s, e, regions)
                      for s, e in zip(cand["start"], cand["end"])]
    log("\n      %d overlapping 9-mers (%d residues - 9 + 1 = %d)"
        % (len(cand), len(seq), len(seq) - 9 + 1))
    log("      Exercise 11 reports %d rows because it works residue by residue."
        % len(seq))
    log("      This one works peptide by peptide. Same protein, different unit.")
    log("      breakdown by region:")
    for region, n in cand["region"].value_counts().items():
        log("        %5d  %s" % (n, region))
    n_signal = int((cand.region == "signal_peptide").sum())
    log("\n      %d fall inside the signal peptide. That fragment is cut off" % n_signal)
    log("      during secretion, which is what SignalP answered in Q1 on slide 18,")
    log("      so those peptides are not part of the mature protein and are dropped.")

    # ---------------------------------------------------------------- 14.2
    log("\n[14.2] [fhaC] Binding and presentation (MHCflurry 2.0)")
    try:
        from mhcflurry import Class1PresentationPredictor
    except Exception as exc:  # noqa: BLE001
        log("      [ERROR] mhcflurry unavailable: %s" % exc)
        log("      pip install mhcflurry")
        log("      mhcflurry-downloads fetch models_class1_presentation")
        return 1

    predictor = Class1PresentationPredictor.load()
    peptides = cand["peptide"].tolist()
    chunks = []
    for allele in HLA_PANEL:
        log("      predicting %s over %d peptides ..." % (allele, len(peptides)))
        # predict()         gives presentation, which folds in antigen processing
        # predict_affinity() gives affinity alone, with its own percentile
        # They are two different questions, so both are requested.
        r = predictor.predict(peptides=peptides, alleles=[allele], verbose=0)
        a = predictor.predict_affinity(peptides=peptides,
                                       alleles={"sample": [allele]}, verbose=0)
        chunks.append(pd.DataFrame({
            "peptide": r["peptide"].to_numpy(),
            "HLA": allele,
            "affinity_nM": r["affinity"].to_numpy(),
            "affinity_rank": a["affinity_percentile"].to_numpy(),
            "presentation_rank": r["presentation_percentile"].to_numpy(),
            "presentation_score": r["presentation_score"].to_numpy(),
        }))
    binding = pd.concat(chunks, ignore_index=True)
    table = cand.merge(binding, on="peptide", how="right")
    log("      %d rows (peptide by allele)" % len(table))

    table["binding_call"] = np.where(
        table.affinity_rank < SB_THRESHOLD, "strong (SB)",
        np.where(table.affinity_rank < WB_THRESHOLD, "weak (WB)", "non binder"))
    log("\n      breakdown using the slide 21 thresholds:")
    for k, v in table["binding_call"].value_counts().items():
        log("        %6d  %s" % (v, k))

    # ---------------------------------------------------------------- 14.3
    log("\n[14.3] [fhaC] Immunogenicity (pretrained weights) and TCR evidence")
    model_path = DATA / "precomputed" / "ex12_immunogenicity_cnn.pt"
    if not model_path.exists():
        log("      [ERROR] missing %s" % model_path.name)
        log("      run first: python ex12_immunogenicity.py")
        return 1

    model = DeepImmunoCNN()
    model.load_state_dict(torch.load(model_path, weights_only=True))
    model.eval()
    pca12, paratopes = load_pca12(), load_paratopes()

    scores = np.zeros(len(table))
    for mhcflurry_allele, deepimmuno_allele in HLA_PANEL.items():
        mask = table.HLA == mhcflurry_allele
        sub = table[mask]
        Xp = np.stack([encode_peptide_pca(p, pca12)
                       for p in sub["peptide"]]).astype(np.float32)
        Xm = np.repeat(encode_hla_pca(deepimmuno_allele, paratopes, pca12)[None, ...],
                       len(sub), axis=0).astype(np.float32)
        with torch.no_grad():
            scores[mask.to_numpy()] = torch.sigmoid(
                model(torch.tensor(Xp), torch.tensor(Xm))).numpy()
    table["immunogenicity_score"] = scores
    log("      immunogenicity score for %d pairs" % len(table))

    table["TCR_evidence"] = "no evidence"
    log("      TCR evidence: 'no evidence' for all of them.")
    log("      That is not laziness. Exercise 13 showed NetTCR is only valid over")
    log("      its 17 training peptides, and fhaC is not among them. Putting a")
    log("      score there would be inventing information.")

    # ---------------------------------------------------------------- 14.4
    log("\n[14.4] [fhaC] The pVACtools style filter chain")
    steps = []
    t = table.copy()
    steps.append(("0. all peptide and allele pairs", len(t)))

    t = t[t.region != "signal_peptide"]
    steps.append(("1. outside the signal peptide (Q1, SignalP)", len(t)))

    t = t[t.affinity_rank < WB_THRESHOLD]
    steps.append(("2. binding, affinity %%rank below %.1f (Q2)" % WB_THRESHOLD, len(t)))

    t = t[t.presentation_rank < WB_THRESHOLD]
    steps.append(("3. presentation, %%rank below %.1f (Q2)" % WB_THRESHOLD, len(t)))

    t = t[t.immunogenicity_score > t.immunogenicity_score.median()]
    steps.append(("4. immunogenicity above the median (Q3)", len(t)))

    log("")
    for name, n in steps:
        log("        %6d   %s" % (n, name))
    log("\n      Each filter is a different question from the module. None of them")
    log("      proves anything. They only shorten the list by an explicit criterion.")

    # Combined score. Binding and presentation go through -log10, immunogenicity
    # is used as is. The weights are an analyst decision, not a model output.
    WEIGHTS = {"binding": 0.35, "presentation": 0.35, "immunogenicity": 0.30}
    t = t.copy()
    for col, dest in [("affinity_rank", "n_binding"),
                      ("presentation_rank", "n_presentation")]:
        v = -np.log10(np.clip(t[col].to_numpy(), 1e-3, None))
        t[dest] = (v - v.min()) / (v.max() - v.min() + 1e-9)
    v = t["immunogenicity_score"].to_numpy()
    t["n_immunogenicity"] = (v - v.min()) / (v.max() - v.min() + 1e-9)
    t["final_score"] = (WEIGHTS["binding"] * t.n_binding
                        + WEIGHTS["presentation"] * t.n_presentation
                        + WEIGHTS["immunogenicity"] * t.n_immunogenicity)

    shortlist = t.sort_values("final_score", ascending=False).reset_index(drop=True)
    cols = ["peptide", "start", "end", "region", "HLA", "affinity_nM",
            "affinity_rank", "presentation_rank", "immunogenicity_score",
            "TCR_evidence", "final_score"]
    log("\n      Weights used: %s" % WEIGHTS)
    log("\n      SHORTLIST, top 15:")
    log(shortlist[cols].head(15).round(3).to_string(index=False))

    shortlist[cols].to_csv(RESULTS / "ex14_shortlist.csv", index=False)
    table.to_csv(RESULTS / "ex14_full_table.csv", index=False)
    log("\n      results/ex14_shortlist.csv")
    log("      results/ex14_full_table.csv")

    log("\n      Allele coverage in the shortlist:")
    for allele, n in shortlist.HLA.value_counts().items():
        log("        %5d  %s" % (n, allele))
    log("      If an allele ends up with no candidates, that patient group has no")
    log("      coverage. That is a design criterion, not a detail.")

    # Figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    names = [s[0].split(". ", 1)[1][:28] for s in steps]
    ax[0].barh(range(len(steps)), [s[1] for s in steps], color="#5b2d8e")
    ax[0].set_yticks(range(len(steps)))
    ax[0].set_yticklabels(names, fontsize=8)
    ax[0].invert_yaxis()
    ax[0].set_xscale("log")
    ax[0].set_xlabel("candidates remaining (log)")
    ax[0].set_title("The filter funnel")

    ax[1].scatter(-np.log10(np.clip(table.affinity_rank, 1e-3, None)),
                  table.immunogenicity_score, s=6, alpha=0.3, c="#9aa5b1")
    ax[1].scatter(-np.log10(np.clip(shortlist.affinity_rank, 1e-3, None)),
                  shortlist.immunogenicity_score, s=14, c="#d1495b", label="shortlist")
    ax[1].set_xlabel("-log10(affinity %rank)")
    ax[1].set_ylabel("immunogenicity score")
    ax[1].set_title("The two axes are not correlated")
    ax[1].legend()

    top = shortlist.head(12)[::-1]
    ax[2].barh(range(len(top)), top.final_score,
               color=["#2a9d8f" if h == "HLA-A*02:01" else "#5b2d8e" for h in top.HLA])
    ax[2].set_yticks(range(len(top)))
    ax[2].set_yticklabels(["%s | %s" % (r.peptide, r.HLA.split("-")[1])
                           for _, r in top.iterrows()], fontsize=8)
    ax[2].set_xlabel("combined score")
    ax[2].set_title("Top 12")
    fig.suptitle("Exercise 14 - From thousands of peptides to a defensible shortlist",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(RESULTS / "ex14_prioritization.png", dpi=140)
    plt.close(fig)
    log("      figure written to results/ex14_prioritization.png")

    # ---------------------------------------------------------------- 14.5
    log("\n[14.5] And TIDE? The other half of the slide")
    log("      TIDE appears nowhere in any of the above, and that is the point.")
    log("      Its input is the expression profile of a TUMOUR and its output is")
    log("      ONE score per tumour, not per peptide. The slide says it outright:")
    log("      'No peptide, no HLA allele, no epitope anywhere in this diagram'.")
    log("")
    log("      The two questions side by side:")
    comparison = pd.DataFrame({
        "aspect": ["unit of analysis", "input", "output", "needs HLA?",
                   "needs the antigen?", "works for fhaC?"],
        "pVACtools (what we did)": [
            "peptide by allele", "protein plus HLA plus predictions",
            "ranked epitope list", "yes", "yes", "yes"],
        "TIDE": [
            "one tumour, one sample", "RNA-seq or NanoString expression",
            "a checkpoint response score", "no", "no",
            "no, fhaC is a pathogen antigen"],
    })
    log(comparison.to_string(index=False))
    comparison.to_csv(RESULTS / "ex14_pvactools_vs_tide.csv", index=False)
    log("\n      results/ex14_pvactools_vs_tide.csv")
    log("      Closing line from the slide: 'Pathogen epitopes and tumour")
    log("      neoantigens are different problems'. A neoantigen needs the")
    log("      patient's own HLA type and sequencing of THEIR tumour. fhaC does not.")

    log("\n" + "=" * 78)
    log("QUESTIONS TO DISCUSS")
    log("=" * 78)
    log("  1. We chose the weights of the combined score ourselves. Change them")
    log("     and watch the top 15 move. What does that say about a shortlist?")
    log("  2. We filtered on the median immunogenicity. Why a median rather than")
    log("     an absolute threshold? What is wrong with each option?")
    log("  3. The TCR_evidence column is empty. What experiment would fill it?")
    log("  4. We did not use conservation across B. pertussis strains. Why would")
    log("     that be the obvious next filter for a vaccine?")
    log("  5. fhaC reaches the immune system as an extracellular protein, so its")
    log("     natural route is MHC class II. What does that make this class I")
    log("     shortlist, and which exercise would you run next?")

    (RESULTS / "ex14_log.txt").write_text("\n".join(buffer), encoding="utf-8")
    print("\n[log] results/ex14_log.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
