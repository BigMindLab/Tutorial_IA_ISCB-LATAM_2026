"""Exercise 13 - Slide 23: does this TCR recognise this peptide?

Question Q4 of the module. NetTCR-2.0 is reimplemented here exactly as the
slide describes it (Montemurro et al., Commun Biol 2021, Fig. 10):

    INPUT    peptide plus CDR3 alpha and CDR3 beta
    ENCODE   BLOSUM50, 20 values per residue, zero padded to 30
    MODEL    1D convolutions on each sequence, 16 filters at kernel sizes
             1, 3, 5, 7 and 9. Max pooled into one 240 value vector, then a
             dense layer of 32, then one sigmoid output.

RUNS WITH NO LICENSED SOFTWARE.

Then the warning on the slide is put to the test:

    "The authors state the model can at present only be applied to the
     peptides in its training data."

Experiments
-----------
  13.1  Look at the data. 17 peptides, and how unbalanced they are.
  13.2  Train on partitions 1 to 4, evaluate on partition 5, which is the
        cross validation the repository itself ships.
  13.3  AUC per peptide. It works where there are thousands of TCRs and it
        fails where there are dozens.
  13.4  THE KEY EXPERIMENT: leave one peptide out. Train without GILGFVFTL
        and evaluate on GILGFVFTL. The AUC collapses. That is extrapolation.
  13.5  Try to score an fhaC peptide. A number comes out and it means
        nothing. The diagnostic at the end shows why.

Usage
-----
    python ../utils/download_data.py
    python ex13_nettcr.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from encoding import (  # noqa: E402
    DATA,
    RESULTS,
    encode_blosum,
    kmers,
    load_blosum50,
    read_fasta,
)

SEED = 2026
CDR3_WIDTH = 30   # the slide says zero padded to length 30
PEPTIDE_WIDTH = 9  # every peptide in this dataset is a 9-mer
torch.manual_seed(SEED)
np.random.seed(SEED)


# ==========================================================================
class NetTCR20(nn.Module):
    """NetTCR-2.0 exactly as slide 23 draws it.

    Three identical branches, one each for peptide, CDR3 alpha and CDR3 beta.
    In each branch, five parallel convolutions with kernel sizes 1, 3, 5, 7
    and 9, 16 filters each. Each convolution reduces to one maximum per filter:

        3 sequences x 5 kernels x 16 filters = 240 values

    Those 240 values go through a dense layer of 32 and then a sigmoid.
    """

    KERNELS = (1, 3, 5, 7, 9)
    FILTERS = 16

    def __init__(self) -> None:
        super().__init__()
        self.convs = nn.ModuleList(
            [
                nn.ModuleList(
                    [nn.Conv1d(20, self.FILTERS, k, padding=k // 2) for k in self.KERNELS]
                )
                for _ in range(3)  # peptide, CDR3a, CDR3b
            ]
        )
        self.dense = nn.Sequential(nn.Linear(240, 32), nn.ReLU(), nn.Dropout(0.3))
        self.out = nn.Linear(32, 1)

    def _branch(self, x: torch.Tensor, i: int) -> torch.Tensor:
        x = x.transpose(1, 2)  # (B, L, 20) becomes (B, 20, L)
        parts = [torch.sigmoid(conv(x)).max(dim=2).values for conv in self.convs[i]]
        return torch.cat(parts, dim=1)  # (B, 80)

    def forward(self, peptide, a, b):
        v = torch.cat([self._branch(peptide, 0), self._branch(a, 1), self._branch(b, 2)], dim=1)
        return self.out(self.dense(v)).squeeze(-1)


def encode_batch(df: pd.DataFrame, blosum) -> tuple[np.ndarray, ...]:
    peptide = np.stack([encode_blosum(p, blosum, PEPTIDE_WIDTH) for p in df["peptide"]])
    a = np.stack([encode_blosum(s, blosum, CDR3_WIDTH) for s in df["CDR3a"]])
    b = np.stack([encode_blosum(s, blosum, CDR3_WIDTH) for s in df["CDR3b"]])
    return peptide, a, b, df["binder"].to_numpy().astype(np.float32)


def train(tr, te, blosum, epochs=25, log=print, tag=""):
    Xp, Xa, Xb, y = encode_batch(tr, blosum)
    Vp, Va, Vb, yv = encode_batch(te, blosum)

    model = NetTCR20()
    opt = torch.optim.Adam(model.parameters(), lr=5e-4)
    pos_w = torch.tensor([(y == 0).sum() / max((y == 1).sum(), 1)], dtype=torch.float32)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_w)

    tp, ta, tb, ty = (torch.tensor(z) for z in (Xp, Xa, Xb, y))
    vp, va, vb = (torch.tensor(z) for z in (Vp, Va, Vb))
    n, batch = len(ty), 128

    for ep in range(epochs):
        model.train()
        order = torch.randperm(n)
        for i in range(0, n, batch):
            k = order[i : i + batch]
            opt.zero_grad()
            criterion(model(tp[k], ta[k], tb[k]), ty[k]).backward()
            opt.step()
        if (ep + 1) % 5 == 0:
            model.eval()
            with torch.no_grad():
                s = torch.sigmoid(model(vp, va, vb)).numpy()
            try:
                log("      %sepoch %2d/%d  AUC=%.3f" % (tag, ep + 1, epochs, roc_auc_score(yv, s)))
            except ValueError:
                log("      %sepoch %2d/%d  (only one class present)" % (tag, ep + 1, epochs))

    model.eval()
    with torch.no_grad():
        s = torch.sigmoid(model(vp, va, vb)).numpy()
    return model, s, yv


def score(model, df, blosum) -> np.ndarray:
    Xp, Xa, Xb, _ = encode_batch(df.assign(binder=0), blosum)
    with torch.no_grad():
        return torch.sigmoid(
            model(torch.tensor(Xp), torch.tensor(Xa), torch.tensor(Xb))
        ).numpy()


# ==========================================================================
def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    buffer: list[str] = []

    def log(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        buffer.append(line)

    log("=" * 78)
    log("EXERCISE 13 - TCR recognition (slide 23)")
    log("=" * 78)
    log("")
    log("  TWO DATASETS, DO NOT CONFUSE THEM")
    log("  ---------------------------------")
    log("  [TRAINING]  16,464 TCR and peptide pairs published with NetTCR-2.0,")
    log("              covering only 17 distinct peptides. fhaC cannot be used")
    log("              for training because no TCR has ever been measured")
    log("              against it. Used in blocks 13.1 to 13.4.")
    log("")
    log("  [fhaC]      YOUR protein. It appears only in block 13.5, and the")
    log("              conclusion there is NEGATIVE: this model must not be")
    log("              used on fhaC. That negative answer is the result.")

    blosum = load_blosum50()
    log("\n[13.0] BLOSUM50 encoding")
    log("      matrix %s, divided by 5 to keep the range comfortable" % (blosum.shape,))
    example = encode_blosum("CASSLAPGATNEKLFF", blosum, CDR3_WIDTH)
    log("      CASSLAPGATNEKLFF becomes %s, %d rows carry signal and the rest is padding"
        % (example.shape, (example != 0).any(axis=1).sum()))

    # ---------------------------------------------------------------- 13.1
    log("\n[13.1] [TRAINING] The NetTCR-2.0 data")
    df = pd.read_csv(DATA / "train_ab_90_alphabeta.csv")
    log("      %d pairs, binders=%d, non binders=%d"
        % (len(df), int(df.binder.sum()), int((df.binder == 0).sum())))
    log("      partitions: %s" % sorted(df.partition.unique()))

    summary = (df.groupby("peptide")["binder"]
                 .agg(positives="sum", total="count")
                 .sort_values("positives", ascending=False))
    log("\n      ONLY %d distinct peptides in the whole dataset:" % df.peptide.nunique())
    log(summary.to_string())
    log("\n      That number is the key to the entire slide. Three peptides hold")
    log("      most of the data and nine have fewer than 20 positive TCRs.")

    # ---------------------------------------------------------------- 13.2
    log("\n[13.2] [TRAINING] Train on partitions 1 to 4, evaluate on partition 5")
    tr = df[df.partition != 5].reset_index(drop=True)
    te = df[df.partition == 5].reset_index(drop=True)
    log("      train=%d  test=%d" % (len(tr), len(te)))
    model, s, yv = train(tr, te, blosum, log=log, tag="")
    auc_global = roc_auc_score(yv, s)
    log("      global AUC on partition 5: %.3f" % auc_global)

    # ---------------------------------------------------------------- 13.3
    log("\n[13.3] [TRAINING] AUC per peptide, which the global average hides")
    te = te.copy()
    te["score"] = s
    rows = []
    for peptide, g in te.groupby("peptide"):
        if g.binder.nunique() < 2:
            rows.append({"peptide": peptide, "n_test": len(g),
                         "positives": int(g.binder.sum()), "AUC": np.nan})
            continue
        rows.append({"peptide": peptide, "n_test": len(g),
                     "positives": int(g.binder.sum()),
                     "AUC": roc_auc_score(g.binder, g.score)})
    per_peptide = pd.DataFrame(rows).sort_values("n_test", ascending=False)
    per_peptide["n_train"] = per_peptide["peptide"].map(
        tr.groupby("peptide").size()).fillna(0).astype(int)
    log(per_peptide.to_string(index=False))
    per_peptide.to_csv(RESULTS / "ex13_auc_per_peptide.csv", index=False)

    valid = per_peptide.dropna(subset=["AUC"])
    if len(valid) > 2:
        r = np.corrcoef(np.log10(valid.n_train + 1), valid.AUC)[0, 1]
        log("\n      correlation between log10(training TCRs) and AUC: r=%.2f" % r)

    # ---------------------------------------------------------------- 13.4
    log("\n[13.4] [TRAINING] THE KEY EXPERIMENT: leave one peptide out")
    target = "GILGFVFTL"   # influenza M1, the peptide with the most data
    log("      Train with no TCR against %s, then evaluate on %s." % (target, target))
    log("      That is exactly the situation of a new peptide, like any from fhaC.")

    tr_without = df[df.peptide != target].reset_index(drop=True)
    te_only = df[(df.peptide == target) & (df.partition == 5)].reset_index(drop=True)
    log("      train=%d (without %s)  test=%d" % (len(tr_without), target, len(te_only)))
    _, s_without, y_without = train(tr_without, te_only, blosum, log=log, tag="[LOPO] ")
    auc_without = roc_auc_score(y_without, s_without)

    base = te[te.peptide == target]
    auc_with = roc_auc_score(base.binder, base.score)
    log("\n      AUC on %s WITH the peptide in training    : %.3f" % (target, auc_with))
    log("      AUC on %s WITHOUT the peptide in training : %.3f" % (target, auc_without))
    log("      (0.500 is chance)")
    log("\n      This is what the slide means when it says the model can only be")
    log("      applied to the peptides in its training data. It is not a bad")
    log("      model. It learned CDR3 motifs specific to a handful of peptides,")
    log("      not a general rule of recognition.")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    v2 = per_peptide.dropna(subset=["AUC"])
    ax[0].scatter(v2.n_train, v2.AUC, s=60, color="#5b2d8e")
    for _, row in v2.iterrows():
        ax[0].annotate(row.peptide, (row.n_train, row.AUC),
                       fontsize=7, xytext=(4, 3), textcoords="offset points")
    ax[0].set_xscale("log")
    ax[0].axhline(0.5, ls="--", c="gray")
    ax[0].set_xlabel("TCR peptide pairs in training (log)")
    ax[0].set_ylabel("AUC on partition 5")
    ax[0].set_title("How much data each peptide needs")

    ax[1].bar(["with the peptide\nin training", "without it\n(extrapolation)"],
              [auc_with, auc_without], color=["#2a9d8f", "#d1495b"])
    ax[1].axhline(0.5, ls="--", c="gray")
    ax[1].set_ylim(0, 1)
    ax[1].set_ylabel("AUC")
    ax[1].set_title("Leave one peptide out: %s" % target)
    for i, val in enumerate([auc_with, auc_without]):
        ax[1].text(i, val + 0.02, "%.3f" % val, ha="center", fontweight="bold")
    fig.suptitle("Exercise 13 - NetTCR-2.0 and the limits of its domain", fontsize=12)
    fig.tight_layout()
    fig.savefig(RESULTS / "ex13_limits.png", dpi=140)
    plt.close(fig)
    log("      figure written to results/ex13_limits.png")

    # ---------------------------------------------------------------- 13.5
    log("\n[13.5] [fhaC] What if we use it on an fhaC peptide?")
    _, seq = read_fasta(DATA / "fhaC.fasta")
    fhac_peptides = kmers(seq, 9)["peptide"].iloc[:5].tolist()
    demo_tcrs = df[df.binder == 1].iloc[:3][["CDR3a", "CDR3b"]]
    probe = pd.DataFrame(
        [{"peptide": p, "CDR3a": r.CDR3a, "CDR3b": r.CDR3b}
         for p in fhac_peptides for _, r in demo_tcrs.iterrows()])
    probe["score"] = score(model, probe, blosum)
    log("      Scoring 5 fhaC peptides against 3 real TCRs from the dataset:")
    log(probe.to_string(index=False))
    log("\n      score range: %.3f to %.3f" % (probe.score.min(), probe.score.max()))

    # Diagnostic: what is the score actually responding to?
    sd_across_peptides = probe.groupby(["CDR3a", "CDR3b"])["score"].std().mean()
    sd_across_tcrs = probe.groupby("peptide")["score"].std().mean()
    log("")
    log("      LOOK AT THIS CAREFULLY. Standard deviation of the score ...")
    log("        ... when the PEPTIDE changes, TCR fixed : %.4f" % sd_across_peptides)
    log("        ... when the TCR changes, peptide fixed : %.4f" % sd_across_tcrs)
    if sd_across_tcrs > 3 * sd_across_peptides:
        log("")
        log("      The score barely moves when the peptide changes but it does move")
        log("      when the TCR changes. On peptides it has never seen, the model is")
        log("      essentially scoring the CDR3 and ignoring the peptide. That is the")
        log("      clearest possible demonstration that it is not answering the")
        log("      question we are asking it.")

    log("")
    log("      The model RETURNS a number for every pair. None of them is")
    log("      interpretable: fhaC is not in the training set, these TCRs were")
    log("      never measured against it, and there is no HLA context in the input.")
    log("      Practical consequence for exercise 14: the TCR column for fhaC is")
    log("      marked 'no evidence' rather than filled with a score.")
    probe.to_csv(RESULTS / "ex13_fhaC_no_evidence.csv", index=False)

    log("\n" + "=" * 78)
    log("QUESTIONS TO DISCUSS")
    log("=" * 78)
    log("  1. The global AUC on partition 5 is high. Looking at the per peptide")
    log("     table, why is reporting that number alone misleading?")
    log("  2. Non binders are built by pairing TCRs with random peptides. What")
    log("     bias does that introduce, and how would you notice it?")
    log("  3. The input does not include the HLA. Two people with different HLA")
    log("     and the same CDR3 get the same score. Is that right?")
    log("  4. If you wanted to use NetTCR for fhaC, what experimental data would")
    log("     you have to generate first?")

    (RESULTS / "ex13_log.txt").write_text("\n".join(buffer), encoding="utf-8")
    print("\n[log] results/ex13_log.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
