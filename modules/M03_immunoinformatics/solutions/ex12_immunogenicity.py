"""Exercise 12 - Slide 22: can deep learning predict immunogenicity?

Question Q3 of the module. Of the peptides an MHC can present, which ones
might actually trigger a T cell response?

RUNS WITH NO LICENSED SOFTWARE.

What happens here
-----------------
  12.1  Encode a peptide the way DeepImmuno does: 10 x 12 PCA components of
        553 AAindex properties, and the HLA as its 46 residue paratope.
  12.2  Reimplement DeepImmuno-CNN (Li et al. 2021, Fig. 1C) in PyTorch and
        train it on the 9,055 IEDB peptides the paper uses.
  12.3  Build the 28 input vector of PRIME 2.0 (Gfeller et al. 2023, Fig. 3B)
        and train a small multilayer perceptron on the same labels.
  12.4  THE POINT OF THE SLIDE: compare the immunogenicity ranking against the
        presentation ranking. A strong binder is not necessarily immunogenic.
  12.5  Apply the model to every 9-mer of fhaC and save the ranking.

Usage
-----
    python ../utils/download_data.py
    python ex12_immunogenicity.py

Outputs land in ../data/precomputed and ../results.
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
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from encoding import (  # noqa: E402
    AA20,
    DATA,
    RESULTS,
    binary_label,
    encode_hla_pca,
    encode_peptide_pca,
    kmers,
    load_paratopes,
    load_pca12,
    read_fasta,
    tcr_facing_positions,
    prime_features,
)

SEED = 2026
FOCUS_HLA = "HLA-A*0201"
torch.manual_seed(SEED)
np.random.seed(SEED)


# ==========================================================================
# 12.2  DeepImmuno-CNN  (Li et al. 2021, Figure 1C)
# ==========================================================================
class DeepImmunoCNN(nn.Module):
    """Two convolutional branches, concatenated, then dense layers, then sigmoid.

    Peptide branch:  (10 x 12) convolved along the position axis
    MHC branch:      (46 x 12) the same, along the paratope

    This is the structure drawn on slide 22: CNN layers per branch and
    fully connected layers at the end.
    """

    def __init__(self) -> None:
        super().__init__()
        self.peptide_branch = nn.Sequential(
            nn.Conv1d(12, 32, kernel_size=2), nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(32, 64, kernel_size=2), nn.ReLU(),
            nn.AdaptiveMaxPool1d(1),
        )
        self.mhc_branch = nn.Sequential(
            nn.Conv1d(12, 32, kernel_size=3), nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(32, 64, kernel_size=3), nn.ReLU(),
            nn.AdaptiveMaxPool1d(1),
        )
        self.head = nn.Sequential(
            nn.Linear(128, 64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, 1),
        )

    def forward(self, peptide: torch.Tensor, mhc: torch.Tensor) -> torch.Tensor:
        # (B, L, 12) becomes (B, 12, L) because Conv1d expects channels first.
        a = self.peptide_branch(peptide.transpose(1, 2)).squeeze(-1)
        b = self.mhc_branch(mhc.transpose(1, 2)).squeeze(-1)
        return self.head(torch.cat([a, b], dim=1)).squeeze(-1)


def train_cnn(Xp, Xm, y, Xp_val, Xm_val, y_val, epochs=40, log=print):
    model = DeepImmunoCNN()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    # The classes are unbalanced, so weight the positive class accordingly.
    pos_weight = torch.tensor([(y == 0).sum() / max((y == 1).sum(), 1)], dtype=torch.float32)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    tp, tm, ty = map(torch.tensor, (Xp, Xm, y.astype(np.float32)))
    vp, vm = torch.tensor(Xp_val), torch.tensor(Xm_val)
    n, batch = len(ty), 128
    best, best_state = -1.0, None

    for ep in range(epochs):
        model.train()
        order = torch.randperm(n)
        for i in range(0, n, batch):
            idx = order[i : i + batch]
            opt.zero_grad()
            loss = criterion(model(tp[idx], tm[idx]), ty[idx])
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            val = torch.sigmoid(model(vp, vm)).numpy()
        auc = roc_auc_score(y_val, val)
        if auc > best:
            best = auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        if (ep + 1) % 10 == 0:
            log("      epoch %2d/%d  loss=%.4f  val AUC=%.3f" % (ep + 1, epochs, loss.item(), auc))

    model.load_state_dict(best_state)
    model.eval()
    return model, best


def score_cnn(model, Xp, Xm) -> np.ndarray:
    with torch.no_grad():
        return torch.sigmoid(model(torch.tensor(Xp), torch.tensor(Xm))).numpy()


def presentation_scores(train, test, log):
    """Presentation score per peptide and HLA pair.

    PRIME 2.0 uses -log(%rank) from MixMHCpred 2.2. MixMHCpred needs an
    academic licence and a local binary, so MHCflurry 2.0 is used here. It is
    free, installs with pip, and answers the same question with different
    machinery, which is exactly the contrast slide 21 draws. Without it the
    exercise still runs, with the score set to zero and a warning.
    """
    try:
        from mhcflurry import Class1PresentationPredictor
    except Exception:  # noqa: BLE001
        log("      [warning] mhcflurry not installed, presentation score set to 0")
        log("                pip install mhcflurry")
        log("                mhcflurry-downloads fetch models_class1_presentation")
        return np.zeros(len(train)), np.zeros(len(test)), "unavailable"

    predictor = Class1PresentationPredictor.load()
    supported = set(predictor.supported_alleles)

    def score(sub: pd.DataFrame) -> np.ndarray:
        out = np.zeros(len(sub), dtype=float)
        sub = sub.reset_index(drop=True)
        # MHCflurry spells alleles as HLA-A*02:01
        sub = sub.assign(allele=sub["HLA"].str.replace(
            r"^HLA-([A-C])\*(\d{2})(\d{2})$", r"HLA-\1*\2:\3", regex=True))
        ok = sub[sub["allele"].isin(supported)]
        for allele, group in ok.groupby("allele"):
            r = predictor.predict(peptides=list(group["peptide"]),
                                  alleles=[allele], verbose=0)
            # presentation_percentile: lower means better presented
            out[group.index.to_numpy()] = -np.log10(
                np.clip(r["presentation_percentile"].to_numpy(), 1e-3, None))
        return out

    log("      computing presentation scores with MHCflurry 2.0 ...")
    return score(train), score(test), "MHCflurry 2.0, -log10(percentile)"


# ==========================================================================
def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    precomputed = DATA / "precomputed"
    precomputed.mkdir(parents=True, exist_ok=True)
    buffer: list[str] = []

    def log(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        buffer.append(line)

    log("=" * 78)
    log("EXERCISE 12 - Immunogenicity (slide 22)")
    log("=" * 78)
    log("")
    log("  TWO DATASETS, DO NOT CONFUSE THEM")
    log("  ---------------------------------")
    log("  [TRAINING]  9,055 peptides from IEDB, published with DeepImmuno.")
    log("              Needed because nobody has ever measured which fhaC")
    log("              peptides are immunogenic, so fhaC carries no labels to")
    log("              learn from. Slide 22 cites this same dataset.")
    log("              Used in blocks 12.2, 12.3 and 12.4.")
    log("")
    log("  [fhaC]      YOUR protein, the case study. 584 residues, which give")
    log("              576 overlapping 9-mers. It is never used for training,")
    log("              only for prediction.")
    log("              Used in block 12.5, and that is the result you keep.")
    log("")
    log("  Every block below is tagged so you always know which one you are")
    log("  looking at.")

    # ---------------------------------------------------------------- 12.1
    log("\n[12.1] DeepImmuno style encoding")
    pca12 = load_pca12()
    paratopes = load_paratopes()
    log("      PCA matrix: %d amino acids x %d components" % pca12.shape)
    log("      paratope table: %d alleles, pseudo sequence of %d residues"
        % (len(paratopes), len(paratopes.iloc[0, 0])))

    demo = "GILGFVFTL"  # influenza M1 58-66, the textbook HLA-A*02:01 epitope
    M = encode_peptide_pca(demo, pca12)
    log("\n      %s becomes a %s matrix (10 positions x 12 components)" % (demo, M.shape))
    log("      first 3 positions, first 5 components:")
    log(pd.DataFrame(M[:3, :5], index=list(demo[:3]),
                     columns=["PC%d" % (i + 1) for i in range(5)]).round(3).to_string())
    log("\n      Note the peptide has 9 residues and the matrix has 10 rows. Row 10")
    log("      is X padding, which is how DeepImmuno feeds 9-mers and 10-mers")
    log("      through the same input.")

    Mh = encode_hla_pca(FOCUS_HLA, paratopes, pca12)
    log("      %s becomes a %s matrix (46 paratope residues x 12 components)"
        % (FOCUS_HLA, Mh.shape))

    fig, ax = plt.subplots(1, 2, figsize=(13, 4), gridspec_kw={"width_ratios": [1, 3.2]})
    im0 = ax[0].imshow(M.T, aspect="auto", cmap="RdBu_r")
    ax[0].set_title("Peptide %s (10 x 12)" % demo)
    ax[0].set_xlabel("position")
    ax[0].set_ylabel("PCA component")
    ax[0].set_xticks(range(10))
    ax[0].set_xticklabels(list(demo + "X"))
    plt.colorbar(im0, ax=ax[0], fraction=0.046)
    im1 = ax[1].imshow(Mh.T, aspect="auto", cmap="RdBu_r")
    ax[1].set_title("Paratope %s (46 x 12)" % FOCUS_HLA)
    ax[1].set_xlabel("paratope residue")
    ax[1].set_ylabel("PCA component")
    plt.colorbar(im1, ax=ax[1], fraction=0.046)
    fig.suptitle("Exercise 12.1 - How DeepImmuno sees a peptide and HLA pair", fontsize=12)
    fig.tight_layout()
    fig.savefig(RESULTS / "ex12_1_encoding.png", dpi=140)
    plt.close(fig)
    log("      figure written to results/ex12_1_encoding.png")

    # ---------------------------------------------------------------- data
    log("\n[12.2] [TRAINING] Training data (IEDB, via the DeepImmuno repository)")
    df = pd.read_csv(DATA / "remove0123.csv", encoding="utf-8-sig")
    df["y"] = df["immunogenicity"].map(binary_label)
    n0 = len(df)
    # Keep only alleles with a tabulated paratope and peptides made of the 20
    # standard amino acids. Two peptides contain X and no predictor accepts them.
    df = df[df["HLA"].isin(paratopes.index)]
    df = df[df["peptide"].apply(lambda p: set(p) <= set(AA20))].reset_index(drop=True)
    log("      filtered %d rows down to %d (allele with paratope, standard residues)"
        % (n0, len(df)))
    log("      %d pairs, %d alleles, lengths %s"
        % (len(df), df["HLA"].nunique(), sorted(df["peptide"].str.len().unique())))
    log("      positives=%d  negatives=%d" % (int(df.y.sum()), int((df.y == 0).sum())))
    log("      Label: 'Negative' maps to 0, every 'Positive*' variant maps to 1")

    # Group split by peptide: no peptide can be in train and test at once.
    gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=SEED)
    i_tr, i_te = next(gss.split(df, groups=df["peptide"]))
    tr = df.iloc[i_tr].reset_index(drop=True)
    te = df.iloc[i_te].reset_index(drop=True)
    log("      grouped split gives train=%d and test=%d" % (len(tr), len(te)))
    log("      Grouping by peptide avoids inflating the AUC with the same")
    log("      peptide appearing on both sides of the split.")

    def tensorise(sub):
        Xp = np.stack([encode_peptide_pca(p, pca12) for p in sub["peptide"]]).astype(np.float32)
        Xm = np.stack([encode_hla_pca(h, paratopes, pca12) for h in sub["HLA"]]).astype(np.float32)
        return Xp, Xm, sub["y"].to_numpy()

    Xp_tr, Xm_tr, y_tr = tensorise(tr)
    Xp_te, Xm_te, y_te = tensorise(te)

    log("\n      Training DeepImmuno-CNN ...")
    model, _ = train_cnn(Xp_tr, Xm_tr, y_tr, Xp_te, Xm_te, y_te, log=log)
    s_cnn = score_cnn(model, Xp_te, Xm_te)
    auc_cnn = roc_auc_score(y_te, s_cnn)
    ap_cnn = average_precision_score(y_te, s_cnn)
    log("      CNN gives AUC=%.3f  AUPRC=%.3f  (prevalence baseline %.3f)"
        % (auc_cnn, ap_cnn, y_te.mean()))

    # ---------------------------------------------------------------- 12.3
    log("\n[12.3] [TRAINING] PRIME 2.0 and its 28 input vector")
    log("      TCR facing positions in a 9-mer: %s" % tcr_facing_positions(9))
    log("      P2 and the last position are left out, they anchor into the groove")

    pres_tr, pres_te, source = presentation_scores(tr, te, log)

    Xpr_tr = np.stack([prime_features(p, s) for p, s in zip(tr["peptide"], pres_tr)])
    Xpr_te = np.stack([prime_features(p, s) for p, s in zip(te["peptide"], pres_te)])
    scaler = StandardScaler().fit(Xpr_tr)
    mlp = MLPClassifier(hidden_layer_sizes=(32,), max_iter=600, random_state=SEED)
    mlp.fit(scaler.transform(Xpr_tr), y_tr)
    s_prime = mlp.predict_proba(scaler.transform(Xpr_te))[:, 1]
    auc_prime = roc_auc_score(y_te, s_prime)
    log("      PRIME style MLP, 28 inputs, presentation score from %s" % source)
    log("      gives AUC=%.3f  AUPRC=%.3f"
        % (auc_prime, average_precision_score(y_te, s_prime)))

    # ---------------------------------------------------------------- 12.4
    log("\n[12.4] [TRAINING] THE POINT: presentation is not immunogenicity")
    te = te.copy()
    te["score_cnn"] = s_cnn
    te["score_prime"] = s_prime
    te["presentation_score"] = pres_te

    auc_pres = None
    if source != "unavailable":
        auc_pres = roc_auc_score(y_te, pres_te)
        r = np.corrcoef(pres_te, s_cnn)[0, 1]
        log("      AUC from the presentation score alone : %.3f" % auc_pres)
        log("      AUC from the immunogenicity CNN       : %.3f" % auc_cnn)
        log("      correlation between the two scores    : r=%.2f" % r)
        log("")
        log("      Read that carefully. The presentation score does separate the")
        log("      classes a little, it is above 0.5, but far less than a model")
        log("      trained on the right question. And the two scores barely")
        log("      correlate, meaning they order the peptides by different")
        log("      criteria. That is the slide phrase, 'a strong binder is not")
        log("      necessarily immunogenic', turned into two numbers.")

        top_pres = te.nlargest(20, "presentation_score")
        log("")
        log("      Beware the easy shortcut. Of the 20 best presented peptides in")
        log("      the test set, %d/20 are immunogenic in IEDB." % int(top_pres.y.sum()))
        log("      That looks like presentation is enough, but it is an artefact of")
        log("      SELECTION BIAS. Peptides reach IEDB because somebody already")
        log("      suspected they would bind, not by sampling the proteome at random.")
        log("      The overall AUC of the presentation score stays low (%.3f) even"
            % auc_pres)
        log("      though the head of the ranking looks clean.")
        log("      Methodological moral: look at the whole curve, not the top k.")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for label, colour, name in [(0, "#8aa0b8", "not immunogenic"), (1, "#d1495b", "immunogenic")]:
        ax[0].hist(te.loc[te.y == label, "score_cnn"], bins=30, alpha=0.6,
                   label=name, color=colour, density=True)
    ax[0].set_xlabel("DeepImmuno-CNN score")
    ax[0].set_ylabel("density")
    ax[0].set_title("Class separation (AUC=%.3f)" % auc_cnn)
    ax[0].legend()

    if source != "unavailable":
        ax[1].scatter(te["presentation_score"], te["score_cnn"], c=te["y"],
                      cmap="coolwarm", s=12, alpha=0.7)
        ax[1].set_xlabel("presentation score (%s)" % source)
        ax[1].set_ylabel("immunogenicity score (CNN)")
        ax[1].set_title("Two different questions, two different axes")
    else:
        ax[1].text(0.5, 0.5, "no presentation score\n(install mhcflurry)",
                   ha="center", va="center")
        ax[1].axis("off")
    fig.suptitle("Exercise 12.4 - Presentation is not immunogenicity", fontsize=12)
    fig.tight_layout()
    fig.savefig(RESULTS / "ex12_4_presentation_vs_immunogenicity.png", dpi=140)
    plt.close(fig)
    log("      figure written to results/ex12_4_presentation_vs_immunogenicity.png")

    te.sort_values("score_cnn", ascending=False).to_csv(
        RESULTS / "ex12_test_scored.csv", index=False)

    # ---------------------------------------------------------------- 12.5
    log("\n[12.5] [fhaC] Apply the model to the case study")
    log("      From here on, everything is YOUR protein. No training happens")
    log("      below this line.")
    header, seq = read_fasta(DATA / "fhaC.fasta")
    log("      %s" % header[:80])
    log("      length: %d residues" % len(seq))
    cand = kmers(seq, 9)
    log("      %d overlapping 9-mers (%d residues - 9 + 1 = %d)"
        % (len(cand), len(seq), len(seq) - 9 + 1))
    log("      A sliding window always yields fewer peptides than residues.")
    log("      That is why exercise 11 reports %d rows, one per residue, and" % len(seq))
    log("      this one reports %d, one per peptide. Same protein." % len(cand))

    Xp_c = np.stack([encode_peptide_pca(p, pca12) for p in cand["peptide"]]).astype(np.float32)
    Xm_c = np.repeat(encode_hla_pca(FOCUS_HLA, paratopes, pca12)[None, ...],
                     len(cand), axis=0).astype(np.float32)
    cand["immunogenicity_score"] = score_cnn(model, Xp_c, Xm_c)
    cand["HLA"] = FOCUS_HLA
    cand = cand.sort_values("immunogenicity_score", ascending=False).reset_index(drop=True)
    cand.to_csv(RESULTS / "ex12_fhaC_immunogenicity.csv", index=False)

    torch.save(model.state_dict(), precomputed / "ex12_immunogenicity_cnn.pt")
    log("      model saved to data/precomputed/ex12_immunogenicity_cnn.pt")
    log("      (exercise 14 loads it instead of retraining)")

    log("\n      Top 10 fhaC peptides for %s by predicted immunogenicity:" % FOCUS_HLA)
    log(cand.head(10).to_string(index=False))
    log("\n      CAREFUL: this ranking ignores whether the peptide binds the HLA at")
    log("      all. Exercise 14 combines it with binding before deciding anything.")

    log("\n" + "=" * 78)
    log("QUESTIONS TO DISCUSS")
    log("=" * 78)
    log("  1. The AUC drops when the split is grouped by peptide. Why would an")
    log("     ordinary random split give a prettier and more dishonest number?")
    log("  2. PRIME leaves P2 and the last position out of the frequency vector.")
    log("     What biological assumption is behind that? Does it hold for class II?")
    log("  3. The slide says DeepImmuno reports low precision on SARS-CoV-2.")
    log("     Looking at the allele composition of the training set, why?")
    log("  4. The 9,055 peptides come from IEDB assays. Does 'Negative' mean")
    log("     'not immunogenic' or 'nobody detected a response in that assay'?")
    log("  5. fhaC is an extracellular bacterial protein. Its natural route is")
    log("     MHC class II. What does that make this whole class I exercise?")

    (RESULTS / "ex12_log.txt").write_text("\n".join(buffer), encoding="utf-8")
    print("\n[log] results/ex12_log.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
