"""Sequence encodings used by the trainable exercises (ex11, ex12, ex13).

Kept apart from common.py, which stays dependency free for the wrapper
exercises. Everything here needs numpy and pandas.

The point of this module is that "encoding a protein" means different things
to different models, and the slides say so explicitly:

  DeepImmuno (slide 22)  each residue becomes 12 PCA components of 553 AAindex
                         properties, and the MHC enters as its paratope residues
  NetTCR     (slide 23)  each residue becomes its BLOSUM50 row of 20 values,
                         zero padded to length 30
  PRIME 2.0  (slide 22)  the peptide collapses into residue frequencies at the
                         positions that face the TCR, plus length
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"

# Amino acid order used by DeepImmuno in after_pca.txt (20 residues plus X).
AA21 = "ARNDCQEGHILKMFPSTWYVX"
AA20 = "ACDEFGHIKLMNPQRSTVWY"


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------
def read_fasta(path: Path) -> tuple[str, str]:
    """Return (header, sequence) for the first record in a FASTA file."""
    header, chunks = "", []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    break
                header = line[1:]
            else:
                chunks.append(line)
    return header, "".join(chunks)


def kmers(sequence: str, k: int = 9) -> pd.DataFrame:
    """Sliding window of k-mers with 1-based start and end positions."""
    rows = [
        {"start": i + 1, "end": i + k, "peptide": sequence[i : i + k]}
        for i in range(len(sequence) - k + 1)
    ]
    df = pd.DataFrame(rows)
    # Predictors only accept the 20 standard amino acids.
    keep = df["peptide"].apply(lambda p: set(p) <= set(AA20))
    return df[keep].reset_index(drop=True)


# ---------------------------------------------------------------------------
# DeepImmuno encoding (slide 22)
# ---------------------------------------------------------------------------
def load_pca12(path: Path | None = None) -> pd.DataFrame:
    """DeepImmuno 21x12 matrix, one row per amino acid, 12 PCA components."""
    path = path or DATA / "after_pca.txt"
    m = pd.read_csv(path, sep=r"\s+", header=None)
    if m.shape != (21, 12):
        raise ValueError("after_pca.txt should be 21x12, got %s" % (m.shape,))
    m.index = list(AA21)
    m.columns = ["PC%d" % (i + 1) for i in range(12)]
    return m


def encode_peptide_pca(peptide: str, pca12: pd.DataFrame, width: int = 10) -> np.ndarray:
    """Peptide to a (width x 12) matrix.

    Peptides shorter than `width` are padded with X, row 21 of the PCA matrix.
    That is how DeepImmuno puts 9-mers and 10-mers through the same input.
    """
    p = peptide.upper()
    if len(p) > width:
        raise ValueError("peptide longer than %d: %s" % (width, peptide))
    p = p + "X" * (width - len(p))
    return np.stack([pca12.loc[aa if aa in pca12.index else "X"].to_numpy() for aa in p])


def load_paratopes(path: Path | None = None) -> pd.DataFrame:
    """HLA to its 46 residue paratope pseudo sequence (DeepImmuno table)."""
    path = path or DATA / "hla2paratopeTable_aligned.txt"
    t = pd.read_csv(path, sep="\t")
    t.columns = ["HLA", "pseudo"]
    return t.set_index("HLA")


def encode_hla_pca(hla: str, paratopes: pd.DataFrame, pca12: pd.DataFrame) -> np.ndarray:
    """HLA to a (46 x 12) matrix. Alignment gaps '-' are encoded as X."""
    pseudo = paratopes.loc[hla, "pseudo"]
    return np.stack(
        [pca12.loc[aa if aa in pca12.index else "X"].to_numpy() for aa in pseudo]
    )


# ---------------------------------------------------------------------------
# NetTCR encoding (slide 23)
# ---------------------------------------------------------------------------
def load_blosum50() -> pd.DataFrame:
    """BLOSUM50 as a 20x20 frame, taken from Biopython rather than hardcoded."""
    from Bio.Align import substitution_matrices

    try:
        m = substitution_matrices.load("BLOSUM50")
    except Exception:  # noqa: BLE001
        print("[warning] BLOSUM50 unavailable, falling back to BLOSUM62")
        m = substitution_matrices.load("BLOSUM62")
    df = pd.DataFrame(np.array(m), index=list(m.alphabet), columns=list(m.alphabet))
    return df.loc[list(AA20), list(AA20)]


def encode_blosum(sequence: str, blosum: pd.DataFrame, width: int = 30) -> np.ndarray:
    """Sequence to a (width x 20) matrix with BLOSUM50 and zero padding.

    Values are divided by 5 to keep them in a comfortable range for the
    network, which is what NetTCR does as well.
    """
    s = sequence.upper()[:width]
    m = np.zeros((width, 20), dtype=np.float32)
    for i, aa in enumerate(s):
        if aa in blosum.index:
            m[i] = blosum.loc[aa].to_numpy() / 5.0
    return m


# ---------------------------------------------------------------------------
# PRIME 2.0 encoding (slide 22): 28 inputs
# ---------------------------------------------------------------------------
def tcr_facing_positions(length: int) -> list[int]:
    """1-based positions that face the TCR: everything except P2 and the last.

    P2 and the C-terminal residue are the anchors that sit in the MHC groove.
    """
    return [p for p in range(1, length + 1) if p != 2 and p != length]


def prime_features(peptide: str, presentation_score: float) -> np.ndarray:
    """The 28 input vector of PRIME 2.0 (slide 22).

      [0]      1 presentation score
      [1:21]   20 residue frequencies at the TCR facing positions
      [21:28]  7 one-hot nodes for peptide length 8 to 14
    """
    p = peptide.upper()
    v = np.zeros(28, dtype=np.float32)
    v[0] = presentation_score

    facing = [p[i - 1] for i in tcr_facing_positions(len(p))]
    if facing:
        for aa in facing:
            if aa in AA20:
                v[1 + AA20.index(aa)] += 1.0 / len(facing)

    if 8 <= len(p) <= 14:
        v[21 + (len(p) - 8)] = 1.0
    return v


def binary_label(value: str) -> int:
    """IEDB labels are Negative or Positive with optional Low/High/Intermediate."""
    return 0 if str(value).strip() == "Negative" else 1
