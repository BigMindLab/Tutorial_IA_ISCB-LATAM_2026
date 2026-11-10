# Table formats, the four binding tools

What each wrapper writes, column by column, so the outputs can be compared
and so missing results can be filled in by hand.

Empty templates with the right headers are in `data/templates/`. Drop a filled
one into `results/` with the matching name and `ex5_compare.py` will pick it up.

---

## NetMHCpan 4.1, class I

File: `results/netmhcpan_class_I.csv`
Template: `data/templates/netmhcpan_class_I.csv`

| Column | Meaning |
|---|---|
| `start` | 1-based position of the peptide in the protein |
| `length` | peptide length, 8 to 11 as configured |
| `peptide` | the peptide itself |
| `allele` | HLA allele, written `HLA-A*02:01` |
| `score_EL` | eluted ligand score, trained on peptides actually eluted from cells |
| `rank_EL` | percentile rank of that score. **This is the one to read** |
| `score_BA` | binding affinity score, trained on measured affinities |
| `rank_BA` | percentile rank of the affinity score |
| `affinity_nM` | predicted affinity in nanomolar, lower is tighter |
| `call` | `SB`, `WB` or `NB` from `rank_EL` |

Thresholds from slide 21: `SB` below 0.5, `WB` below 2.

The EL and BA split is the point of slide 21. BA is trained on measured
affinities and EL on peptides eluted from real cells. They answer related but
different questions, and neither asks about T cells.

Example row:

```csv
start,length,peptide,allele,score_EL,rank_EL,score_BA,rank_BA,affinity_nM,call
9,9,RPGLVGRAL,HLA-B*07:02,0.8421,0.0310,0.7355,0.0620,18.4,SB
```

---

## NetMHCIIpan 4.0, class II

File: `results/netmhciipan_class_II.csv`
Template: `data/templates/netmhciipan_class_II.csv`

| Column | Meaning |
|---|---|
| `start` | 1-based position in the protein |
| `peptide` | the 15-mer |
| `allele` | written the NetMHCIIpan way, `DRB1_0401` |
| `score_EL` | eluted ligand score |
| `rank_EL` | percentile rank |
| `call` | `SB`, `WB` or `NB` |

Thresholds are **different** from class I: `SB` below 2, `WB` below 10.

There is no `length` column because every peptide is a 15-mer, and no
`affinity_nM` because the wrapper reads the eluted ligand output. The groove is
open at both ends, so the 15-mer is scored through a 9-mer binding core.

Example row:

```csv
start,peptide,allele,score_EL,rank_EL,call
103,RLFDPAPLVQDYLNR,DRB1_0401,0.5120,1.84,SB
```

---

## MixMHCpred 2.2, class I

File: `results/mixmhcpred_class_I.csv`

| Column | Meaning |
|---|---|
| `start` | 1-based position |
| `length` | peptide length |
| `peptide` | the peptide |
| `allele` | `HLA-A*02:01` |
| `score` | raw score, **negative values are normal** |
| `rank` | percentile rank, lower is better |
| `call` | `SB`, `WB` or `NB` |

Real rows from the fhaC run:

```csv
start,length,peptide,allele,score,rank,call
1,8,MTDATNRF,HLA-A*02:01,-3.578068,49.6938,NB
1,8,MTDATNRF,HLA-A*24:02,-3.611811,27.2604,NB
1,8,MTDATNRF,HLA-B*07:02,-4.104978,50.6827,NB
```

MixMHCpred is position weight matrices from motif deconvolution of mass
spectrometry ligands, not a neural network. Same question, different
machinery, which is exactly the contrast slide 21 draws.

---

## MHCflurry 2.0, class I

File: `results/mhcflurry_class_I.csv`

| Column | Meaning |
|---|---|
| `start`, `length`, `peptide`, `allele` | as above |
| `affinity_nM` | predicted affinity |
| `affinity_percentile` | rank of the affinity |
| `processing_score` | antigen processing, what happens before the groove |
| `presentation_score` | affinity and processing combined |
| `presentation_percentile` | rank of the presentation score |
| `call` | `SB`, `WB` or `NB` |

Real rows from the fhaC run:

```csv
start,length,peptide,allele,affinity_nM,affinity_percentile,processing_score,presentation_score,presentation_percentile,call
9,9,RPGLVGRAL,HLA-B*07:02,13.536504,0.002352,0.779360,0.988729,0.004721,SB
465,9,LLNSYQITV,HLA-A*02:01,11.682791,0.037398,0.489238,0.974643,0.017651,SB
140,9,YATSIVTFV,HLA-A*02:01,25.921557,0.131372,0.571021,0.955700,0.043591,SB
```

MHCflurry is the only one of the four with a processing term, which is why it
has both an affinity rank and a presentation rank.

---

## What the four have in common

Every tool gives `start`, `peptide`, `allele`, some score, a percentile rank
and a call. That shared core is what `ex5_compare.py` merges on, so a peptide
called `SB` by more than one tool can be spotted.

**Always compare on the rank, never on the raw score.** The scores are on
different scales. MixMHCpred is negative, MHCflurry presentation runs 0 to 1,
NetMHCpan has two separate scores. The percentile rank is the only column that
means the same thing across all four.

---

## Filling a template by hand

If a tool was run elsewhere and only the printed output survives, copy the
template, paste the rows and keep the column order. The comparison script does
not care how the file was produced, only that the headers match.

```bash
cp data/templates/netmhcpan_class_I.csv results/netmhcpan_class_I.csv
```

Then edit it, and run:

```bash
python3 solutions/ex5_compare.py
```
