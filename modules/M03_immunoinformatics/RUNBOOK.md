# RUNBOOK, module M03

What to type, what goes in, what comes out, and what each number means.

One protein throughout: **fhaC**, 584 residues, already in `data/fhaC.fasta`.

---

## Before anything

```bash
cd modules/M03_immunoinformatics
```

```bash
bash setup.sh
```

That installs the packages, downloads the MHCflurry weights and fetches the
public datasets. Around five minutes, mostly the weights.

It does **not** install NetMHCpan, NetMHCIIpan, SignalP or DeepTMHMM. Those
need an academic licence from DTU and a manual download. Every exercise still
runs, because the outputs of those tools for fhaC are in `data/precomputed`.

---

## Part A, the tools you already have installed

These four are the ones the team has licences for. Run them on your own
machine, not in Codespaces, because the licensed binaries cannot be installed
there.

### A1. SignalP 6.0, does fhaC enter the secretory pathway?

```bash
python3 solutions/ex7_signalp6.py data/fhaC.fasta other
```

| | |
|---|---|
| **Input** | `data/fhaC.fasta`, plus the organism group. For fhaC it is `other`, because SignalP 6.0 only separates eukaryotes from everything else and *B. pertussis* is Gram negative. |
| **Output** | `results/signalp6_summary.csv` with the predicted type, the per class probabilities and the cleavage position. The per residue n, h and c region probabilities land in `results/signalp6/`. |
| **What it means** | A signal peptide means the protein is secreted, so those residues are cut off and are not part of the mature protein. It says nothing about antigenicity. |
| **What to look at** | The cleavage position. UniProt annotates it at residue 30 for fhaC. Anything within a residue or two is agreement. |

Three ways to get SignalP, in order of reliability:

1. **Local install** with the DTU academic licence. This is the one to use.
2. **The DTU web server** at services.healthtech.dtu.dk. Submitting works and
   the job queues and runs, but it returned a failed run for fhaC on two
   separate attempts while this module was being built, with a backend error
   and no output files produced. Worth a try, not worth planning around.
3. **Skip it.** DeepTMHMM already labels the signal peptide, and the committed
   run puts the boundary at residue 32 against UniProt's 30. What SignalP adds
   on top is the signal peptide *type*, one of the five on slide 18, and the
   per residue n, h and c region probabilities. The cleavage site itself is
   already answered twice over.

### A2. DeepTMHMM, is it a membrane protein and which way round?

```bash
python3 solutions/ex8_deeptmhmm.py data/fhaC.fasta
```

| | |
|---|---|
| **Input** | `data/fhaC.fasta` |
| **Output** | `results/deeptmhmm_summary.txt`, plus the tool's own `predicted_topologies.3line` and `TMRs.gff3` |
| **What it means** | Every residue gets a label. S signal peptide, I inside, O outside, M alpha helical membrane, B beta barrel membrane, P periplasm. |
| **What to look at** | The number of B segments. fhaC is an outer membrane transporter, so expect a beta barrel. A verified run is already in `data/precomputed/deeptmhmm_fhaC.3line` and gives **16 strands**. |

Two gotchas. BioLib now needs an account and no longer runs locally. The DTU
academic package runs offline but pins `fair-esm==0.4.0`, and it crashes on the
final plot with a modern matplotlib **after** writing the results, so your
output is already on disk when that happens.

### A3. Score those two against curated annotation

```bash
python3 solutions/ex11_topology_reference.py --signalp results/signalp6/prediction_results.txt --deeptmhmm data/precomputed/deeptmhmm_fhaC.3line
```

| | |
|---|---|
| **Input** | the two outputs above. Both flags are optional, and with neither you get the reference alone. |
| **Output** | `results/topology_reference.csv`, one row per residue, plus a printed report |
| **What it means** | The reference is the curated UniProt annotation of P35077, partly from a solved structure. This is where you find out whether to believe the predictors. |
| **What to look at** | Strand count predicted versus annotated, 16 against 16. The cleavage site, 32 against 30. Agreement on the strand or not strand question, 88.4 percent. |

Do not be alarmed by the raw per residue agreement of 47.3 percent. UniProt
does not label inside versus outside, so every residue DeepTMHMM correctly
calls periplasm counts as a mismatch. The 88.4 percent figure is the fair one.

**If your protein is not in UniProt**, this exercise says so and tells you what
to do instead. Exercise 14 falls back to DeepTMHMM automatically.

### A4. NetMHCpan 4.1, class I binding

```bash
export NETMHCPAN_BIN=/path/to/netMHCpan
```

```bash
python3 solutions/ex1_netmhcpan.py data/fhaC.fasta
```

| | |
|---|---|
| **Input** | `data/fhaC.fasta`, tiled into 8 to 11-mers. Alleles come from `CLASS_I_ALLELES` in `utils/common.py`. |
| **Output** | `results/netmhcpan_class_I.csv` with affinity, %rank and an SB, WB or NB call |
| **What it means** | How well one peptide sits in the groove of one HLA allele. It does not model the proteasome, TAP or the cell. |
| **What to look at** | The %rank, not the nM value. Slide 21: strong binder below 0.5, weak binder below 2. |

### A5. NetMHCIIpan 4.0, class II binding

```bash
export NETMHCIIPAN_BIN=/path/to/netMHCIIpan
```

```bash
python3 solutions/ex2_netmhciipan.py data/fhaC.fasta
```

| | |
|---|---|
| **Input** | `data/fhaC.fasta`, tiled into 15-mers. Alleles from `CLASS_II_ALLELES`. |
| **Output** | `results/netmhciipan_class_II.csv` |
| **What it means** | Class II presents to CD4 helper T cells, and the groove is open at both ends, so a 15-mer is scored through a 9-mer core. |
| **What to look at** | Different thresholds. Class II strong binder below 2, weak binder below 10. |

**This is the most important exercise for fhaC and the easiest to skip.** fhaC
is an extracellular bacterial protein, so its natural route is endocytosis and
class II presentation. Everything we do in class I is, strictly, a hypothesis
about cross presentation.

---

## Part B, runs anywhere including Codespaces

Nothing licensed, nothing trained.

### B1. The whole story in one command

```bash
python3 solutions/ex00_fhaC_demo.py
```

| | |
|---|---|
| **Input** | `data/fhaC.fasta` and the pretrained weights in `data/precomputed` |
| **Output** | the shortlist, the full table and the figures, all in `results/` |
| **Runtime** | a few minutes, almost all MHCflurry |

This is the one to show. It chains topology and prioritisation and prints
which part is which.

### B2. MHCflurry, the free stand in for NetMHCpan

```bash
python3 solutions/ex3_mhcflurry.py data/fhaC.fasta
```

| | |
|---|---|
| **Output** | `results/mhcflurry_class_I.csv`, 6906 rows, 60 strong binders and 179 weak |
| **Why it is here** | Same question, different machinery, which is the contrast slide 21 draws with MixMHCpred. Compare its top hits against NetMHCpan when you have both. |

### B3. Compare the binding tools

```bash
python3 solutions/ex5_compare.py
```

| | |
|---|---|
| **Input** | whatever binding outputs exist in `results/` |
| **Output** | `results/class_I_merged.csv`, `results/class_I_consensus.csv` and `results/summary.txt` |
| **What to look at** | Which peptides more than one tool calls a strong binder. Agreement between independent methods is worth far more than any single ranking. |

---

## Part C, what is inside the models

These two train. Run them when there is time for the methodology discussion,
not when the goal is the biology of fhaC. Every block of output is tagged
`[TRAINING]` or `[fhaC]`.

### C1. Immunogenicity, slide 22

```bash
python3 solutions/ex12_immunogenicity.py
```

| | |
|---|---|
| **Input** | 9055 IEDB peptides from the DeepImmuno repository, then fhaC |
| **Output** | `results/ex12_fhaC_immunogenicity.csv`, two figures, and the model weights |
| **The result that matters** | AUC 0.853 for the immunogenicity model against 0.582 for the presentation score alone, correlation 0.13 between them |
| **What it means** | Binding and immunogenicity are different questions. A strong binder is not necessarily immunogenic. |

### C2. TCR recognition, slide 23

```bash
python3 solutions/ex13_nettcr.py
```

| | |
|---|---|
| **Input** | 16464 TCR and peptide pairs from NetTCR-2.0, covering only 17 peptides |
| **Output** | `results/ex13_auc_per_peptide.csv` and the limits figure |
| **The result that matters** | AUC 0.897 on GILGFVFTL with the peptide in training, 0.569 without it |
| **What it means** | The model only works on the peptides it was trained on. For fhaC the honest answer is no evidence. |

---

## Part D, the TCR input

NetTCR is the only tool here that needs something a protein sequence can never
provide: CDR3 alpha and CDR3 beta from real T cell sequencing.

```bash
python3 solutions/ex10_nettcr_prepare.py data/tcrs_example.csv 20
```

| | |
|---|---|
| **Input** | `data/tcrs_example.csv`, 20 real published CDR3 pairs, plus the PRIME candidates |
| **Output** | `results/nettcr_input.csv`, ready for NetTCR itself |
| **Important** | Those TCRs were measured against other peptides, never against fhaC. The `measured_against` column says which. They are there to show the format. |

See `data/README_tcrs.md` for how to swap in your own, and what experiment
would be needed to fill that column for real.

---

## What every output file is

| File | What it holds |
|---|---|
| `results/topology_reference.csv` | one row per residue of fhaC with its reference label |
| `results/mhcflurry_class_I.csv` | every peptide and allele pair with affinity, %rank and call |
| `results/class_I_consensus.csv` | peptides that more than one binding tool agrees on |
| `results/ex12_fhaC_immunogenicity.csv` | all 576 nine-mers ranked by predicted immunogenicity |
| `results/ex13_auc_per_peptide.csv` | how well NetTCR does per peptide, and on how much data |
| `results/ex14_shortlist.csv` | **the deliverable**, 54 candidates with every score behind them |
| `results/ex14_full_table.csv` | all 2880 pairs before filtering, for anyone who wants to re-filter |
| `results/*_log.txt` | the full printed output of each exercise, including the discussion questions |

---

## Why the counts differ between exercises

Same protein, different unit of analysis. This confuses everyone once.

| Unit | Count |
|---|---|
| residues | 584 |
| overlapping 9-mers | 576 |
| 9-mers by 5 HLA alleles | 2880 |
| after all filters | 54 |

576 is 584 minus 9 plus 1. A sliding window always gives fewer peptides than
residues, because the last 8 residues cannot start another full 9-mer.

---

## Running it in GitHub Codespaces

The repository already has a devcontainer, so Codespaces works. Open it, then:

```bash
cd modules/M03_immunoinformatics && bash setup.sh
```

```bash
python3 solutions/ex00_fhaC_demo.py
```

What works there: everything in parts B, C and D. What does not: the licensed
tools in part A, because they cannot be installed in a Codespace. Their outputs
for fhaC are committed in `data/precomputed`, so the exercises that consume
them still run.

On the default two CPU machine the demo takes a few minutes and the two
training exercises take noticeably longer. If you plan to run part C live,
pick a four core machine when creating the Codespace.
