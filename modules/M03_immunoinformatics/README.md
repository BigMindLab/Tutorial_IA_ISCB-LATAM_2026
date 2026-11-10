# M03: Immunoinformatics

**Time:** 10:45-11:30 | **Type:** Hands-on

From a protein sequence to a prioritised shortlist of T cell epitope candidates
with SignalP 6.0, DeepTMHMM, NetMHCpan, NetMHCIIpan, MHCflurry, MixMHCpred,
DeepImmuno, PRIME, NetTCR and pVACtools.
Case study: **fhaC**, the filamentous hemagglutinin transporter of
*Bordetella pertussis* (UniProtKB
[P35077](https://www.uniprot.org/uniprotkb/P35077), 584 aa), the same protein
used by the structural module.

| Script | Content |
|----------|---------|
| `solutions/ex00_fhaC_demo.py` | The whole workflow on fhaC, one command, no training |
| `solutions/ex11_topology_reference.py` | Exercise 1: where does fhaC go, scored against curated annotation |
| `solutions/ex14_prioritization.py` | Exercise 2: from peptides to a ranked shortlist |
| `solutions/ex12_immunogenicity.py` | Complementary: inside DeepImmuno and PRIME |
| `solutions/ex13_nettcr.py` | Complementary: inside NetTCR, and where it breaks |
| `notebooks/immunoinformatics_colab.ipynb` | All of the above in Colab, nothing to install |

## Learning goals

By the end of this module, participants will be able to:

- State what each model is trained to predict, and what it is not.
- Encode a protein or peptide the way each model requires.
- Read binding affinity, %rank and presentation scores correctly.
- Keep binding, presentation, immunogenicity and TCR recognition apart.
- Combine predictions into a ranked, defensible candidate list.

## Session plan

| Activity | Time |
|---|---:|
| Theory, the five questions | 15 min |
| Exercise 1: topology of fhaC and what it rules out | 10 min |
| Exercise 2: binding, immunogenicity and the shortlist | 14 min |
| Discussion | 6 min |
| **Total** | **45 min** |

---

## The five questions

Every arrow narrows the list and no arrow proves a T cell response.

| Slide | Question | Tools | Exercise |
|---|---|---|---|
| 18 | Does the protein enter the secretory pathway? | SignalP 6.0 | `ex7_signalp6.py`, `ex11_topology_reference.py` |
| 19 | Is it a membrane protein, and which way round? | DeepTMHMM | `ex8_deeptmhmm.py`, `ex11_topology_reference.py` |
| 20, 21 | Which peptides can MHC bind? | NetMHCpan, NetMHCIIpan, MHCflurry, MixMHCpred | `ex1` to `ex6` |
| 22 | Which are immunogenic? | DeepImmuno, PRIME 2.0 | `ex9_deepimmuno.py`, `ex12_immunogenicity.py` |
| 23 | Which TCR recognises them? | NetTCR | `ex10_nettcr_prepare.py`, `ex13_nettcr.py` |
| 24 | Which do we take forward? | pVACtools, TIDE | `ex14_prioritization.py` |

---

## Which script do I run?

**To show the biology of fhaC, run one script:**

```bash
python solutions/ex00_fhaC_demo.py
```

That is the whole workflow of slides 18 to 24 on fhaC, and **nothing is
trained while it runs**. It loads pretrained weights from disk, exactly the
way MHCflurry loads its own. Every number it prints comes from the case study
protein.

**To look inside the models, run the two that train:**

```bash
python solutions/ex12_immunogenicity.py
```

```bash
python solutions/ex13_nettcr.py
```

These exist because two points cannot be made by prediction alone. Exercise 12
shows that presentation and immunogenicity are different questions, as two AUC
numbers on the same peptides. Exercise 13 shows that NetTCR must not be used on
fhaC, by watching its AUC collapse. Both train on published datasets and then
apply the result to fhaC, and every block of their output is tagged `[TRAINING]`
or `[fhaC]` so the two never get confused.

### A note on counting

The same protein is reported three different ways, and it trips people up:

| Unit | Count | Where |
|---|---|---|
| residues | 584 | exercise 11 |
| overlapping 9-mers | 576 | exercises 12 and 14 |
| 9-mers by 5 HLA alleles | 2880 | exercise 14 |

576 is 584 minus 9 plus 1. A sliding window always yields fewer peptides than
residues, because the last 8 residues cannot start another full 9-mer. It is
one protein throughout.

### Why nothing trains on fhaC

Nobody has ever measured which fhaC peptides are immunogenic, or which TCR
recognises them. fhaC carries no labels, so there is nothing in it to learn
from. The models are trained on published labelled data and then applied to
fhaC, which is what the tools themselves do. Slide 22 says DeepImmuno was
trained on 9,056 IEDB peptides.

---

## Two tracks, and why

Several tools on these slides need an academic licence and a local install, so
a workshop that depends on them stalls before it starts. The exercises are
therefore split in two.

**Track A, wrappers for licensed tools.** `ex1` to `ex10` drive NetMHCpan,
NetMHCIIpan, MixMHCpred, PRIME, SignalP, DeepTMHMM and DeepImmuno. Participants
who installed the software on Linux run these and get real output from the
tools the slides actually show.

**Track B, reimplementations that always run.** `ex11` to `ex14` need nothing
but pip packages. They rebuild DeepImmuno-CNN and NetTCR-2.0 from the
architecture drawn on the slides, train them on the authors' own published
data, and reproduce the exact point each slide is making. These work the moment
the repository is cloned, in Colab or on a laptop.

Both tracks run on fhaC, so their outputs can be compared directly. They agree
on several peptides, which is worth showing in the room.

---

## Exercises

### Track A, needs the licensed tool installed

| Script | Tool | Status |
|---|---|---|
| `ex1_netmhcpan.py` | NetMHCpan 4.1, class I | wrapper, not executed here |
| `ex2_netmhciipan.py` | NetMHCIIpan 4.0, class II | wrapper, not executed here |
| `ex3_mhcflurry.py` | MHCflurry 2.0 | **runs, verified on fhaC** |
| `ex4_mixmhcpred.py` | MixMHCpred 2.2 | wrapper, output included |
| `ex5_compare.py` | compares the binding tools | runs on whatever output exists |
| `ex6_prime.py` | PRIME 2.0 | wrapper, output included |
| `ex7_signalp6.py` | SignalP 6.0 | wrapper, not executed here |
| `ex8_deeptmhmm.py` | DeepTMHMM | wrapper, see note below |
| `ex9_deepimmuno.py` | DeepImmuno-CNN | input preparation, flags verified |
| `ex10_nettcr_prepare.py` | NetTCR | input preparation, needs your own CDR3s |

### Track B, always runs

| Script | Slide | What it builds | Trains? |
|---|---|---|---|
| `ex00_fhaC_demo.py` | 18 to 24 | The whole workflow on fhaC, one command | no |
| `ex11_topology_reference.py` | 18, 19 | A curated topology reference from UniProt to score SignalP and DeepTMHMM against | no |
| `ex12_immunogenicity.py` | 22 | DeepImmuno-CNN and the PRIME 28 input vector, trained on 9,055 IEDB peptides | yes |
| `ex13_nettcr.py` | 23 | NetTCR-2.0, plus the leave one peptide out experiment that finds its limit | yes |
| `ex14_prioritization.py` | 24 | The pVACtools style filter chain over fhaC, and the contrast with TIDE | no |

---

## Reference results

These came out of a full run. Decimals move a little with the PyTorch version
and the batch shuffle. The comparisons do not move.

**Exercise 12, slide 22**

| Measure | Value |
|---|---|
| AUC of DeepImmuno-CNN, split grouped by peptide | 0.853 |
| AUC from the presentation score alone | 0.582 |
| Correlation between the two scores | r = 0.13 |

Two models ordering the same peptides almost independently. That is
"a strong binder is not necessarily immunogenic" expressed as two numbers.

**Exercise 13, slide 23**

| Measure | Value |
|---|---|
| Global AUC on the held out partition | 0.868 |
| AUC on GILGFVFTL, with the peptide in training | 0.897 |
| AUC on GILGFVFTL, without it (leave one peptide out) | 0.569 |

Per peptide the AUC falls to between 0.20 and 0.65 for peptides with around
100 training TCRs. With 17 peptides in total the model could not have learned
anything else.

**Exercise 14, slide 24**

| Filter step | Peptide and allele pairs |
|---|---|
| all | 2880 |
| outside the signal peptide | 2770 |
| binding, %rank below 2 | 163 |
| presentation, %rank below 2 | 109 |
| immunogenicity above the median | **54** |

Top candidate `RLFDPAPLV` at positions 103 to 111, inside the POTRA domain,
HLA-A\*02:01, 14.9 nM, %rank 0.059. The `TCR_evidence` column stays empty for
all 54 on purpose, and exercise 13 explains why filling it would be inventing
information.

**Cross check worth showing.** MixMHCpred, MHCflurry and the combined ranking
converge on several of the same peptides, among them `RPGLVGRAL`, `LLNSYQITV`
and `RLFDPAPLV`. Three independent tools agreeing is a far better argument than
any single ranking.

---

## Running it

See [QUICKSTART.md](QUICKSTART.md) for installation, including the licensed
tools. The short version for track B:

```bash
pip install -r requirements.txt
```

```bash
mhcflurry-downloads fetch models_class1_presentation
```

```bash
python utils/download_data.py
```

```bash
python solutions/ex11_topology_reference.py
```

```bash
python solutions/ex12_immunogenicity.py
```

```bash
python solutions/ex13_nettcr.py
```

```bash
python solutions/ex14_prioritization.py
```

Exercise 14 loads the model exercise 12 saves, so that order matters. On CPU
the whole of track B takes roughly 25 to 35 minutes, most of it MHCflurry
scoring thousands of peptides.

The Colab notebook in `notebooks/` runs the same thing with nothing installed
locally.

---

## Notes and gotchas found while building this

**DeepImmuno input has no header row.** `deepimmuno-cnn.py` reads the input
with `pd.read_csv(upload, sep=',', header=None)` and names the columns itself.
A header line gets scored as a peptide called "peptide". `ex9` writes the file
without one.

**DeepImmuno must run from its own repository root.** It loads `./data/` and
`./models/` by relative path and fails from anywhere else. Its `--intdir` flag
takes a file, not a directory.

**DeepTMHMM no longer runs locally through BioLib.** `biolib run` now requires
an account and local execution was dropped. The academic package from DTU does
run offline, but it pins `fair-esm==0.4.0` and the saved alphabet object is
incompatible with newer releases, so that exact version matters.

**fhaC reaches the immune system as an extracellular protein.** Its natural
route is endocytosis and MHC class II presentation to CD4 T cells. Every class
I result here is, strictly, a hypothesis about cross presentation. That is why
`ex2_netmhciipan.py` is not optional for this antigen, and it is a good
discussion to open on slide 20.

---

## Data and credits

| File | Source | Reference |
|---|---|---|
| `fhaC.fasta`, `P35077.json` | UniProtKB | FhaC, *B. pertussis* Tohama I |
| `after_pca.txt`, `hla2paratopeTable_aligned.txt`, `remove0123.csv` | [DeepImmuno](https://github.com/frankligy/DeepImmuno) | Li *et al.*, *Brief Bioinform* 22(6):bbab160 (2021) |
| `train_ab_90_alphabeta.csv` | [NetTCR-2.0](https://github.com/mnielLab/NetTCR-2.0) | Montemurro *et al.*, *Commun Biol* 4, 1060 (2021) |
| binding and presentation predictions | [MHCflurry 2.0](https://github.com/openvax/mhcflurry) | O'Donnell *et al.*, *Cell Syst* (2020) |

Architectures reimplemented from DeepImmuno-CNN (Li *et al.* 2021, Fig. 1C),
PRIME 2.0 (Gfeller *et al.*, *Cell Syst* 14:72-83, 2023, Fig. 3B) and
NetTCR-2.0 (Montemurro *et al.* 2021, Fig. 10). Underlying immunology from
Pishesha, Harmand and Ploegh, *Nat Rev Immunol* 22:751-764 (2022).

Licensed software is never redistributed here. Participants install
NetMHCpan, NetMHCIIpan, SignalP and DeepTMHMM under their own DTU academic
licence, and MixMHCpred and PRIME from their own repositories.

**Academic and teaching use.** None of these predictions is validated for
clinical decisions.
