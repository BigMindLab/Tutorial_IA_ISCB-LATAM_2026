# QUICKSTART, module M03.1

Two paths. Pick the one that matches what you have.

- **Path 1, Colab.** Nothing to install. Covers exercises 11 to 14, which is
  slides 18, 19, 22, 23 and 24. Start here if you are short on time.
- **Path 2, local install on Linux.** Adds the licensed tools, so you also get
  NetMHCpan, NetMHCIIpan, SignalP 6.0, DeepTMHMM, MixMHCpred and PRIME running
  on your own machine.

Everything uses **fhaC**, which ships in `data/fhaC.fasta`.

---

## Path 1, Colab

Open `notebooks/immunoinformatics_colab.ipynb` in Google Colab and run the
cells top to bottom. CPU is enough. No GPU is needed anywhere in this module.

Budget 25 to 35 minutes, most of it MHCflurry scoring thousands of peptides.

---

## Path 2, local

### Base environment

Python 3.10 or newer for track B.

```bash
python3 -m venv venv
```

```bash
source venv/bin/activate
```

```bash
pip install -r requirements.txt
```

```bash
mhcflurry-downloads fetch models_class1_presentation
```

```bash
python utils/download_data.py
```

Check it works:

```bash
python solutions/ex11_topology_reference.py
```

That should print the fhaC signal peptide at 1 to 30 and 16 transmembrane beta
strands. If it does, track B is ready.

### Licensed tools

Each one is obtained by you, under your own licence. None of them is
redistributed in this repository.

**NetMHCpan 4.1 and NetMHCIIpan 4.0.** Request an academic licence from DTU
Health Tech and follow their installer. Then point the wrappers at the
binaries:

```bash
export NETMHCPAN_BIN=/path/to/netMHCpan
```

```bash
export NETMHCIIPAN_BIN=/path/to/netMHCIIpan
```

**SignalP 6.0.** Also DTU, academic licence and manual download. It provides a
`signalp6` command. For fhaC the organism argument is `other`, not `eukarya`,
because SignalP 6.0 only separates eukaryotes from everything else.

**DeepTMHMM.** Two options, and both have a catch.

The BioLib route no longer runs locally:

```bash
pip install pybiolib
```

```bash
biolib run DTU/DeepTMHMM --fasta data/fhaC.fasta
```

This now requires a BioLib account, and `machine="local"` was removed, so the
job runs in their cloud.

The academic package from DTU does run offline. Two things to know before you
start. It pins `fair-esm==0.4.0`, and the alphabet object saved inside the
package is not compatible with newer releases, so that exact version matters.
It also calls `plt.style.use('seaborn-whitegrid')`, a style name matplotlib
removed, so on a modern matplotlib it crashes at the very last step **after**
writing `predicted_topologies.3line` and `TMRs.gff3`. Your results are already
on disk when that happens.

A verified DeepTMHMM run on fhaC ships in
`data/precomputed/deeptmhmm_fhaC.3line`, so the exercise works even if you
never get the tool installed.

**MixMHCpred 2.2.** Clone the `v2.2` tag, not the default branch, which is 3.0:

```bash
git clone --branch v2.2 https://github.com/GfellerLab/MixMHCpred.git
```

Compile it, then:

```bash
export MIXMHCPRED=/path/to/MixMHCpred
```

**PRIME 2.0.** Same group, same pattern:

```bash
git clone --branch v2.0 https://github.com/GfellerLab/PRIME.git
```

```bash
export PRIME=/path/to/PRIME
```

**DeepImmuno.** Clone it and keep the path, because it has to be run from its
own root:

```bash
git clone https://github.com/frankligy/DeepImmuno.git
```

---

## Running the exercises

Order matters in one place only: exercise 14 loads the model exercise 12 saves.

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

Once SignalP and DeepTMHMM are installed, feed their output back into
exercise 11 to score them against the curated reference:

```bash
python solutions/ex11_topology_reference.py --signalp results/signalp6/prediction_results.txt --deeptmhmm biolib_results/predicted_topologies.3line
```

For the licensed binding tools:

```bash
python solutions/ex1_netmhcpan.py data/fhaC.fasta
```

```bash
python solutions/ex2_netmhciipan.py data/fhaC.fasta
```

```bash
python solutions/ex3_mhcflurry.py data/fhaC.fasta
```

```bash
python solutions/ex4_mixmhcpred.py data/fhaC.fasta
```

```bash
python solutions/ex5_compare.py
```

Alleles are edited in `utils/common.py`, in the `CLASS_I_ALLELES` and
`CLASS_II_ALLELES` lists.

---

## Known problems

**`ModuleNotFoundError: No module named 'common'`.** Run the scripts from
inside `solutions/`, or let the path shim at the top of each file do its job.
The shared helpers live in `utils/`.

**DeepImmuno scores a peptide called "peptide".** Your input CSV has a header
row. It must not have one, because DeepImmuno reads it with `header=None`.
`ex9_deepimmuno.py` already writes the file correctly.

**DeepImmuno cannot find its models.** Run it from the DeepImmuno repository
root. It loads `./data/` and `./models/` by relative path.

**MHCflurry refuses a peptide.** It only accepts the 20 standard amino acids.
The helpers filter out anything else before predicting.

**Results land somewhere unexpected.** Everything is written to `results/` at
the module root, no matter which directory you launched from.
