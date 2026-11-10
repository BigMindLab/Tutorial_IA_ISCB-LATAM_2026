# M03: Immunoinformatics

**ISCB-LATAM 2026, Tutorial IP3**

This module follows a protein sequence through topology analysis, MHC binding
and presentation, immunogenicity scoring, and epitope prioritisation. The case
study is **FhaC**, the filamentous hemagglutinin transporter of *Bordetella
pertussis* (UniProtKB P35077, 584 residues).

## Learning goals

- Distinguish protein localisation, MHC binding, presentation, immunogenicity
  and TCR recognition.
- Encode peptides and HLA or TCR sequences for different models.
- Interpret affinity, percentile ranks and prediction scores.
- Combine predictions into a candidate list and identify missing evidence.

## Install and run

Use Linux or Ubuntu/WSL, Bash, and **Python 3.10-3.12** with `venv` and `pip`.
The exercises run on CPU. Internet access is required during setup.
See [QUICKSTART.md](QUICKSTART.md) for prerequisites, Conda and troubleshooting.

From the repository root:

```bash
cd modules/M03_immunoinformatics
bash setup.sh
source .venv/bin/activate
python solutions/ex00_fhaC_demo.py
```

If your terminal is already in this module directory, omit the `cd` command.
`setup.sh` creates or reuses `.venv`, installs the pinned packages, downloads
MHCflurry weights and public datasets, and checks imports and model inference.
Activate the environment in each new terminal. If the demo is started with a
Python that lacks its dependencies, it automatically retries with this
module's `.venv/bin/python` when available. Other exercises require the
selected environment to be active or its interpreter to be specified.

To select the module's Python without activating the environment:

```bash
.venv/bin/python solutions/ex00_fhaC_demo.py
```

The demo runs topology analysis and prioritisation with pretrained models.
**Exercises 12 and 13 are optional training exercises, not prerequisites for
the demo.** The installation does not include the external predictors listed
below, pVACtools or TIDE.

## Exercises

Run these commands from the module directory with the environment active.

| Command | Purpose | Trains a model? |
|---|---|---|
| `python solutions/ex00_fhaC_demo.py` | Run topology analysis and epitope prioritisation | No |
| `python solutions/ex11_topology_reference.py --deeptmhmm data/precomputed/deeptmhmm_fhaC.3line` | Compare the bundled DeepTMHMM prediction with UniProt annotation | No |
| `python solutions/ex3_mhcflurry.py data/fhaC.fasta` | Predict class I binding and presentation with MHCflurry | No |
| `python solutions/ex12_immunogenicity.py` | Train the tutorial immunogenicity model and compare scores | Yes |
| `python solutions/ex13_nettcr.py` | Train the tutorial TCR model and evaluate peptide generalisation | Yes |
| `python solutions/ex14_prioritization.py` | Rank candidates using topology, MHCflurry and pretrained immunogenicity scores | No |
| `python solutions/ex16_consensus.py` | Compare available binding tables and class I/class II regions | No |

Exercise 12 replaces `data/precomputed/ex12_immunogenicity_cnn.pt` with its
newly trained model. Subsequent runs of the demo and exercise 14 use those
weights. Training time depends on the CPU and the dataset.

Exercise 15 optionally requests predictions from the IEDB API:

```bash
python solutions/ex15_netmhcpan_iedb.py
python solutions/ex3_mhcflurry.py data/fhaC.fasta
python solutions/ex16_consensus.py
```

Exercise 15 requires internet access and an available remote service. Exercise
16 reads the NetMHCpan and NetMHCIIpan tables in `results/`; MHCflurry and
MixMHCpred tables can also be read from `data/precomputed/`. Bundled tables
apply to FhaC only.

## Outputs

Scripts write their outputs to `results/` and replace files with the same
names. The repository includes example tables, figures and logs.

| File | Contents |
|---|---|
| `results/topology_reference.csv` | Per-residue UniProt topology annotation |
| `results/mhcflurry_class_I.csv` | Class I binding and presentation predictions |
| `results/ex14_full_table.csv` | Scores for every candidate peptide/HLA pair |
| `results/ex14_shortlist.csv` | Candidates retained by the prioritisation filters |
| `results/ex14_prioritization.png` | Filtering steps and candidate ranking |
| `results/consensus_class_I.csv` | Consensus from exercise 16 |

For FhaC, exercise 14 generates 576 overlapping 9-mers and scores them against
five HLA alleles, giving 2880 peptide/HLA pairs before filtering. The number
retained depends on the scores, thresholds and model weights. Bundled results
are examples rather than exact numerical targets for every environment.

The tutorial models in exercises 12 and 13 are educational implementations;
they are not the original DeepImmuno or NetTCR packages. Exercise 14 implements
a prioritisation workflow without running pVACtools or TIDE. Predictions do
not establish experimental immunogenicity or TCR recognition.

## External predictors

These exercises require separate software or additional inputs. Follow the
installation and licence instructions for the release you use; do not mix its
requirements into the tutorial environment.

| Script | Additional requirement |
|---|---|
| `ex1_netmhcpan.py` | NetMHCpan; set `NETMHCPAN_BIN` |
| `ex2_netmhciipan.py` | NetMHCIIpan; set `NETMHCIIPAN_BIN` |
| `ex4_mixmhcpred.py` | MixMHCpred; set `MIXMHCPRED` |
| `ex6_prime.py` | PRIME; set `PRIME` |
| `ex7_signalp6.py` | SignalP 6.0 command on `PATH` |
| `ex8_deeptmhmm.py` | An existing DeepTMHMM `.3line` output, or BioLib |
| `ex9_deepimmuno.py` | Prepares input for the original DeepImmuno package |
| `ex10_nettcr_prepare.py` | Prepares peptide/CDR3 pairs; requires PRIME results and TCR sequences |

See [QUICKSTART.md](QUICKSTART.md#external-predictors) for commands and
[RUNBOOK.md](RUNBOOK.md) for the individual exercises.

## Colab

Open [notebooks/immunoinformatics_colab.ipynb](notebooks/immunoinformatics_colab.ipynb)
in Google Colab and run the cells in order. Use a runtime with Python
3.10-3.12. The notebook obtains the module by cloning or ZIP upload, then runs
`setup.sh` in the remote runtime. Local installation is not required; package
and data downloads still take place in Colab.

## Data and references

The FhaC sequence and annotation are in `data/fhaC.fasta` and
`data/P35077.json`. `utils/download_data.py` fetches the encoding tables and
training datasets from their upstream sources.

| Source | Use |
|---|---|
| [UniProtKB P35077](https://www.uniprot.org/uniprotkb/P35077) | FhaC sequence and curated annotation |
| [MHCflurry](https://github.com/openvax/mhcflurry) | Class I binding and presentation |
| [DeepImmuno](https://github.com/frankligy/DeepImmuno) | Encoding tables and labelled immunogenicity data |
| [NetTCR-2.0](https://github.com/mnielLab/NetTCR-2.0) | Paired TCR/peptide training data |

The allele panel for the binding wrappers is in `utils/common.py`.
Exercise 14 uses its own five-allele panel, `HLA_PANEL`, in
`solutions/ex14_prioritization.py`.
