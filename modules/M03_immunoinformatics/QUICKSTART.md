# QUICKSTART: M03 Immunoinformatics

## Prerequisites

- Linux or Ubuntu/WSL and Bash. On Windows, run the commands in the Ubuntu/WSL
  terminal.
- Python **3.10, 3.11 or 3.12**, with `venv` and `pip`, or Conda/Mamba.
- Internet access to download packages, MHCflurry models and public datasets.
- Disk space for the Python environment, model cache and downloaded inputs.

A GPU is not required. Commands below assume the module directory is the
working directory. From the repository root, enter it with:

```bash
cd modules/M03_immunoinformatics
```

## Install with a virtual environment

Check the Python version:

```bash
python3 --version
```

On Ubuntu, if the selected Python cannot create a virtual environment, install
its matching venv package. For Ubuntu's Python 3.10, for example:

```bash
sudo apt install python3.10-venv
```

Then install and run:

```bash
bash setup.sh
source .venv/bin/activate
python solutions/ex00_fhaC_demo.py
```

The installer creates or reuses `.venv`, installs pinned dependencies, runs
`pip check`, downloads MHCflurry presentation weights and public inputs, and
checks model inference. It does not install external predictors such as
SignalP or NetMHCpan. Linux x86_64 uses the CPU PyTorch wheel.

If you need to choose an installed Python when creating a new environment:

```bash
PYTHON=python3.11 bash setup.sh --env-dir .venv-py311
source .venv-py311/bin/activate
python solutions/ex00_fhaC_demo.py
```

An existing environment keeps its own interpreter. Use a new `--env-dir` to
select a different Python version. Relative environment paths are resolved
from the module directory.

## Install with Conda or Mamba

Use this option when a compatible Python or `venv` is unavailable. From the
module directory:

```bash
conda env create -f environment.yml
conda activate m03-immunoinformatics
bash setup.sh --use-current-env
python solutions/ex00_fhaC_demo.py
```

`environment.yml` selects Python 3.11 and installs the same requirements as the
virtual environment. `--use-current-env` reuses the active environment instead
of creating `.venv`. It requires an active virtual environment or Conda prefix.
Use `mamba env create -f environment.yml` if you manage environments with Mamba.

For later sessions, activate the environment before running exercises:

```bash
conda activate m03-immunoinformatics
python solutions/ex00_fhaC_demo.py
```

## Run after installation

For the virtual environment, activate it in each new terminal:

```bash
source .venv/bin/activate
python solutions/ex00_fhaC_demo.py
```

If the current Python lacks the demo dependencies, the demo automatically
restarts with the module's `.venv/bin/python` when it exists. This includes
starting from Miniconda's base Python. An incomplete `.venv` reports the
missing packages and stops; run setup to finish its installation.

You can also select the module's Python explicitly:

```bash
.venv/bin/python solutions/ex00_fhaC_demo.py
```

If you used a custom `--env-dir`, use that directory's interpreter instead.
With the environment active, verify packages, downloaded inputs and inference:

```bash
python -m pip check
python utils/check_environment.py --check-data
```

The demo runs exercises 11 and 14 using the bundled pretrained immunogenicity
checkpoint and downloaded MHCflurry models. It does not train a model.
Exercise 14 can also be run directly:

```bash
python solutions/ex14_prioritization.py
```

To inspect topology using the bundled DeepTMHMM output:

```bash
python solutions/ex11_topology_reference.py --deeptmhmm data/precomputed/deeptmhmm_fhaC.3line
```

The training exercises are optional and can be run independently:

```bash
python solutions/ex12_immunogenicity.py
python solutions/ex13_nettcr.py
```

Exercise 12 overwrites `data/precomputed/ex12_immunogenicity_cnn.pt`; later
runs of exercise 14 and the demo use the newly trained weights. Save a copy of
that checkpoint first if you want to retain the bundled model.

All generated tables, figures and logs are written to `results/`, replacing
files with the same names. The bundled outputs are examples; prediction
scores and shortlist sizes can vary with the model weights and environment.

## Setup options

```bash
bash setup.sh --help
```

| Option | Effect |
|---|---|
| `--env-dir PATH` | Create or reuse a venv at `PATH`; default `.venv` |
| `--use-current-env` | Install into the active venv or Conda environment |
| `--skip-models` | Skip MHCflurry weight download and pretrained-model checks |
| `--skip-data` | Skip public input download and input checks |

The two environment options are mutually exclusive. Skipping downloads does
not make missing inputs optional for exercises that require them. To install
packages only:

```bash
bash setup.sh --skip-models --skip-data
source .venv/bin/activate
python utils/check_environment.py --packages-only
```

Rerun `bash setup.sh` without skip flags to complete the installation. In the
Conda workflow, retain `--use-current-env` when rerunning setup. MHCflurry keeps
its weights in its user data cache; inspect the location with:

```bash
mhcflurry-downloads info
```

## External predictors

Install each external predictor according to its own release instructions
and licence. Keep incompatible dependencies in separate environments. The
commands below configure existing installations; they do not install them.

### NetMHCpan and NetMHCIIpan

Set the paths to the executables and run the wrappers:

```bash
export NETMHCPAN_BIN=/path/to/netMHCpan
export NETMHCIIPAN_BIN=/path/to/netMHCIIpan
python solutions/ex1_netmhcpan.py data/fhaC.fasta
python solutions/ex2_netmhciipan.py data/fhaC.fasta
```

The optional IEDB route uses a remote service instead of local executables:

```bash
python solutions/ex15_netmhcpan_iedb.py
```

This requires internet access and service availability. It is not needed for
the core demo.

### MixMHCpred and PRIME

The wrappers target MixMHCpred 2.2 and PRIME 2.0. With those releases installed:

```bash
export MIXMHCPRED=/path/to/MixMHCpred
export PRIME=/path/to/PRIME
python solutions/ex4_mixmhcpred.py data/fhaC.fasta
python solutions/ex6_prime.py data/fhaC.fasta
```

### SignalP and DeepTMHMM

Run the SignalP wrapper with the `signalp6` command on `PATH`, for example
from its own environment:

```bash
python solutions/ex7_signalp6.py data/fhaC.fasta other
```

To analyse the bundled DeepTMHMM output, no DeepTMHMM installation is needed:

```bash
python solutions/ex8_deeptmhmm.py data/fhaC.fasta data/precomputed/deeptmhmm_fhaC.3line
```

If the `.3line` argument is omitted, exercise 8 invokes `biolib run
DTU/DeepTMHMM`; that route needs a configured BioLib installation and access
to its service. Outputs for a different protein must be generated for that
protein; the bundled FhaC prediction cannot be reused.

After generating new topology predictions, return to the tutorial environment
and provide their actual paths to exercise 11:

```bash
python solutions/ex11_topology_reference.py --signalp results/signalp6/prediction_results.txt --deeptmhmm /path/to/predicted_topologies.3line
```

Both flags are optional. Without them, exercise 11 reports the UniProt
reference alone.

### Original DeepImmuno and NetTCR packages

Exercise 9 writes an input file and prints a command for the original
DeepImmuno package:

```bash
python solutions/ex9_deepimmuno.py data/fhaC.fasta
```

Install that package in a separate environment and run the printed command
from its repository root. Exercise 9 does not execute DeepImmuno itself.

Exercise 10 combines `results/prime_class_I.csv` with supplied TCR sequences:

```bash
python solutions/ex10_nettcr_prepare.py data/tcrs_example.csv 20
```

Run exercise 6 first to create the PRIME table. This step only prepares
peptide/CDR3 pairs; check the input requirements of your NetTCR release before
using them. The bundled TCR examples are not experimental evidence for FhaC.
See [data/README_tcrs.md](data/README_tcrs.md).

## Colab

Open `notebooks/immunoinformatics_colab.ipynb` in Google Colab, select a
Python 3.10-3.12 runtime and run the cells in order. The notebook obtains the
module by clone or ZIP upload, then installs its dependencies remotely with
`setup.sh`. Local package installation is not required.

## Troubleshooting

| Problem | Action |
|---|---|
| `ModuleNotFoundError: No module named 'torch'` | Select `.venv/bin/python` or activate the tutorial Conda environment; if dependencies are missing there, rerun `bash setup.sh --use-current-env` |
| `externally-managed-environment` | Run `bash setup.sh` and activate `.venv` before using `python -m pip` |
| `ensurepip` or `pip` is missing | Install the venv package matching your Python, or use Conda; select a fresh `--env-dir` if creation was incomplete |
| Unsupported Python / no compatible wheel | Use Python 3.10-3.12; the Conda specification selects 3.11 |
| Missing models or public inputs | Rerun setup without skip flags, keeping `--use-current-env` for Conda |
| `signalp6`, `biolib` or another external command is missing | Activate the predictor's environment or configure its executable path |
| DeepImmuno cannot find `data/` or `models/` | Run the original package from its own repository root |

Confirm the selected interpreter when diagnosing environment problems:

```bash
python -c "import sys; print(sys.executable)"
```

See [README.md](README.md) for the module overview and
[RUNBOOK.md](RUNBOOK.md) for exercise inputs and outputs.
