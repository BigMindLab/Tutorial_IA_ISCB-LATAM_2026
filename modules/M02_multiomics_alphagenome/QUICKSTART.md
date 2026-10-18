# AlphaGenome Quickstart

This guide explains how to request an AlphaGenome API key and make it available
to the M02 notebooks without adding it to the repository.

The AlphaGenome API is available for non-commercial use, subject to the
[AlphaGenome Terms of Use](https://deepmind.google.com/science/alphagenome/terms).

## 1. Request an API key

1. Visit the [AlphaGenome API page](https://deepmind.google.com/science/alphagenome/api).
2. Sign in and follow the instructions to request API access.
3. Keep the key private. Do not paste it into a notebook, commit it to GitHub,
   or share it in an issue or pull request.

## 2. Add the key to Google Colab

1. Open the notebook in Colab.
2. Select **Secrets** (the key icon) in the left sidebar.
3. Add a secret named `ALPHA_GENOME_API_KEY` and paste your key as its value.
4. Enable notebook access for this secret.

## 3. Add the key to GitHub Codespaces

1. Open your GitHub account **Settings**.
2. Go to **Codespaces** secrets and select **New secret**.
3. Name the secret `ALPHA_GENOME_API_KEY`, enter your API key, and allow access
   to this repository.
4. If your Codespace is already running, stop and restart it so the new secret
   becomes available.

In Python, read the secret from the environment:

```python
import os

api_key = os.environ.get("ALPHA_GENOME_API_KEY")
In a notebook, retrieve the secret with:

```python
from google.colab import userdata

api_key = userdata.get("ALPHA_GENOME_API_KEY")
```

## 4. If you do not have an API key

You can review the notebooks and run cells that do not make AlphaGenome API
requests. Live predictions require a valid API key, so API-dependent cells
cannot return predictions without one. Request a key before the workshop if
you want to run those steps yourself.
