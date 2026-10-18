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

In a notebook, retrieve the secret with:

```python
from google.colab import userdata

api_key = userdata.get("ALPHA_GENOME_API_KEY")
