"""Shared AlphaGenome helpers for the M02 notebooks.

Every notebook gets its client and its results the same way:

    import sys
    sys.path.append("../utils")
    from ag_helpers import get_client, cached_call

    model = get_client()  # None means offline mode
    out = cached_call("ex1_tal1", lambda: model.predict_interval(...))

With a working API key, `cached_call` runs the request live. Without one, or
if the request fails, it loads the precomputed copy saved in
`data/precomputed/<name>/`. Set the environment variable `AG_PRECOMPUTE=1`
to save every live result there (this is how `scripts/precompute_M02.py`
builds the offline copies).

The API key is never printed or written to disk.
"""

import dataclasses
import getpass
import json
import os
import re
import shutil
from pathlib import Path
from typing import Any, Callable

import grpc
import numpy as np
import pandas as pd

from alphagenome.data import genome, junction_data, track_data
from alphagenome.models import dna_client, dna_output

KEY_NAME = "ALPHA_GENOME_API_KEY"
PRECOMPUTED = Path(__file__).resolve().parents[1] / "data" / "precomputed"
CONNECT_TIMEOUT = 20  # seconds; without it dna_client.create waits forever offline

_client = None
_client_ready = False  # True once get_client has chosen live or offline mode


# ---------------------------------------------------------------- API key and client


def get_api_key(prompt: bool = True) -> str | None:
    """Return the AlphaGenome API key, or None if there is none.

    Looks for a Colab secret, then an environment variable (Codespaces
    secret), then asks with a hidden prompt. Press Enter at the prompt to
    work offline.
    """
    key = None
    try:
        from google.colab import userdata

        key = userdata.get(KEY_NAME)
    except Exception:  # not in Colab, secret missing, or access not granted
        pass
    key = key or os.environ.get(KEY_NAME)
    if not key and prompt:
        try:
            key = getpass.getpass("AlphaGenome API key (press Enter for offline mode): ")
        except Exception:  # no interactive input, e.g. a batch run
            pass
    return (key or "").strip() or None


def get_client(prompt: bool = True):
    """Return an AlphaGenome client, or None for offline mode.

    The choice is remembered, so later calls (and `cached_call`) reuse it.
    Restart the kernel to pick up a key added later.
    """
    global _client, _client_ready
    if _client_ready:
        return _client
    key = get_api_key(prompt)
    if key is None:
        print("No API key found: offline mode, results load from data/precomputed/.")
    else:
        try:
            _client = dna_client.create(key, timeout=CONNECT_TIMEOUT)
            print("Live mode: requests go to the AlphaGenome API "
                  "(the key is checked on the first one).")
        except Exception as e:  # no network or server unreachable
            print(f"Could not reach the AlphaGenome API ({_reason(e)}): offline mode.")
    _client_ready = True
    return _client


def cached_call(name: str, fn: Callable[[], Any], save: bool | None = None) -> Any:
    """Run `fn()` live, or load the precomputed result called `name`.

    Live mode returns `fn()` and saves it when `save` is True (default: only
    when AG_PRECOMPUTE=1). If the API call fails and a precomputed copy
    exists, that copy is returned instead. Offline mode loads the copy.
    """
    path = _path(name)
    if get_client() is not None:
        try:
            result = fn()
        except grpc.RpcError as e:  # bad key, quota, network
            if not path.exists():
                raise
            print(f"API call failed ({_reason(e)}): using precomputed '{name}'.")
            return load_result(name)
        if save is None:
            save = os.environ.get("AG_PRECOMPUTE") == "1"
        if save:
            save_result(name, result)
            print(f"Saved '{name}' to data/precomputed/.")
        return result
    if not path.exists():
        raise FileNotFoundError(
            f"No precomputed result '{name}': this step needs an API key "
            "(see QUICKSTART.md)."
        )
    print(f"Offline: loaded precomputed '{name}'.")
    return load_result(name)


def _reason(e: Exception) -> str:
    code = getattr(e, "code", None)
    return code().name if isinstance(e, grpc.RpcError) and callable(code) else type(e).__name__


# ---------------------------------------------------------------- save and load
#
# Each result is a folder with meta.json plus data files:
#   TrackData / JunctionData   values.npz + metadata.parquet
#   Output                     one subfolder per output type (rna_seq/, dnase/, ...)
#   VariantOutput              reference/ and alternate/
#   DataFrame                  table.parquet
#   ndarray                    array.npz
#   dict / list                one subfolder per item
#   anything JSON can hold     stored in meta.json


def save_result(name: str, obj: Any) -> Path:
    """Write `obj` to data/precomputed/<name>/, replacing any old copy."""
    path = _path(name)
    tmp = path.with_name(path.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    _write(tmp, obj)
    shutil.rmtree(path, ignore_errors=True)
    tmp.rename(path)
    return path


def load_result(name: str) -> Any:
    """Read a result written by `save_result`."""
    return _read(_path(name))


def _path(name: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name.startswith("."):
        raise ValueError(f"Use letters, digits, '_', '-' or '.' in result names, not {name!r}.")
    return PRECOMPUTED / name


def _write(d: Path, obj: Any) -> None:
    d.mkdir(parents=True)
    meta: dict[str, Any] = {"type": type(obj).__name__}
    if isinstance(obj, (track_data.TrackData, junction_data.JunctionData)):
        np.savez_compressed(d / "values.npz", values=obj.values)
        _to_parquet(obj.metadata, d / "metadata.parquet")
        meta["interval"] = _interval_dict(obj.interval)
        meta["uns"] = obj.uns if _is_json(obj.uns) else None  # predictions leave it empty
        if isinstance(obj, track_data.TrackData):
            meta["resolution"] = obj.resolution
        else:
            meta["junctions"] = [str(j) for j in obj.junctions]
    elif isinstance(obj, dna_output.VariantOutput):
        _write(d / "reference", obj.reference)
        _write(d / "alternate", obj.alternate)
    elif isinstance(obj, dna_output.Output):
        meta["fields"] = [f.name for f in dataclasses.fields(obj)
                          if getattr(obj, f.name) is not None]
        for field in meta["fields"]:
            _write(d / field, getattr(obj, field))
    elif isinstance(obj, pd.DataFrame):
        _to_parquet(obj, d / "table.parquet")
    elif isinstance(obj, np.ndarray):
        np.savez_compressed(d / "array.npz", array=obj)
    elif isinstance(obj, (dict, list, tuple)) and not _is_json(obj):
        items = obj.items() if isinstance(obj, dict) else enumerate(obj)
        meta["type"] = "dict" if isinstance(obj, dict) else "list"
        meta["keys"] = []
        for i, (key, value) in enumerate(items):
            meta["keys"].append(key)
            _write(d / str(i), value)
    elif _is_json(obj):
        meta["type"], meta["value"] = "json", obj
    else:
        raise TypeError(
            f"Cannot save {type(obj).__name__}. Convert it first, e.g. "
            "variant_scorers.tidy_scores(...) for variant scores."
        )
    (d / "meta.json").write_text(json.dumps(meta, indent=1))


def _read(d: Path) -> Any:
    if not (d / "meta.json").exists():
        raise FileNotFoundError(f"No precomputed result in {d}.")
    meta = json.loads((d / "meta.json").read_text())
    kind = meta["type"]
    if kind in ("TrackData", "JunctionData"):
        with np.load(d / "values.npz") as f:
            values = f["values"]
        if values.dtype == np.float16:  # allowed on disk to save space
            values = values.astype(np.float32)
        common = dict(
            values=values,
            metadata=pd.read_parquet(d / "metadata.parquet"),
            interval=genome.Interval(**meta["interval"]) if meta["interval"] else None,
            uns=meta["uns"],
        )
        if kind == "TrackData":
            return track_data.TrackData(resolution=meta["resolution"], **common)
        junctions = np.array([genome.Interval.from_str(j) for j in meta["junctions"]])
        return junction_data.JunctionData(junctions=junctions, **common)
    if kind == "VariantOutput":
        return dna_output.VariantOutput(reference=_read(d / "reference"),
                                        alternate=_read(d / "alternate"))
    if kind == "Output":
        return dna_output.Output(**{field: _read(d / field) for field in meta["fields"]})
    if kind == "DataFrame":
        return pd.read_parquet(d / "table.parquet")
    if kind == "ndarray":
        with np.load(d / "array.npz") as f:
            return f["array"]
    if kind == "dict":
        return {key: _read(d / str(i)) for i, key in enumerate(meta["keys"])}
    if kind == "list":
        return [_read(d / str(i)) for i in range(len(meta["keys"]))]
    return meta["value"]


def _to_parquet(df: pd.DataFrame, path: Path) -> None:
    # Columns holding Interval, Variant or scorer objects are stored as text.
    df = df.copy()
    plain = (str, bool, int, float, np.number)
    for col in df.columns[df.dtypes == object]:
        if not df[col].map(lambda v: v is None or isinstance(v, plain)).all():
            df[col] = df[col].map(lambda v: None if v is None else str(v))
    df.to_parquet(path)


def _interval_dict(interval: genome.Interval | None) -> dict | None:
    if interval is None:
        return None
    return {f: getattr(interval, f) for f in ("chromosome", "start", "end", "strand", "name")}


def _is_json(obj: Any) -> bool:
    try:
        json.dumps(obj)
        return True
    except (TypeError, ValueError):
        return False
