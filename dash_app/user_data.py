"""
A visitor's own uploaded ensembles: CSV parsing/validation and the same
numbers analysis.py computes for the bundled studies (rates, time-series
bands, 2015->2100 change).

Uploads are per visitor: the parsed data lives in the browser (a dcc.Store)
and comes back with each request, and is never merged into the shared run
table -- data.RUNS/CUM are module globals shared by every visitor of the
deployed app. Several uploads can be shown at once; each is always drawn as
its own group, so it is never pooled with bundled runs (or other uploads) and
its runs can simply be weighted equally.
"""

import base64
import binascii
import io
from types import SimpleNamespace

import numpy as np
import pandas as pd

from analysis import TS_END_YEAR, TS_QUANTILES, _fast_slope, weighted_quantiles
from data import YEAR_GRID, proj_start

MAX_BYTES = 5_000_000
MAX_RUNS = 1000
YEAR_RANGE = (1800, 2500)
MAX_ENSEMBLES = 5
# One color per upload slot, chosen to stand apart from the study colors.
USER_COLORS = ["#e41a1c", "#4daf4a", "#999933", "#9c755f", "#01665e"]  # >= 18 CIEDE2000 from every study color and each other
USER_COLOR = USER_COLORS[0]
DEFAULT_LABEL = "Your ensemble"


def parse_upload(contents, filename):
    """dcc.Upload contents -> (payload dict for the browser store, summary
    text). Raises ValueError with a message meant for the visitor."""
    if not contents or "," not in contents:
        raise ValueError("The file could not be read.")
    try:
        raw = base64.b64decode(contents.split(",", 1)[1], validate=False)
    except (binascii.Error, ValueError):
        raise ValueError("The file could not be read.")
    if len(raw) > MAX_BYTES:
        raise ValueError(f"The file is {len(raw) / 1e6:.1f} MB; the limit is {MAX_BYTES / 1e6:.0f} MB.")
    try:
        # sep=None sniffs commas, semicolons or tabs.
        df = pd.read_csv(io.BytesIO(raw), sep=None, engine="python")
    except Exception:  # pandas raises several parser/decoder types
        raise ValueError("The file isn't a readable text CSV.")
    if df.shape[1] < 2:
        raise ValueError("Expected a year column followed by at least one simulation column.")

    years = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    df = df.loc[years.notna()].copy()
    years = years[years.notna()].to_numpy(dtype=float)
    if len(years) < 2:
        raise ValueError(f"The first column ({df.columns[0]!r}) must hold years; fewer than 2 numeric years found.")
    if years.min() < YEAR_RANGE[0] or years.max() > YEAR_RANGE[1]:
        raise ValueError(f"Years must fall between {YEAR_RANGE[0]} and {YEAR_RANGE[1]}; "
                         f"found {years.min():g} to {years.max():g}. Is the first column the year?")
    if len(np.unique(years)) != len(years):
        raise ValueError("Some years appear more than once in the first column.")

    vals = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    n_bad = int((vals.isna() & df.iloc[:, 1:].notna()).to_numpy().sum())
    keep = vals.notna().sum() >= 2
    vals = vals.loc[:, keep]
    if vals.shape[1] == 0:
        raise ValueError("No simulation column has at least 2 numeric values.")
    if vals.shape[1] > MAX_RUNS:
        raise ValueError(f"Found {vals.shape[1]} simulations; the limit is {MAX_RUNS}.")

    order = np.argsort(years)
    years, vals = years[order], vals.to_numpy(dtype=float)[order]
    payload = {
        "filename": str(filename or "upload.csv"),
        "names": [str(c) for c in df.columns[1:][keep.to_numpy()]],
        "years": years.tolist(),
        # one list per simulation; NaN -> None so it survives JSON
        "values": [[None if not np.isfinite(v) else round(float(v), 4) for v in col] for col in vals.T],
    }
    notes = []
    if n_bad:
        notes.append(f"{n_bad} non-numeric cell{'s' if n_bad != 1 else ''} ignored")
    if (~keep).sum():
        notes.append(f"{int((~keep).sum())} column{'s' if (~keep).sum() != 1 else ''} without data skipped")
    summary = (f"Read {len(payload['names'])} simulation{'s' if len(payload['names']) != 1 else ''}, "
               f"{years.min():g}–{years.max():g}" + (f" ({'; '.join(notes)})" if notes else "") + ".")
    return payload, summary


def entries(store):
    """The store's list of uploads (a lone dict is a pre-multi-upload
    session's single upload)."""
    if not store:
        return []
    return [store] if isinstance(store, dict) else [e for e in store if e and e.get("values")]


def from_store(store):
    """Browser store -> list of namespaces of numpy arrays (possibly empty)."""
    out = []
    for i, e in enumerate(entries(store)):
        vals = np.array([[np.nan if v is None else v for v in col] for col in e["values"]], dtype=float)
        out.append(SimpleNamespace(label=e.get("label") or DEFAULT_LABEL, ice_sheet=e["ice_sheet"],
                                   color=e.get("color") or USER_COLORS[i % len(USER_COLORS)],
                                   names=e["names"], years=np.asarray(e["years"], dtype=float), vals=vals))
    return out


def unique_label(name, taken):
    """`name`, or `name (2)`, `name (3)`... if another upload already uses it
    (labels double as legend entries and category names)."""
    base = (name or "").strip() or DEFAULT_LABEL
    label, k = base, 2
    while label in taken:
        label, k = f"{base} ({k})", k + 1
    return label


def free_color(taken):
    return next((c for c in USER_COLORS if c not in taken), USER_COLORS[len(taken) % len(USER_COLORS)])


def rates(u, lo, hi):
    """Per-simulation rate (Gt/yr) over [lo, hi] inclusive -- the same
    window rule analysis.run_rates applies to the bundled runs."""
    out = np.full(len(u.vals), np.nan)
    inwin = (u.years >= lo) & (u.years <= hi)
    for i, v in enumerate(u.vals):
        m = inwin & np.isfinite(v)
        if m.sum() >= 2:
            out[i] = _fast_slope(u.years[m], v[m])
    return out


def _grid(u):
    """(n_sims x YEAR_GRID) matrix on integer years, like data.CUM."""
    mat = np.full((len(u.vals), len(YEAR_GRID)), np.nan)
    whole = u.years == np.round(u.years)
    idx = (u.years[whole] - YEAR_GRID[0]).astype(int)
    ok = (idx >= 0) & (idx < len(YEAR_GRID))
    mat[:, idx[ok]] = u.vals[:, whole][:, ok]
    return mat


def timeseries(u, lo, hi):
    """One analysis.timeseries-style entry for the upload (equal weights),
    or None when no simulation has a value at lo or lo + 1."""
    grid = _grid(u)
    years = np.arange(YEAR_GRID[0], TS_END_YEAR + 1)  # all of it, like analysis.timeseries
    base = grid[:, lo - YEAR_GRID[0]].copy()
    missing = ~np.isfinite(base)
    base[missing] = grid[missing, lo + 1 - YEAR_GRID[0]]
    ok = np.isfinite(base)
    if not ok.any():
        return None
    sub = grid[ok][:, years - YEAR_GRID[0]] - base[ok][:, None]
    q = np.full((len(TS_QUANTILES), len(years)), np.nan)
    for j in range(len(years)):
        col = sub[:, j]
        if np.isfinite(col).any():
            q[:, j] = weighted_quantiles(col, np.ones_like(col), TS_QUANTILES)
    return {"category": u.label, "years": years, "q": q, "n_runs": int(ok.sum()), "color": u.color}


def change_2015_2100(u):
    """Per-simulation 2015->2100 change (Gt), with the same 2016 fallback
    baseline analysis.change_2015_2100 uses; NaN where unavailable."""
    grid = _grid(u)
    end = grid[:, 2100 - YEAR_GRID[0]]
    base = grid[:, proj_start - YEAR_GRID[0]].copy()
    missing = ~np.isfinite(base)
    base[missing] = grid[missing, proj_start + 1 - YEAR_GRID[0]]
    return end - base
