"""
Numerical layer: per-window rates, run weights, KDEs, time-series bands,
2100 values, and weighted ANOVA. Pure numpy/scipy (no statsmodels -- keeps
the Render free-tier memory budget intact); checked against the notebook's
statsmodels results (see the verification notes in the redesign plan).
"""

from functools import lru_cache

import numpy as np
import scipy.stats
from scipy.stats import linregress

from data import (
    CUM, DEFAULT_OBS, ICE_SHEET_AREA_M2, IMBIE, IS_ISMIP6, MASS_COL, RUN_VALS, RUN_YEARS, RUNS, YEAR_GRID, proj_start,
)

_ICE = RUNS["ice_sheet"].to_numpy()
_PUB = RUNS["publication"].to_numpy()
_GROUP = RUNS["group"].to_numpy()


def _fast_slope(x, y):
    """Closed-form OLS slope (== scipy.stats.linregress(x, y).slope)."""
    if len(x) < 2:
        return np.nan
    xd = x - x.mean()
    ssxx = np.dot(xd, xd)
    return np.nan if ssxx == 0 else np.dot(xd, y - y.mean()) / ssxx


@lru_cache(maxsize=64)
def run_rates(lo, hi):
    """Per-run rate (Gt/yr) over [lo, hi] inclusive -- the same window
    semantics as the original app's _rows_for_window. NaN where a run has
    fewer than 2 points in the window."""
    out = np.full(len(RUNS), np.nan)
    for i, (y, v) in enumerate(zip(RUN_YEARS, RUN_VALS)):
        m = (y >= lo) & (y <= hi)
        if m.sum() >= 2:
            out[i] = _fast_slope(y[m], v[m])
    out.setflags(write=False)
    return out


@lru_cache(maxsize=128)
def imbie_rate(ice_sheet, lo, hi, obs=DEFAULT_OBS):
    """(slope, stderr) of IMBIE cumulative mass over Jan lo .. Dec hi, from
    observation product `obs` (a data.OBS_PRODUCTS key)."""
    df = IMBIE[obs][ice_sheet]
    m = (df["Year"] >= lo) & (df["Year"] < hi + 1)
    r = linregress(df.loc[m, "Year"], df.loc[m, MASS_COL])
    return r.slope, r.stderr


def checked_mask(checked):
    return np.isin(_PUB, list(checked))


def run_weights(valid):
    """Weight per run (0 where not valid). ISMIP6 runs count 1 each; every
    run of another study counts typical_model_n / n_study, so a whole study
    weighs about as much as one typical ISMIP6 institution -- the same rule
    the original app's KDEs used. Without ISMIP6 in a panel, everything is 1."""
    w = np.zeros(len(RUNS))
    for ice_sheet in ("AIS", "GIS"):
        m = valid & (_ICE == ice_sheet)
        ism = m & IS_ISMIP6
        n_ism = ism.sum()
        typical = max(1.0, n_ism / max(1, len(set(_GROUP[ism])))) if n_ism else None
        w[ism] = 1.0
        for pub in set(_PUB[m & ~IS_ISMIP6]):
            sel = m & (_PUB == pub)
            w[sel] = typical / sel.sum() if typical else 1.0
    return w


def weighted_quantiles(values, weights, qs):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ok = np.isfinite(values) & (weights > 0)
    values, weights = values[ok], weights[ok]
    if len(values) == 0:
        return np.full(len(qs), np.nan)
    if len(values) == 1:
        return np.full(len(qs), values[0])
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
    return np.interp(qs, cdf, values)


def kde_curve(values, weights, xs):
    """Weighted Gaussian KDE normalized to a peak of 1 (each curve is
    self-normalized -- see the user's 2026-10-01 decision); flat if
    degenerate."""
    if len(values) < 2 or np.std(values) == 0:
        return np.zeros_like(xs)
    try:
        dens = scipy.stats.gaussian_kde(values, weights=weights)(xs)
    except np.linalg.LinAlgError:
        return np.zeros_like(xs)
    peak = dens.max()
    return dens / peak if peak > 0 else dens


def category_masks(valid, ice_sheet, dim):
    """[(category, bool mask)] for runs valid in this panel; one pooled
    entry when dim is None."""
    m = valid & (_ICE == ice_sheet)
    if dim is None:
        return [("All simulations", m)] if m.any() else []
    col = RUNS[dim].to_numpy()
    return [(c, m & (col == c)) for c in sorted(set(col[m]))]


# ── time series ───────────────────────────────────────────────────────────

TS_PAD_YEARS = 6
TS_QUANTILES = (0.05, 0.25, 0.5, 0.75, 0.95)


def timeseries(valid, ice_sheet, dim, lo, hi):
    """Per category: years and weighted 5/25/50/75/95 % of cumulative change
    rebased to 0 at `lo`, over [lo - 6, hi + 6]. A run with no value at
    `lo` but one at `lo + 1` (Edwards 2021 starts in 2016) is rebased there
    instead -- the same one-year fallback change_2015_2100 uses; runs with
    neither are left out."""
    years = np.arange(lo - TS_PAD_YEARS, hi + TS_PAD_YEARS + 1)
    cols = years - YEAR_GRID[0]
    base = CUM[:, lo - YEAR_GRID[0]].copy()
    missing = ~np.isfinite(base)
    base[missing] = CUM[missing, lo + 1 - YEAR_GRID[0]]
    ok = valid & np.isfinite(base)
    w = run_weights(ok)
    out = []
    for cat, m in category_masks(ok, ice_sheet, dim):
        sub = CUM[m][:, cols] - base[m][:, None]
        q = np.full((len(TS_QUANTILES), len(years)), np.nan)
        n = np.isfinite(sub).sum(axis=0)
        for j in range(len(years)):
            if n[j]:
                q[:, j] = weighted_quantiles(sub[:, j], w[m], TS_QUANTILES)
        out.append({"category": cat, "years": years, "q": q, "n_runs": int(m.sum())})
    return out


def imbie_timeseries(ice_sheet, lo, hi, obs=DEFAULT_OBS):
    """IMBIE cumulative change rebased to 0 at Jan `lo`, with 1-sigma
    uncertainty of that change. IMBIE accumulates its cumulative uncertainty
    in quadrature (checked: the file's end value matches annual errors summed
    in quadrature), so the uncertainty of change since `lo` is
    sqrt(|U(t)^2 - U(lo)^2|)."""
    full = IMBIE[obs][ice_sheet]
    df = full[(full["Year"] >= lo - TS_PAD_YEARS) & (full["Year"] < hi + TS_PAD_YEARS + 1)]
    at0 = full.loc[full["Year"] == lo]
    if at0.empty:
        return None
    u = df["Cumulative ice sheet mass change uncertainty (Gt)"].to_numpy()
    u0 = at0["Cumulative ice sheet mass change uncertainty (Gt)"].iloc[0]
    return {
        "years": df["Year"].to_numpy(),
        "value": df[MASS_COL].to_numpy() - at0[MASS_COL].iloc[0],
        "sigma": np.sqrt(np.abs(u ** 2 - u0 ** 2)),
    }


# ── 2100 ──────────────────────────────────────────────────────────────────

END_YEAR = 2100


def change_2015_2100():
    """Per-run change from 2015 to 2100 (Gt), NaN where unavailable.
    Edwards 2021 starts in 2016 (its series is zeroed there), so a run
    without a 2015 value falls back to 2016 as its baseline."""
    end = CUM[:, END_YEAR - YEAR_GRID[0]]
    base = CUM[:, proj_start - YEAR_GRID[0]].copy()
    fallback = CUM[:, proj_start + 1 - YEAR_GRID[0]]
    base[~np.isfinite(base)] = fallback[~np.isfinite(base)]
    return (end - base).astype(float)


CHANGE_2100 = change_2015_2100()


# ── bias & ANOVA ──────────────────────────────────────────────────────────

def bias(valid, lo, hi, obs=DEFAULT_OBS, area_normalized=False):
    """Run rate minus the IMBIE rate for the run's ice sheet (Gt/yr), or
    per unit area (mm/yr water equivalent: 1 Gt/yr / area[m^2] * 1e12 ==
    mm/yr w.e.) so AIS and GIS can be pooled."""
    rates = run_rates(lo, hi)
    b = np.full(len(RUNS), np.nan)
    for ice_sheet in ("AIS", "GIS"):
        m = valid & (_ICE == ice_sheet)
        b[m] = rates[m] - imbie_rate(ice_sheet, lo, hi, obs)[0]
        if area_normalized:
            b[m] *= 1e12 / ICE_SHEET_AREA_M2[ice_sheet]
    return b


def _design(cols, levels):
    """Treatment-coded dummies for each categorical column."""
    blocks = []
    for col, lv in zip(cols, levels):
        blocks.append(np.column_stack([(col == l).astype(float) for l in lv[1:]]) if len(lv) > 1 else
                      np.zeros((len(col), 0)))
    return blocks


def _wls_rss(X, y, sw):
    Xw, yw = X * sw[:, None], y * sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    r = yw - Xw @ beta
    return float(r @ r), int(np.linalg.matrix_rank(Xw))


def anova(y, w, factors):
    """Weighted (WLS) ANOVA of y on categorical factors.

    Returns (one_way, combined):
      one_way  -- per factor, fit alone: R^2, F, p, n categories (sorted by R^2)
      combined -- all factors together, Type II (drop-one-term) SS, df, F, p,
                  plus the model's R^2. With main effects only, Type II equals
                  the drop-one-term comparison, as in statsmodels' anova_lm.
    """
    ok = np.isfinite(y) & (w > 0)
    y, w = y[ok], w[ok]
    factors = {k: np.asarray(v)[ok] for k, v in factors.items()}
    factors = {k: v for k, v in factors.items() if len(set(v)) >= 2}
    n = len(y)
    sw = np.sqrt(w)
    ones = np.ones((n, 1))
    ybar = np.average(y, weights=w)
    tss = float(np.sum(w * (y - ybar) ** 2))

    one_way = []
    for k, col in factors.items():
        lv = sorted(set(col))
        X = np.hstack([ones] + _design([col], [lv]))
        rss, rank = _wls_rss(X, y, sw)
        df_m, df_r = rank - 1, n - rank
        if df_m < 1 or df_r < 1 or tss == 0:
            continue
        f = ((tss - rss) / df_m) / (rss / df_r)
        one_way.append({"characteristic": k, "r2": 1 - rss / tss, "F": f,
                        "p": float(scipy.stats.f.sf(f, df_m, df_r)), "n_categories": len(lv), "n": n})
    one_way.sort(key=lambda r: -r["r2"])

    combined = {"terms": [], "r2": np.nan, "n": n}
    if factors and tss > 0:
        keys = list(factors)
        levels = [sorted(set(factors[k])) for k in keys]
        blocks = _design([factors[k] for k in keys], levels)
        X_full = np.hstack([ones] + blocks)
        rss_full, rank_full = _wls_rss(X_full, y, sw)
        df_r = n - rank_full
        combined["r2"] = 1 - rss_full / tss
        for i, k in enumerate(keys):
            X_red = np.hstack([ones] + [b for j, b in enumerate(blocks) if j != i])
            rss_red, rank_red = _wls_rss(X_red, y, sw)
            df_t = rank_full - rank_red
            if df_t < 1 or df_r < 1:
                combined["terms"].append({"characteristic": k, "SS": 0.0, "df": 0, "F": np.nan, "p": np.nan})
                continue
            ss = rss_red - rss_full
            f = (ss / df_t) / (rss_full / df_r)
            combined["terms"].append({"characteristic": k, "SS": ss, "df": df_t, "F": f,
                                      "p": float(scipy.stats.f.sf(f, df_t, df_r))})
    return one_way, combined


def low_bias_enrichment(b, w, col, frac=0.10):
    """For the `frac` lowest-|bias| runs (by weight), each category's
    weighted share among them vs. its weighted share overall."""
    ok = np.isfinite(b) & (w > 0)
    b, w, col = np.abs(b[ok]), w[ok], np.asarray(col)[ok]
    if len(b) == 0:
        return [], np.nan
    cut = weighted_quantiles(b, w, [frac])[0]
    top = b <= cut
    rows = []
    for c in sorted(set(col)):
        m = col == c
        rows.append({"category": c,
                     "share_low_bias": w[m & top].sum() / w[top].sum(),
                     "share_all": w[m].sum() / w.sum()})
    return rows, cut
