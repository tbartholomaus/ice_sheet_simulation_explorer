"""Tests for user_data.py.

Run from dash_app/:  ../.venv/bin/python -m pytest
"""

import base64
import os
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import linregress

import user_data as U

SAMPLE_CSV = os.path.join(os.path.dirname(__file__), "..", "test_user_ensembles", "test_user_ensemble.csv")


def _upload(years, vals):
    """Just the fields U.rates reads from a from_store() entry."""
    return SimpleNamespace(years=np.asarray(years, dtype=float), vals=np.asarray(vals, dtype=float))


def _linregress_rates(u, lo, hi):
    """Reference: scipy.stats.linregress on each simulation's non-missing
    values inside [lo, hi], one simulation at a time."""
    out = np.full(len(u.vals), np.nan)
    for i, v in enumerate(u.vals):
        m = (u.years >= lo) & (u.years <= hi) & np.isfinite(v)
        if m.sum() >= 2:
            out[i] = linregress(u.years[m], v[m]).slope
    return out


# ── rates ─────────────────────────────────────────────────────────────────

def test_rates_exact_slopes():
    years = np.arange(2010, 2021)
    u = _upload(years, [-50 * years, 3 + 0 * years, 12.5 * years])
    np.testing.assert_allclose(U.rates(u, 2015, 2020), [-50, 0, 12.5], atol=1e-9)


def test_rates_skip_missing_values():
    # Slope 10 everywhere it's defined; the gaps must be dropped, not read as 0.
    years = np.arange(2015, 2021, dtype=float)
    v = 10 * years
    v[[1, 3]] = np.nan
    np.testing.assert_allclose(U.rates(_upload(years, [v]), 2015, 2020), [10])


def test_rates_fewer_than_two_values_in_window_gives_nan():
    years = np.arange(2015, 2021, dtype=float)
    one_left = np.full(6, np.nan)
    one_left[2] = 5.0                       # 1 value in window
    none_left = np.full(6, np.nan)          # all missing
    two_left = np.full(6, np.nan)
    two_left[[0, 5]] = [0.0, 10.0]          # exactly 2: still fitted
    out = U.rates(_upload(years, [one_left, none_left, two_left]), 2015, 2020)
    assert np.isnan(out[0])
    assert np.isnan(out[1])
    assert out[2] == pytest.approx(2)


def test_rates_window_bounds_are_inclusive():
    u = _upload([2014, 2015, 2016, 2017], [[100, 0, 4, -100]])
    np.testing.assert_allclose(U.rates(u, 2015, 2016), [4])


def test_rates_match_linregress_with_gaps_and_sub_annual_steps():
    rng = np.random.default_rng(0)
    years = np.arange(1990, 2101, 0.5)
    vals = rng.normal(-200, 300, (300, len(years))).cumsum(axis=1)
    vals[rng.random(vals.shape) < 0.1] = np.nan
    u = _upload(years, vals)
    for lo, hi in [(2015, 2023), (2000, 2020), (2010, 2016)]:
        np.testing.assert_allclose(U.rates(u, lo, hi), _linregress_rates(u, lo, hi),
                                   rtol=1e-9, atol=1e-9, equal_nan=True)


def test_rates_match_linregress_on_sample_upload():
    # End to end through the same parsing the upload dialog uses.
    with open(SAMPLE_CSV, "rb") as f:
        contents = "data:text/csv;base64," + base64.b64encode(f.read()).decode()
    payload, _ = U.parse_upload(contents, "test_user_ensemble.csv")
    payload.update(label="sample", ice_sheet="AIS")
    u = U.from_store([payload])[0]
    np.testing.assert_allclose(U.rates(u, 2015, 2023), _linregress_rates(u, 2015, 2023),
                               rtol=1e-9, atol=1e-9, equal_nan=True)
