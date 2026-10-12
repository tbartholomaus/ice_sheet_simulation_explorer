"""Tests for analysis.py.

Run from dash_app/:  ../.venv/bin/python -m pytest
"""

import numpy as np
import pytest
from scipy.stats import linregress

import analysis as A
from data import IMBIE, RUN_ID, RUN_VALS, RUN_YEARS, RUNS


# ── stacked_rates / run_rates ─────────────────────────────────────────────

# Gt/yr. Some runs have a true slope of exactly 0 (e.g. DeConto & Pollard
# 2016 values that rise and fall symmetrically over 2000-2020), where both
# methods return different ~1e-15 rounding noise that rtol alone rejects.
ATOL = 1e-9


def _rates(runs, lo, hi):
    """stacked_rates on [(years, vals), ...], stacked end to end the way data.py
    stacks the bundled runs."""
    run_id = np.concatenate([np.full(len(y), i) for i, (y, _) in enumerate(runs)])
    years = np.concatenate([np.asarray(y, dtype=float) for y, _ in runs])
    vals = np.concatenate([np.asarray(v, dtype=float) for _, v in runs])
    return A.stacked_rates(lo, hi, run_id, years, vals, len(runs))


def _linregress_slopes(runs, lo, hi):
    """Reference: scipy.stats.linregress on each run's windowed points, one
    run at a time (too slow for the app, ~4 s for the bundled runs)."""
    out = np.full(len(runs), np.nan)
    for i, (y, v) in enumerate(runs):
        y, v = np.asarray(y, dtype=float), np.asarray(v, dtype=float)
        m = (y >= lo) & (y <= hi)
        if m.sum() >= 2:
            out[i] = linregress(y[m], v[m]).slope
    return out


def test_stacked_rates_exact_slopes_with_different_time_steps():
    # Annual, half-yearly, and a run starting later -- each exactly linear.
    annual = np.arange(2010, 2021)
    half = np.arange(2010, 2020.5, 0.5)
    late = np.arange(2016, 2021)
    runs = [(annual, -100 * (annual - 2010)),
            (half, 3 + 25 * (half - 2010)),
            (late, -7.5 * late)]
    out = _rates(runs, 2012, 2020)
    np.testing.assert_allclose(out, [-100, 25, -7.5])


def test_stacked_rates_points_outside_window_are_ignored():
    # Slope 10 inside 2015-2020, but a jump outside it that would change a
    # fit over the whole run.
    y = np.arange(2000, 2031)
    v = np.where(y < 2015, 5000, 10 * y)
    v = np.where(y > 2020, -9999, v)
    out = _rates([(y, v)], 2015, 2020)
    np.testing.assert_allclose(out, [10])


def test_stacked_rates_window_bounds_are_inclusive():
    y = np.array([2015, 2016])
    out = _rates([(y, [0, 4])], 2015, 2016)
    np.testing.assert_allclose(out, [4])


def test_stacked_rates_fewer_than_two_points_gives_nan():
    runs = [(np.arange(2015, 2021), np.arange(6.0)),  # 6 points in window
            ([2010, 2018, 2030], [0, 1, 2]),          # 1 point in window
            ([1990, 1991], [0, 1]),                   # none in window
            ([2018, 2019], [0, 1])]                   # exactly 2: still fitted
    out = _rates(runs, 2015, 2020)
    assert out[0] == pytest.approx(1)
    assert np.isnan(out[1])
    assert np.isnan(out[2])
    assert out[3] == pytest.approx(1)


def test_stacked_rates_trailing_runs_without_points_still_get_an_entry():
    # bincount would otherwise size its output by the largest run id seen.
    runs = [(np.arange(2015, 2021), np.arange(6.0)), ([1990, 1991], [0, 1])]
    out = _rates(runs, 2015, 2020)
    assert out.shape == (2,)
    assert np.isnan(out[1])


def test_stacked_rates_matches_linregress_on_random_ragged_data():
    rng = np.random.default_rng(0)
    runs = []
    for _ in range(200):
        start = rng.integers(1990, 2020)
        step = rng.choice([0.5, 1.0, 10.0])
        y = np.arange(start, 2101, step)
        runs.append((y, rng.normal(0, 500, len(y)).cumsum()))
    for lo, hi in [(2015, 2023), (2000, 2020), (2010, 2016), (2019, 2024)]:
        np.testing.assert_allclose(_rates(runs, lo, hi), _linregress_slopes(runs, lo, hi),
                                   rtol=1e-9, atol=ATOL, equal_nan=True)


@pytest.mark.parametrize("lo, hi", [(2015, 2023), (2000, 2020), (2010, 2016)])
def test_run_rates_matches_linregress_on_bundled_data(lo, hi):
    bounds = np.flatnonzero(np.diff(RUN_ID)) + 1
    runs = list(zip(np.split(RUN_YEARS, bounds), np.split(RUN_VALS, bounds)))
    assert len(runs) == len(RUNS)
    out = A.run_rates(lo, hi)
    assert out.shape == (len(RUNS),)
    np.testing.assert_allclose(out, _linregress_slopes(runs, lo, hi), rtol=1e-9, atol=ATOL, equal_nan=True)



def test_run_rates_repeat_window_comes_from_cache():
    A.run_rates.cache_clear()
    first = A.run_rates(2012, 2021)
    assert A.run_rates(2012, 2021) is first
    assert A.run_rates.cache_info().hits == 1


def test_run_rates_cached_array_is_read_only():
    # Every request shares the cached array, so no caller may modify it.
    out = A.run_rates(2015, 2023)
    with pytest.raises(ValueError):
        out[0] = 0.0

# ── imbie_rate / imbie_timeseries ─────────────────────────────────────────

def test_imbie_rate_window_is_jan_lo_through_dec_hi(monkeypatch):
    # Monthly values, slope 12 Gt/yr inside 2015-2016 only. Jan 2017 is
    # outside the window (< hi + 1); Dec 2016 is inside.
    year = 2014 + np.arange(48) / 12
    mass = np.where((year >= 2015) & (year < 2017), 12 * year, 0.0)
    monkeypatch.setitem(IMBIE, "test", {"AIS": {"year": year, "mass": mass, "sigma": np.zeros_like(year)}})
    slope, stderr = A.imbie_rate("AIS", 2015, 2016, "test")
    assert slope == pytest.approx(12)
    assert stderr == pytest.approx(0, abs=1e-9)


def test_imbie_timeseries_rebases_to_window_start():
    d = IMBIE["imbie3"]["AIS"]
    ts = A.imbie_timeseries("AIS", 2010, 2023, "imbie3")
    i = np.flatnonzero(d["year"] == 2010)[0]
    assert ts["value"][i] == 0
    assert ts["sigma"][i] == 0
    np.testing.assert_allclose(ts["value"], d["mass"] - d["mass"][i])
    np.testing.assert_allclose(ts["sigma"], np.sqrt(np.abs(d["sigma"] ** 2 - d["sigma"][i] ** 2)))
    np.testing.assert_array_equal(ts["years"], d["year"])


def test_imbie_timeseries_none_when_start_year_not_in_record():
    assert A.imbie_timeseries("AIS", 1960, 2023, "imbie3") is None


def test_imbie_timeseries_does_not_modify_shared_record():
    d = IMBIE["imbie3"]["GIS"]
    before = {k: v.copy() for k, v in d.items()}
    A.imbie_timeseries("GIS", 2015, 2023, "imbie3")
    for k in before:
        np.testing.assert_array_equal(d[k], before[k])
