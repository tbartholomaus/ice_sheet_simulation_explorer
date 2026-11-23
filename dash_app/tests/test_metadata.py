"""Tests for the simulation metadata (metadata.py + data.py's run table).

These pin the decisions of the metadata review (ismip6_check_261123/
metadata_review/01_Seroussi2020_ISMIP6_AIS.md) so a later edit can't
silently undo them.

Run from dash_app/:  python -m pytest
"""

import numpy as np
import pandas as pd
import pytest

import analysis as A
import metadata as M
from data import BIAS_ELIGIBLE, PUBLICATION_LABEL, RUNS

AIS = RUNS[RUNS["publication"] == PUBLICATION_LABEL["AIS"]]


def _by_config(col):
    return AIS.groupby(["group", "model"])[col].first().to_dict()


# ── structure ─────────────────────────────────────────────────────────────

def test_every_run_has_an_exact_metadata_row():
    for ice, group, model in RUNS[["ice_sheet", "group", "model"]].drop_duplicates().itertuples(index=False):
        M.model_meta(ice, group, model)  # raises KeyError if missing


def test_no_fallback_for_unknown_keys():
    with pytest.raises(KeyError):
        M.model_meta("AIS", "ULB", "FETISH1")  # the key the old substring fallback used to land on
    with pytest.raises(KeyError):
        M.model_meta("GIS", "ULB", "fETISh_16km")  # AIS model never borrowed for GIS


def test_ais_and_gis_rows_are_independent():
    # These (group, model) pairs exist for both ice sheets; their rows must differ.
    for group, model in [("IMAU", "IMAUICE1"), ("IMAU", "IMAUICE2"), ("LSCE", "GRISLI"), ("NCAR", "CISM"), ("VUW", "PISM")]:
        a, g = M.model_meta("AIS", group, model), M.model_meta("GIS", group, model)
        assert a is not g
        assert a["review_status"].startswith("reviewed") and g["review_status"] == "not yet reviewed"


def test_all_188_ais_runs_are_present():
    assert len(AIS) == 188
    assert AIS[["group", "model"]].drop_duplicates().shape[0] == 16


# ── ISMIP6-AIS experiments (Seroussi et al. 2020, Table 1) ────────────────

def test_tier2_experiments_are_cmip5_rcps():
    sc = AIS.groupby("exp")["scenario"].first()
    assert (sc[["expA1", "expA2", "expA3", "expA5", "expA6", "expA7"]] == "RCP8.5").all()
    assert (sc[["expA4", "expA8"]] == "RCP2.6").all()
    assert not AIS["scenario"].str.startswith("SSP").any()


def test_ocean_sensitivity_and_forcing():
    sens = AIS.groupby("exp")["ocean_sensitivity"].first()
    assert sens["exp09"] == "High" and sens["exp10"] == "Low" and sens["exp13"] == "PIGL"
    assert set(sens.drop(["exp09", "exp10", "exp13"])) == {"Medium"}
    forcing = AIS.groupby("exp")["ocean_forcing"].first()
    assert set(forcing[["exp01", "exp02", "exp03", "exp04", "exp11", "expA1", "expA2", "expA3", "expA4"]]) == {"Open"}
    assert set(forcing.drop(["exp01", "exp02", "exp03", "exp04", "exp11", "expA1", "expA2", "expA3", "expA4"])) == {"Standard"}


def test_noresm_is_noresm1_m():
    assert "NorESM" not in set(AIS["climate_model"])
    assert set(AIS.loc[AIS["exp"] == "exp05", "climate_model"]) == {"NorESM1-M"}


def test_ice_shelf_collapse_runs_are_excluded_from_bias():
    collapse = RUNS["ice_shelf_fracture"].to_numpy()
    assert collapse.sum() == 14
    assert set(RUNS.loc[collapse, "exp"]) == {"exp11", "exp12"}
    assert not BIAS_ELIGIBLE[collapse].any()
    assert BIAS_ELIGIBLE[~collapse].all()


# ── ISMIP6-AIS model classes (decisions 2026-11-23) ───────────────────────

def test_initialization_classes():
    init = _by_config("initialization")
    spin, hyb, da = "Spin-up / equilibrium", "Hybrid (thickness-nudged)", "Data assimilation (velocity inversion)"
    for k in [("AWI", "PISM1"), ("IMAU", "IMAUICE1"), ("IMAU", "IMAUICE2"), ("PIK", "PISM1"), ("PIK", "PISM2"),
              ("VUB", "AISMPALEO"), ("VUW", "PISM")]:
        assert init[k] == spin, k
    for k in [("ILTS_PIK", "SICOPOLIS"), ("LSCE", "GRISLI"), ("NCAR", "CISM"), ("ULB", "fETISh_16km"), ("ULB", "fETISh_32km")]:
        assert init[k] == hyb, k
    for k in [("DOE", "MALI"), ("JPL1", "ISSM"), ("UCIJPL", "ISSM"), ("UTAS", "ElmerIce")]:
        assert init[k] == da, k


def test_sliding_law_classes():
    sl = _by_config("sliding_law")
    expect = {
        ("AWI", "PISM1"): "Coulomb-limited", ("DOE", "MALI"): "Linear", ("ILTS_PIK", "SICOPOLIS"): "Budd",
        ("IMAU", "IMAUICE1"): "Coulomb-limited", ("IMAU", "IMAUICE2"): "Coulomb-limited", ("JPL1", "ISSM"): "Budd",
        ("LSCE", "GRISLI"): "Linear", ("NCAR", "CISM"): "Coulomb-limited", ("PIK", "PISM1"): "Coulomb-limited",
        ("PIK", "PISM2"): "Coulomb-limited", ("UCIJPL", "ISSM"): "Budd", ("ULB", "fETISh_16km"): "Power law (Weertman)",
        ("ULB", "fETISh_32km"): "Power law (Weertman)", ("UTAS", "ElmerIce"): "Linear", ("VUB", "AISMPALEO"): "Budd",
        ("VUW", "PISM"): "Coulomb-limited",
    }
    assert sl == expect


def test_stress_balance_and_resolution():
    sb = _by_config("stress_balance")
    assert {k for k, v in sb.items() if v == "Higher-order"} == {("DOE", "MALI"), ("NCAR", "CISM"), ("UCIJPL", "ISSM"),
                                                                 ("UTAS", "ElmerIce")}
    res = _by_config("resolution")
    assert res[("JPL1", "ISSM")] == "≤5 km" and res[("VUB", "AISMPALEO")] == ">10–20 km"
    assert res[("ULB", "fETISh_32km")] == ">20 km" and res[("PIK", "PISM1")] == ">5–10 km"
    assert set(res.values()) <= set(M.CATEGORY_ORDER["resolution"])


@pytest.mark.parametrize("r, expect", [(1, "≤5 km"), (5, "≤5 km"), (8, ">5–10 km"), (10, ">5–10 km"),
                                       (16, ">10–20 km"), (20, ">10–20 km"), (32, ">20 km")])
def test_resolution_bin_edges_belong_to_lower_bin(r, expect):
    assert M.resolution_bin(r) == expect


# ── unreviewed sources are labelled, never pooled with reviewed classes ───

def test_unreviewed_labels_are_marked():
    gis = RUNS[RUNS["publication"] == PUBLICATION_LABEL["GIS"]]
    assert gis["initialization"].str.endswith(M.UNREVIEWED_SUFFIX).all()
    assert (gis["stress_balance"] == M.NOT_REVIEWED).all()
    reviewed = set().union(*M.CATEGORY_ORDER.values())
    assert not (set(gis["initialization"]) | set(gis["sliding_law"])) & reviewed


# ── weighting: every (group, model) configuration counts on its own ───────

def test_external_study_weighs_one_ismip6_configuration():
    valid = (RUNS["ice_sheet"] == "AIS").to_numpy() & RUNS["publication"].isin(
        [PUBLICATION_LABEL["AIS"], "Coulon 2024"]).to_numpy()
    w = A.run_weights(valid)
    coulon = (RUNS["publication"] == "Coulon 2024").to_numpy()
    assert np.isclose(w[coulon].sum(), 188 / 16)
