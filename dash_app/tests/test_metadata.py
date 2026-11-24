"""Tests for the simulation metadata (metadata.py + data.py's run table).

These pin the encodings in dash_app/metadata/ (decided in the metadata
review, ismip6_check_261123/metadata_review/) so a later edit can't
silently change them, and check that the metadata files and the bundled
data agree.

Run from dash_app/:  python -m pytest
"""

import os

import numpy as np
import pandas as pd
import pytest

import analysis as A
import metadata as M
from data import BIAS_ELIGIBLE, PUBLICATION_LABEL, RUNS

AIS = RUNS[RUNS["publication"] == PUBLICATION_LABEL["AIS"]]
GIS = RUNS[RUNS["publication"] == PUBLICATION_LABEL["GIS"]]


def _by_config(col, runs=None):
    return (AIS if runs is None else runs).groupby(["group", "model"])[col].first().to_dict()


# ── structure ─────────────────────────────────────────────────────────────

def test_every_run_has_an_exact_metadata_row():
    for ice, group, model in RUNS[["ice_sheet", "group", "model"]].drop_duplicates().itertuples(index=False):
        M.model_meta(ice, group, model)  # raises KeyError if missing


def test_no_fallback_for_unknown_keys():
    with pytest.raises(KeyError):
        M.model_meta("AIS", "ULB", "FETISH1")  # a similar-looking but nonexistent key
    with pytest.raises(KeyError):
        M.model_meta("GIS", "ULB", "fETISh_16km")  # AIS model never borrowed for GIS


def test_ais_and_gis_rows_are_independent():
    # These (group, model) pairs exist for both ice sheets; their rows must differ.
    for group, model in [("IMAU", "IMAUICE1"), ("IMAU", "IMAUICE2"), ("LSCE", "GRISLI"), ("NCAR", "CISM"), ("VUW", "PISM")]:
        a, g = M.model_meta("AIS", group, model), M.model_meta("GIS", group, model)
        assert a is not g and a["study"] == "Seroussi2020" and g["study"] == "Goelzer2020"
    # e.g. IMAU-ICE: plastic Mohr-Coulomb (SIA+SSA) for AIS, Weertman (SIA-only) for GIS
    assert M.model_meta("AIS", "IMAU", "IMAUICE1")["sliding_law"] == "Coulomb-limited"
    assert M.model_meta("GIS", "IMAU", "IMAUICE1")["sliding_law"] == "Power law (Weertman)"


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
              ("VUB", "AISMPALEO")]:
        assert init[k] == spin, k
    for k in [("ILTS_PIK", "SICOPOLIS"), ("LSCE", "GRISLI"), ("NCAR", "CISM"), ("ULB", "fETISh_16km"), ("ULB", "fETISh_32km"),
              ("VUW", "PISM")]:
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
    unrev = RUNS[RUNS["publication"] == "Rahlves 2025"]
    assert unrev["initialization"].str.endswith(M.UNREVIEWED_SUFFIX).all()
    assert (unrev["stress_balance"] == M.NOT_REVIEWED).all()
    reviewed = set().union(*M.CATEGORY_ORDER.values())
    assert not (set(unrev["initialization"]) | set(unrev["sliding_law"])) & reviewed


# ── weighting: every (group, model) configuration counts on its own ───────

def test_external_study_weighs_one_ismip6_configuration():
    valid = (RUNS["ice_sheet"] == "AIS").to_numpy() & RUNS["publication"].isin(
        [PUBLICATION_LABEL["AIS"], "Coulon 2024"]).to_numpy()
    w = A.run_weights(valid)
    coulon = (RUNS["publication"] == "Coulon 2024").to_numpy()
    assert np.isclose(w[coulon].sum(), 188 / 16)


# ── the metadata files ────────────────────────────────────────────────────

def test_every_study_has_a_models_and_an_experiments_file():
    for stem in M.STUDIES:
        for kind in ("models", "experiments"):
            assert os.path.isfile(os.path.join(M.META_DIR, f"{stem}_{kind}.csv")), (stem, kind)


def test_seroussi_experiments_are_verbatim_table_1():
    e = pd.read_csv(os.path.join(M.META_DIR, "Seroussi2020_experiments.csv"), comment="#", dtype=str)
    assert list(e.columns) == ["Experiment", "AOGCM", "Scenario", "Ocean forcing", "Ocean sensitivity",
                               "Ice shelf fracture", "Tier"]
    assert set(e["Ocean sensitivity"]) == {"Medium", "Medium*", "High*", "Low*", "PIGL"}


def test_every_bundled_run_has_an_experiment_row():
    for ice, group, model, exp in RUNS[["ice_sheet", "group", "model", "exp"]].drop_duplicates().itertuples(index=False):
        M.run_meta(ice, group, model, exp)  # raises KeyError if missing


@pytest.mark.parametrize("fn, stem, prefix", [
    ("external_sources_rahlves2025_gis.csv.gz", "Rahlves2025", ""),
    ("external_sources_coulon2024_ais.csv.gz", "Coulon2024", ""),
    ("external_sources_aschwanden2022_gis.csv.gz", "Aschwanden2019", ""),
    ("external_sources_goelzer2025_gis.csv.gz", "Goelzer2025", ""),
    ("external_sources_deconto2021_ais.csv.gz", "DeConto2021", ""),
    ("external_sources_edwards2021_ais.csv.gz", "Edwards2021", ""),
    ("external_sources_edwards2021_ais_risk.csv.gz", "Edwards2021", "RISK_"),
    ("external_sources_edwards2021_gis.csv.gz", "Edwards2021", ""),
    ("external_sources_ar6_ais.csv.gz", "IPCCAR6", ""),
    ("external_sources_ar6_gis.csv.gz", "IPCCAR6", ""),
])
def test_bundled_data_agrees_with_experiments_file(fn, stem, prefix):
    """The bundled data still carries climate_model/scenario columns written
    by utilities/external_sources.py; if the data are regenerated they must
    stay consistent with the study's experiments file."""
    d = pd.read_csv(os.path.join(os.path.dirname(M.META_DIR), "data", fn), dtype=str, keep_default_na=False,
                    usecols=["IS", "Exp", "climate_model", "scenario"]).drop_duplicates("Exp")
    for ice, exp, gcm, scen in d[["IS", "Exp", "climate_model", "scenario"]].itertuples(index=False):
        e = M.EXP_META[(stem, ice, prefix + exp)]
        assert (e["climate_model"], e["scenario"]) == (gcm, scen), (stem, exp)


# ── display order of ordinal characteristics ─────────────────────────────

def test_bias_plots_order_resolution_finest_first():
    import figures as F
    valid = (RUNS["publication"] == PUBLICATION_LABEL["AIS"]).to_numpy()
    b = A.bias(valid, 2015, 2023)
    w = A.run_weights(valid)
    fig = F.bias_figure(b, w, valid, "resolution", "Gt/yr")
    assert list(fig.layout.yaxis.ticktext) == M.CATEGORY_ORDER["resolution"]  # drawn top-down
    rows, _ = A.low_bias_enrichment(b[valid], w[valid], RUNS["resolution"].to_numpy()[valid])
    bars = F.enrichment_figure(rows, "resolution", 0.10)
    assert list(bars.data[0].y) == M.CATEGORY_ORDER["resolution"][::-1]  # drawn bottom-up


# ── ISMIP6-GIS (Goelzer et al. 2020), reviewed 2026-11-23 ─────────────────

def test_all_165_gis_runs_are_present():
    assert len(GIS) == 165
    assert GIS[["group", "model"]].drop_duplicates().shape[0] == 21


def test_goelzer_experiments_are_verbatim_table_1():
    e = pd.read_csv(os.path.join(M.META_DIR, "Goelzer2020_experiments.csv"), comment="#", dtype=str)
    assert list(e.columns) == ["Exp ID", "GCM", "RCP", "Ocean sensitivity", "Table 3 group"]
    assert e.set_index("Exp ID").loc["exp06", "GCM"] == "NorESM-M"  # verbatim; encoded NorESM1-M


def test_gis_scenarios_and_sensitivity():
    sc = GIS.groupby("exp")["scenario"].first()
    assert sc["exp07"] == "RCP2.6" and set(sc.drop("exp07")) == {"RCP8.5"}
    sens = GIS.groupby("exp")["ocean_sensitivity"].first()
    assert sens["exp09"] == "High" and sens["exp10"] == "Low" and set(sens.drop(["exp09", "exp10"])) == {"Medium"}


def test_gcm_names_are_shared_across_ice_sheets():
    names = set(AIS["climate_model"]) | set(GIS["climate_model"])
    assert not names & {"NorESM", "NorESM-M", "CSIRO-MK3", "CSIRO-Mk3.6", "ACCESS1.3"}
    assert {"NorESM1-M", "HadGEM2-ES", "IPSL-CM5A-MR", "CSIRO-Mk3-6-0"} <= set(AIS["climate_model"]) & set(GIS["climate_model"])


def test_gis_open_format_configurations():
    forcing = GIS.groupby(["group", "model"])["ocean_forcing"].agg(lambda s: set(s))
    open_ = {("BGC", "BISICLES"), ("UAF", "PISM2"), ("UCIJPL", "ISSM2"), ("VUW", "PISM")}
    for k, v in forcing.items():
        assert v == ({"Open"} if k in open_ else {"Standard"}), k


def test_gis_model_classes():
    sl, init, sb = _by_config("sliding_law", GIS), _by_config("initialization", GIS), _by_config("stress_balance", GIS)
    issm = [("AWI", "ISSM1"), ("AWI", "ISSM2"), ("AWI", "ISSM3"), ("GSFC", "ISSM"), ("JPL", "ISSM"), ("JPL", "ISSMPALEO"),
            ("UCIJPL", "ISSM1"), ("UCIJPL", "ISSM2")]
    assert all(sl[k] == "Budd" for k in issm)
    assert sl[("UAF", "PISM1")] == sl[("UAF", "PISM2")] == sl[("MUN", "GSM1")] == sl[("NCAR", "CISM")] == "Coulomb-limited"
    assert sl[("MUN", "GSM2")] == sl[("BGC", "BISICLES")] == sl[("LSCE", "GRISLI")] == "Linear"
    assert sl[("VUB", "GISM")] == sl[("ILTS_PIK", "SICOPOLIS1")] == "Budd"
    assert {k for k, v in sb.items() if v == "Higher-order"} == {
        ("AWI", "ISSM1"), ("AWI", "ISSM2"), ("AWI", "ISSM3"), ("BGC", "BISICLES"), ("JPL", "ISSM"), ("JPL", "ISSMPALEO"),
        ("NCAR", "CISM"), ("UCIJPL", "ISSM1"), ("UCIJPL", "ISSM2"), ("VUB", "GISM")}
    assert {k for k, v in init.items() if v == "Spin-up / equilibrium"} == {
        ("IMAU", "IMAUICE1"), ("IMAU", "IMAUICE2"), ("MUN", "GSM1"), ("MUN", "GSM2")}
    assert init[("VUW", "PISM")] == init[("UAF", "PISM1")] == "Hybrid (thickness-nudged)"
    assert init[("JPL", "ISSMPALEO")] == "Data assimilation (velocity inversion)"
