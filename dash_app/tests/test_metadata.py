"""Tests for the simulation metadata (metadata.py + data.py's run table).

These pin the encodings in dash_app/metadata/ (decided in the metadata
review, docs/metadata_review/) so a later edit can't
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
    s20 = (RUNS["publication"] == PUBLICATION_LABEL["AIS"]).to_numpy()
    assert (collapse & s20).sum() == 14
    assert set(RUNS.loc[collapse & s20, "exp"]) == {"exp11", "exp12"}
    # Elsewhere only Seroussi 2024's expAE11-14 impose collapse (its own test)
    assert set(RUNS.loc[collapse, "publication"]) == {PUBLICATION_LABEL["AIS"], "Seroussi 2024 (ISMIP6 AIS 2300)"}
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
    """Rows not yet reviewed carry marked labels that never pool with reviewed
    classes; reviewed rows never carry them."""
    reviewed = set().union(*M.CATEGORY_ORDER.values())
    for key, row in M.MODEL_META.items():
        vals = {row[f] for f in M.CLASS_FIELDS if f != "ice_model"}
        if row["review_status"] == "not yet reviewed":
            assert not vals & reviewed, key
        elif row["review_status"] == "reviewed":
            assert vals <= reviewed, key


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
    # Table 1's columns, plus the regional climate model documented in the file header (Sect. 2.2)
    assert list(e.columns) == ["Experiment", "AOGCM", "Scenario", "Ocean forcing", "Ocean sensitivity",
                               "Ice shelf fracture", "Tier", "Regional climate model"]
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
    # The data keep the source's own spelling; translate it with the study's value maps.
    spec = M.COLUMNS[M.COLUMNS["study"] == stem]
    maps = {f: M._value_map(v) for f, v in zip(spec["app_field"], spec["value_map"]) if v}
    for ice, exp, gcm, scen in d[["IS", "Exp", "climate_model", "scenario"]].itertuples(index=False):
        e = M.EXP_META[(stem, ice, prefix + exp)]
        expect = (maps.get("climate_model", {}).get(gcm, gcm), maps.get("scenario", {}).get(scen, scen))
        assert (e["climate_model"], e["scenario"]) == expect, (stem, exp)


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


def test_bias_plots_order_ocean_sensitivity_pigl_first():
    """PIGL on top, then High, Medium, Low, whatever the bias; off-scale
    categories (DeConto & Pollard 2016's study-specific OCFAC) below."""
    import figures as F
    valid = (RUNS["ice_sheet"] == "AIS").to_numpy() & RUNS["publication"].isin(
        [PUBLICATION_LABEL["AIS"], "DeConto & Pollard 2016"]).to_numpy()
    b = A.bias(valid, 2015, 2023)
    w = A.run_weights(valid)
    expect = ["PIGL", "High", "Medium", "Low", "Study-specific (see experiment)"]
    fig = F.bias_figure(b, w, valid, "ocean_sensitivity", "Gt/yr")
    assert list(fig.layout.yaxis.ticktext) == expect  # drawn top-down
    rows, _ = A.low_bias_enrichment(b[valid], w[valid], RUNS["ocean_sensitivity"].to_numpy()[valid])
    bars = F.enrichment_figure(rows, "ocean_sensitivity", 0.10)
    assert list(bars.data[0].y) == [c for c in expect if c in {r["category"] for r in rows}][::-1]
    gis = (RUNS["publication"] == "Goelzer 2025 (PROTECT GIS)").to_numpy()
    cats = set(RUNS["ocean_sensitivity"].to_numpy()[gis])
    assert F.ordinal_order("ocean_sensitivity", cats) == ["Very high", "High", "Medium", "Low", "Very low"]


# ── ISMIP6-GIS (Goelzer et al. 2020), reviewed 2026-11-23 ─────────────────

def test_all_165_gis_runs_are_present():
    assert len(GIS) == 165
    assert GIS[["group", "model"]].drop_duplicates().shape[0] == 21


def test_goelzer_experiments_are_verbatim_table_1():
    e = pd.read_csv(os.path.join(M.META_DIR, "Goelzer2020_experiments.csv"), comment="#", dtype=str)
    assert list(e.columns) == ["Exp ID", "GCM", "RCP", "Ocean sensitivity", "Table 3 group", "Regional climate model"]
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


# ── PROTECT Greenland (Goelzer et al. 2025), reviewed 2026-11-23 ─────────

G25 = RUNS[RUNS["publication"] == "Goelzer 2025 (PROTECT GIS)"]


def test_goelzer2025_configurations_are_independent():
    assert len(G25) == 1270
    cfg = G25[["group", "model"]].drop_duplicates()
    assert len(cfg) == 39
    assert set(G25["ice_model"]) == {"Elmer/Ice", "IMAU-ICE", "CISM", "GISM"}  # code bases, not configurations


def test_goelzer2025_model_classes():
    sl = _by_config("sliding_law", G25)
    assert sl[("Goelzer2025-IGE", "ElmerIce2")] == "Linear"
    assert sl[("Goelzer2025-IGE", "ElmerIce3")] == "Power law (Weertman)"
    assert all(v == "Coulomb-limited" for (g, m), v in sl.items() if g in ("Goelzer2025-IMAU", "Goelzer2025-NORCE"))
    assert sl[("Goelzer2025-VUB", "GISMHOMv1")] == sl[("Goelzer2025-VUB", "GISMSIAv1")] == "Budd"
    sb = _by_config("stress_balance", G25)
    assert sb[("Goelzer2025-VUB", "GISMHOMv1")] == "Higher-order" and sb[("Goelzer2025-VUB", "GISMSIAv1")] == "Shallow"
    res = _by_config("resolution", G25)
    assert res[("Goelzer2025-IMAU", "IMAUICE1")] == ">20 km" and res[("Goelzer2025-IMAU", "IMAUICE5")] == ">5–10 km"
    assert res[("Goelzer2025-NORCE", "CISM16c-MAR312-p50")] == ">10–20 km" and res[("Goelzer2025-NORCE", "CISM02-MAR39-p50")] == "≤5 km"


def test_goelzer2025_forcing_fields():
    pct = G25["exp"].str.split("_").str[3]
    # A lower kappa percentile means more retreat (Goelzer et al. 2025: "25th (high) and 75th (low)")
    expect = {"p05": "Very high", "p25": "High", "p50": "Medium", "p75": "Low", "p95": "Very low"}
    assert (G25["ocean_sensitivity"] == pct.map(expect)).all()
    assert (G25["regional_climate_model"] == G25["exp"].str.split("_").str[2]).all()
    assert (G25["ocean_forcing"] == "Standard").all()
    assert {"CESM2-Leo", "CESM2", "UKESM1-0-LL-Robin", "UKESM1-0-LL", "ACCESS1-3"} <= set(G25["climate_model"])


def test_regional_climate_model_for_ismip6():
    assert set(AIS["regional_climate_model"]) == {"None (AOGCM output used directly)"}
    assert set(GIS["regional_climate_model"]) == {"MARv3.9"}


# ── DeConto & Pollard 2016 (review 04) ───────────────────────────────────

DP16 = RUNS[RUNS["group"] == "DeConto2016"]


def test_deconto2016_runs_are_dp16_as_published():
    """The "(BC)" rows of Kopp et al. (2017) Data Set S1: DP16's +3 degC ocean
    correction in the Amundsen and Bellingshausen seas (decision DP1). Their
    RCP8.5 2050 median, +52 mm, is Kopp et al.'s "about +5 cm under DP16"; the
    uncorrected rows give -28 mm."""
    import data as D
    df = D._load_deconto2016()
    at = df[(df["Year"] == 2050) & df["Exp"].str.startswith("RCP85_")]
    mm = at[D.MASS_COL] * D.GT_TO_MM_SLE
    assert len(at) == 29 and round(float(np.median(mm))) == 52
    assert len(DP16) == 87 and set(DP16["scenario"]) == {"RCP2.6", "RCP4.5", "RCP8.5"}


def test_deconto2016_model_and_forcing():
    assert (DP16["review_status"] == "reviewed").all()
    m = M.model_meta("AIS", "DeConto2016", "PSU3D-ICE")
    assert (m["stress_balance"], m["sliding_law"], m["initialization"], m["resolution"]) == (
        "Shallow", "Power law (Weertman)", "Hybrid (thickness-nudged)", ">5–10 km")
    assert set(DP16["climate_model"]) == {"CCSM4 (ocean); GENESIS v3 (atmosphere)"}
    assert set(DP16["regional_climate_model"]) == {"RegCM3"}
    assert set(DP16["ocean_forcing"]) == {"Open"}
    assert set(DP16["ocean_sensitivity"]) == {"Study-specific (see experiment)"}
    assert not DP16["ice_shelf_fracture"].any()
    e = M.EXP_META[("DeConto2016", "AIS", "RCP26_m23")]
    assert e["extras"]["OCFAC"] == "10" and e["extras"]["VCLIF (km/yr)"] == "0"


# ── DeConto et al. 2021 (review 05) ──────────────────────────────────────

DC21 = RUNS[RUNS["group"] == "DeConto2021"]


def test_deconto2021_is_the_calibrated_ensemble():
    """109 members x 4 scenarios; parameters are CALVLIQ and VCLIFF, whose
    means are Extended Data Table 1's calibrated averages (107, 7.7)."""
    assert len(DC21) == 436 and set(DC21["scenario"]) == {"+1.5 °C", "+2.0 °C", "+3.0 °C", "RCP8.5"}
    ex = [M.EXP_META[("DeConto2021", "AIS", e)]["extras"] for e in DC21["exp"]]
    calv = np.array([float(x["CALVLIQ"]) for x in ex])
    vcl = np.array([float(x["VCLIFF (km/yr)"]) for x in ex])
    assert round(calv.mean()) == 107 and round(vcl.mean(), 1) == 7.7


def test_deconto2021_model_and_forcing():
    assert (DC21["review_status"] == "reviewed").all()
    m = M.model_meta("AIS", "DeConto2021", "PSU3D-ICE")
    assert (m["stress_balance"], m["sliding_law"], m["initialization"], m["resolution"]) == (
        "Shallow", "Power law (Weertman)", "Hybrid (thickness-nudged)", ">5–10 km")
    assert set(DC21["climate_model"]) == {"CCSM4 (ocean); RegCM3 CO2 snapshots (atmosphere)"}
    assert set(DC21["regional_climate_model"]) == {"RegCM3"}
    assert set(DC21["ocean_forcing"]) == {"Open"} and set(DC21["ocean_sensitivity"]) == {"Medium"}
    assert not DC21["ice_shelf_fracture"].any() and BIAS_ELIGIBLE[DC21.index].all()  # decision DC3: kept
    # The protocol (incl. the partial IMBIE calibration) is documented in the file, not the hover
    e = pd.read_csv(os.path.join(M.META_DIR, "DeConto2021_experiments.csv"), comment="#", dtype=str)
    assert e["Protocol"].str.contains("matching IMBIE").all()
    assert not any("Protocol:" in h for h in pd.concat([DC21["hover"], RUNS.loc[RUNS["group"] == "DeConto2016", "hover"]]))


def test_hover_lines_stay_narrow():
    """Long hover lines are wrapped so the hover box doesn't cover the plot."""
    longest = max(len(line) for h in RUNS["hover"] for line in h.split("<br>") if "<" not in line)
    assert longest <= 60


# ── Coulon et al. 2024 (review 06) ───────────────────────────────────────

C24 = RUNS[RUNS["group"] == "Coulon2024"]


def test_coulon2024_model_and_forcing():
    assert len(C24) == 100 and (C24["review_status"] == "reviewed").all()
    m = M.model_meta("AIS", "Coulon2024", "Kori-ULB")
    assert (m["stress_balance"], m["sliding_law"], m["initialization"], m["resolution"]) == (
        "Shallow", "Power law (Weertman)", "Hybrid (thickness-nudged)", ">10–20 km")
    # one code base with Seroussi 2020's ULB rows (decision C1)
    assert m["ice_model"] == M.model_meta("AIS", "ULB", "fETISh_16km")["ice_model"]
    assert set(C24["climate_model"]) == {"MRI-ESM2-0", "UKESM1-0-LL", "CESM2-WACCM", "IPSL-CM6A-LR"}
    assert set(C24["ocean_forcing"]) == {"Open"}
    assert set(C24["ocean_sensitivity"]) == {"Study-specific (see experiment)"}
    assert set(C24["regional_climate_model"]) == {"None (AOGCM output used directly)"}
    assert not C24["ice_shelf_fracture"].any()


def test_coulon2024_keeps_the_2014_mass_change():
    """Each member's 2014->2015 change is its 2014 net mass balance from the
    archive, not zero (decision C4)."""
    d = pd.read_csv(os.path.join(os.path.dirname(M.META_DIR), "data", "external_sources_coulon2024_ais.csv.gz"))
    wide = d.pivot(index="Year", columns="Exp", values="Cumulative ice sheet mass change (Gt)")
    step = wide.loc[2015] - wide.loc[2014]
    assert (step != 0).all() and -200 < step.median() < -60


def test_gis_high_retreat_sensitivity_loses_more_ice():
    """The encoded sensitivity follows the physics: for the same configuration
    and forcing, High (and Very high) runs lose more ice by 2100 than Low
    (Very low) ones."""
    d = pd.read_csv(os.path.join(os.path.dirname(M.META_DIR), "data", "external_sources_goelzer2025_gis.csv.gz"))
    y = d[d["Year"].isin([2015, 2100])].pivot_table(index=["Group", "Model", "Exp"], columns="Year",
                                                     values="Cumulative ice sheet mass change (Gt)").reset_index()
    y["loss"] = y[2015] - y[2100]
    y["forcing"] = y["Exp"].str.rsplit("_", n=1).str[0]
    y["level"] = [M.EXP_META[("Goelzer2025", "GIS", e)]["ocean_sensitivity"] for e in y["Exp"]]
    w = y.pivot_table(index=["Group", "Model", "forcing"], columns="level", values="loss")
    hl = w.dropna(subset=["High", "Low"])
    assert len(hl) > 100 and (hl["High"] > hl["Low"]).mean() > 0.95
    vv = w.dropna(subset=["Very high", "Very low"])
    assert (vv["Very high"] > vv["Very low"]).mean() > 0.9


# ── Rahlves et al. 2025 (review 07) ──────────────────────────────────────

R25 = RUNS[RUNS["group"] == "Rahlves2025"]


def test_rahlves2025_configurations_and_runs():
    """Two configurations (decision R1); 16 Table A1 projections x 3 retreat
    sensitivities each; no control runs (decision R2)."""
    assert set(R25["model"]) == {"CISM (ERA5-init)", "CISM (ESM-init)"}
    assert (R25.groupby("model").size() == 48).all()
    assert "Control" not in set(R25["scenario"])
    assert set(R25["ocean_sensitivity"]) == {"Low", "Medium", "High"}
    # ACCESS1.3 is CMIP5: RCP8.5 in both branches
    assert set(R25.loc[R25["climate_model"] == "ACCESS1-3", "scenario"]) == {"RCP8.5"}
    assert set(R25["climate_model"]) == {"ACCESS1-3", "CNRM-CM6-1", "UKESM1-0-LL-Robin", "CESM2-Leo", "CNRM-ESM2-1",
                                         "MPI-ESM1-2-HR", "IPSL-CM6A-LR", "NorESM2-MM", "CESM2", "UKESM1-0-LL"}


def test_rahlves2025_model_classes():
    for model in ("CISM (ERA5-init)", "CISM (ESM-init)"):
        m = M.model_meta("GIS", "Rahlves2025", model)
        assert (m["stress_balance"], m["sliding_law"], m["initialization"], m["resolution"]) == (
            "Higher-order", "Budd", "Hybrid (thickness-nudged)", "≤5 km")
    assert set(R25["ocean_forcing"]) == {"Standard"} and set(R25["regional_climate_model"]) == {"MARv3.12"}


# ── Aschwanden et al. 2019 (review 08) ───────────────────────────────────

ASCH = RUNS[RUNS["group"] == "Aschwanden2019"]


def test_aschwanden2019_runs_and_model():
    """A19's large ensemble: 489 + 491 + 478 real PISM runs under RCP2.6/4.5/8.5."""
    assert len(ASCH) == 1458 and (ASCH["review_status"] == "reviewed").all()
    assert ASCH["scenario"].value_counts().to_dict() == {"RCP4.5": 491, "RCP2.6": 489, "RCP8.5": 478}
    m = M.model_meta("GIS", "Aschwanden2019", "PISM")
    assert (m["stress_balance"], m["sliding_law"], m["initialization"], m["resolution"]) == (
        "Shallow", "Coulomb-limited", "Hybrid (thickness-nudged)", "≤5 km")
    # same initial state and sliding law as UAF-PISM in ISMIP6-GIS
    uaf = M.model_meta("GIS", "UAF", "PISM1")
    assert (m["initialization"], m["sliding_law"]) == (uaf["initialization"], uaf["sliding_law"])


def test_aschwanden2019_forcing_fields():
    assert set(ASCH["climate_model"]) == {"GISS-E2-H", "GISS-E2-R", "IPSL-CM5A-LR", "MPI-ESM-LR"}
    assert set(ASCH["regional_climate_model"]) == {"None (AOGCM output used directly)"}
    assert set(ASCH["ocean_forcing"]) == {"Open"}
    assert set(ASCH["ocean_sensitivity"]) == {"Study-specific (see experiment)"}
    assert not ASCH["ice_shelf_fracture"].any() and BIAS_ELIGIBLE[ASCH.index].all()
    e = M.EXP_META[("Aschwanden2019", "GIS", "les000_rcp26")]["extras"]
    assert set(e) >= {"OCM", "OCS", "q"} and e["OCM"] in {"low", "mid", "high"}
    assert "PRS (%/K)" not in e  # file-only parameters stay out of the hover


# ── Seroussi et al. 2024, ISMIP6 Antarctica 2300 (review 09) ─────────────

S24 = RUNS[RUNS["publication"] == "Seroussi 2024 (ISMIP6 AIS 2300)"]


def test_seroussi2024_runs_and_experiments():
    """Table 4's counts without the repeat-forcing duplicates expAE06/08/09
    (decision S3); collapse runs expAE11-14 out of the bias analysis."""
    assert len(S24) == 329 and S24[["group", "model"]].drop_duplicates().shape[0] == 43
    assert not set(S24["exp"]) & {"expAE06", "expAE08", "expAE09", "ctrlAE"}
    assert S24.groupby("exp").size().to_dict() == {
        "expAE01": 43, "expAE02": 43, "expAE03": 43, "expAE04": 43, "expAE05": 43, "expAE07": 29, "expAE10": 29,
        "expAE11": 14, "expAE12": 14, "expAE13": 14, "expAE14": 14}
    collapse = S24["exp"].isin(["expAE11", "expAE12", "expAE13", "expAE14"])
    assert (S24["ice_shelf_fracture"] == collapse).all() and not BIAS_ELIGIBLE[S24.index[collapse]].any()
    assert set(S24["climate_model"]) == {"NorESM1-M", "CCSM4", "HadGEM2-ES", "CESM2-WACCM", "UKESM1-0-LL"}


def test_seroussi2024_encodings():
    m = lambda g, mod: M.model_meta("AIS", "Seroussi2024-" + g, mod)  # noqa: E731
    assert m("NORCE", "CISM2-MAR364-ERA-t1")["resolution"] == "≤5 km"  # 4 km, decision S4
    assert m("VUW", "PISM1")["initialization"] == "Hybrid (thickness-nudged)"  # thickness resetting, rule D1
    assert m("PIK", "PISM")["initialization"] == "Spin-up / equilibrium"
    assert m("DOE", "MALI_8km_Ant95")["sliding_law"] == "Power law (Weertman)"
    assert m("UCM", "Yelmo")["sliding_law"] == "Coulomb-limited" and m("VUB", "AISMPALEO")["sliding_law"] == "Budd"
    _, e = M.run_meta("AIS", "Seroussi2024-DOE", "MALI_8km_Ant95", "expAE02")
    assert (e["ocean_forcing"], e["ocean_sensitivity"]) == ("Standard", "High")
    _, e = M.run_meta("AIS", "Seroussi2024-PIK", "PISM", "expAE02")
    assert (e["ocean_forcing"], e["ocean_sensitivity"]) == ("Open", "Medium")  # PICO


def test_seroussi2024_is_weighted_as_ismip6():
    """Each configuration counts on its own (decision S5): every run weighs 1,
    as for Seroussi 2020, while a non-ISMIP6 study weighs as one configuration."""
    valid = (RUNS["ice_sheet"] == "AIS").to_numpy() & RUNS["publication"].isin(
        [PUBLICATION_LABEL["AIS"], "Seroussi 2024 (ISMIP6 AIS 2300)", "Coulon 2024"]).to_numpy()
    w = A.run_weights(valid)
    assert np.allclose(w[S24.index], 1.0)
    typical = (188 + 329) / (16 + 43)
    assert np.isclose(w[(RUNS["publication"] == "Coulon 2024").to_numpy()].sum(), typical)
