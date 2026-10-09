"""
Data layer: loads every bundled dataset once at import time and flattens all
simulations (ISMIP6 + every extra source) into one run table plus aligned
per-run arrays, so the analysis/figure code never has to re-group pandas
frames per request.

Everything is read from dash_app/data/ -- including IMBIE3 and IMBIE2, bundled there
rather than fetched from ramadda.data.bas.ac.uk at startup, since that host
being down (it returned site-wide 503s on 2026-10-02) would otherwise stop the
app from starting at all.
"""

import os

import numpy as np
import pandas as pd
import plotly.express as px

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

proj_start = 2015  # utilities/helper.py
GT_TO_MM_SLE = -1.0 / 362.5  # Gt -> mm sea-level equivalent; negative because ice loss raises sea level
ICE_SHEET_AREA_M2 = {"AIS": 12.1e12, "GIS": 1.15e12}  # for area-normalizing AIS+GIS bias (see analysis.py)
ICE_SHEETS = ["AIS", "GIS"]
ICE_SHEET_NAMES = {"AIS": "Antarctic Ice Sheet", "GIS": "Greenland Ice Sheet"}
MASS_COL = "Cumulative ice sheet mass change (Gt)"


def _load_imbie2026(region):
    """IMBIE3 (Otosaka et al., 2026, Scientific Data,
    https://doi.org/10.1038/s41597-026-08088-0) headline mass balance for
    "antarctica"/"greenland", re-zeroed at proj_start. Same transform as
    utilities/imbie2026_loader.py, totals only."""
    df = pd.read_csv(os.path.join(DATA_DIR, f"imbie3_{region}_Gt_partitioned.csv.gz"), comment="#")
    date = pd.to_datetime(df["Date"])
    df["Year"] = date.dt.year + (date.dt.month - 1) / 12
    imbie = df.rename(columns={
        "Cumulative mass balance anomaly (Gt)": MASS_COL,
        "Cumulative mass balance anomaly uncertainty (Gt)": "Cumulative ice sheet mass change uncertainty (Gt)",
    })[["Year", MASS_COL, "Cumulative ice sheet mass change uncertainty (Gt)"]].copy()
    imbie[MASS_COL] -= imbie.loc[imbie["Year"] == proj_start, MASS_COL].values
    return imbie


def _load_imbie2023(region):
    """IMBIE2 (Otosaka et al., 2023, ESSD, https://doi.org/10.5194/essd-15-1597-2023;
    data DOI 10.5285/77B64C55-7166-4A06-9DEF-2E400398E452, Open Government
    Licence v3.0), 1992-2020, re-zeroed at proj_start. Raw files as served
    by the BAS Polar Data Centre (URLs in utilities/imbie2023_loader.py).
    Its cumulative uncertainty accumulates in quadrature like IMBIE3's
    (checked: sqrt(cumsum(rate_unc**2 / 12)) reproduces it to 1e-4 Gt), so
    analysis.imbie_timeseries treats both products alike."""
    df = pd.read_csv(os.path.join(DATA_DIR, f"imbie2_{region}_2021_Gt.csv.gz"))
    imbie = df.rename(columns={
        "Cumulative mass balance (Gt)": MASS_COL,
        "Cumulative mass balance uncertainty (Gt)": "Cumulative ice sheet mass change uncertainty (Gt)",
    })[["Year", MASS_COL, "Cumulative ice sheet mass change uncertainty (Gt)"]].copy()
    imbie[MASS_COL] -= imbie.loc[imbie["Year"] == proj_start, MASS_COL].values
    return imbie


# Observational benchmarks, selectable in the app. last_year = the last full
# calendar year of data, which caps the rate window.
OBS_PRODUCTS = {
    "imbie3": {"label": "IMBIE3", "citation": "Otosaka et al. (2026)",
               "url": "https://doi.org/10.1038/s41597-026-08088-0", "last_year": 2023},
    "imbie2": {"label": "IMBIE2", "citation": "Otosaka et al. (2023)",
               "url": "https://doi.org/10.5194/essd-15-1597-2023", "last_year": 2020},
}
DEFAULT_OBS = "imbie3"
IMBIE = {
    "imbie3": {"AIS": _load_imbie2026("antarctica"), "GIS": _load_imbie2026("greenland")},
    "imbie2": {"AIS": _load_imbie2023("antarctica"), "GIS": _load_imbie2023("greenland")},
}


# ─────────────────────────────────────────────────────────────────────────────
# Experiment metadata look-up tables
# Sources:
#   GIS – Goelzer et al. (2020) The Cryosphere 14, 3071-3096
#   AIS – Seroussi et al. (2020) The Cryosphere 14, 3033-3070
# ─────────────────────────────────────────────────────────────────────────────

# GIS experiment descriptions (Goelzer et al. 2020, Table 1 & Appendix)
gis_exp_meta = {
    # "exp01":  {"climate_model": "MIROC5",       "scenario": "RCP8.5", "protocol": "Open",     "ocean_sensitivity": "Medium"},
    # "exp02":  {"climate_model": "MIROC5",       "scenario": "RCP8.5", "protocol": "Open",     "ocean_sensitivity": "Low"},
    # "exp03":  {"climate_model": "MIROC5",       "scenario": "RCP2.6", "protocol": "Open",     "ocean_sensitivity": "Medium"},
    # "exp04":  {"climate_model": "MIROC5",       "scenario": "RCP2.6", "protocol": "Open",     "ocean_sensitivity": "Low"},
    "exp05":  {"climate_model": "MIROC5",       "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp06":  {"climate_model": "NorESM",       "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp07":  {"climate_model": "MIROC5",       "scenario": "RCP2.6", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp08":  {"climate_model": "HadGEM2-ES",   "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp09":  {"climate_model": "MIROC5",       "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "High"},
    "exp10":  {"climate_model": "MIROC5",       "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Low"},
    # "exp11":  {"climate_model": "ACCESS1.3",    "scenario": "RCP8.5", "protocol": "Open",     "ocean_sensitivity": "Medium"},
    # "exp12":  {"climate_model": "ACCESS1.3",    "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    # "exp13":  {"climate_model": "CESM2",        "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "High"},
    "expa01": {"climate_model": "IPSL-CM5A-MR", "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "expa02": {"climate_model": "CSIRO-Mk3.6",  "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "expa03": {"climate_model": "ACCESS1.3",    "scenario": "RCP8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
}

# AIS experiment descriptions (Seroussi et al. 2020, Table 1)
ais_exp_meta = {
    "exp01": {"climate_model": "NorESM",         "scenario": "RCP8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp02": {"climate_model": "MIROC-ESM-CHEM", "scenario": "RCP8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp03": {"climate_model": "NorESM",         "scenario": "RCP2.6",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp04": {"climate_model": "CCSM4",          "scenario": "RCP8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp05": {"climate_model": "NorESM",         "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp06": {"climate_model": "MIROC-ESM-CHEM", "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp07": {"climate_model": "NorESM",         "scenario": "RCP2.6",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp08": {"climate_model": "CCSM4",          "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp09": {"climate_model": "NorESM",         "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "PIGL medium"},
    "exp10": {"climate_model": "NorESM",         "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "PIGL high"},
    "exp11": {"climate_model": "CCSM4",          "scenario": "RCP8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp12": {"climate_model": "CCSM4",          "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp13": {"climate_model": "NorESM",         "scenario": "RCP8.5",  "protocol": "Standard", "basal_melt_param": "PIGL very high"},
    "expA1": {"climate_model": "HadGEM2-ES",     "scenario": "SSP5-8.5", "protocol": "Open",     "basal_melt_param": "Standard"},
    "expA2": {"climate_model": "CSIRO-MK3",      "scenario": "SSP5-8.5", "protocol": "Open",     "basal_melt_param": "Standard"},
    "expA3": {"climate_model": "IPSL-CM5A-MR",   "scenario": "SSP1-2.6", "protocol": "Open",     "basal_melt_param": "Standard"},
    "expA4": {"climate_model": "IPSL-CM5A-MR",   "scenario": "SSP1-2.6", "protocol": "Open",     "basal_melt_param": "Standard"},
    "expA5": {"climate_model": "HadGEM2-ES",     "scenario": "SSP5-8.5", "protocol": "Standard", "basal_melt_param": "Standard"},
    "expA6": {"climate_model": "CSIRO-MK3",      "scenario": "SSP5-8.5", "protocol": "Standard", "basal_melt_param": "Standard"},
    "expA7": {"climate_model": "IPSL-CM5A-MR",   "scenario": "SSP1-2.6", "protocol": "Standard", "basal_melt_param": "Standard"},
    "expA8": {"climate_model": "IPSL-CM5A-MR",   "scenario": "SSP1-2.6", "protocol": "Standard", "basal_melt_param": "Standard"},
}

# ─────────────────────────────────────────────────────────────────────────────
# Ice sheet model (Group + Model) metadata
# Sources: Goelzer et al. 2020 Table A1; Seroussi et al. 2020 Table A1
# ─────────────────────────────────────────────────────────────────────────────
ism_meta = {
    # Group          Model        ice_sheet_model  sliding_law                  initialization
    ("AWI",      "ISSM1"):     {"ice_model": "ISSM",      "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("AWI",      "ISSM2"):     {"ice_model": "ISSM",      "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("DMI",      "PISM"):      {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic (Budd)",    "initialization": "Spin-up"},
    ("ILTS_PIK", "SICOPOLIS"): {"ice_model": "SICOPOLIS", "sliding_law": "Weertman",                 "initialization": "Spin-up"},
    ("IMAU",     "IMAUICE1"):  {"ice_model": "IMAU-ICE",  "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("IMAU",     "IMAUICE2"):  {"ice_model": "IMAU-ICE",  "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("JPL1",     "ISSM"):      {"ice_model": "ISSM",      "sliding_law": "Budd / Schoof",            "initialization": "Data assimilation"},
    ("LSCE",     "GRISLI"):    {"ice_model": "GRISLI",    "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("MIROC",    "ICIES1"):    {"ice_model": "IcIES",     "sliding_law": "Weertman",                 "initialization": "Spin-up"},
    ("MIROC",    "ICIES2"):    {"ice_model": "IcIES",     "sliding_law": "Weertman",                 "initialization": "Spin-up"},
    ("NCAR",     "CISM"):      {"ice_model": "CISM",      "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
    ("PIK",      "PISM1"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("PIK",      "PISM2"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("UAF",      "PISM1"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic (Coulomb)", "initialization": "Spin-up"},
    ("UAF",      "PISM2"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic (Coulomb)", "initialization": "Spin-up"},
    ("UCIJPL",   "ISSM"):      {"ice_model": "ISSM",      "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
    ("ULB",      "FETISH1"):   {"ice_model": "Kori-ULB/f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
    ("ULB",      "FETISH2"):   {"ice_model": "Kori-ULB/f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
    ("UNN",      "ElmerIce"):  {"ice_model": "Elmer/Ice", "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
    ("VUW",      "PISM"):      {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic (Coulomb)", "initialization": "Spin-up"},
    # Tim edits to Greenland
    ("BGC",      "BISICLES"):  {"ice_model": "BISICLES",  "sliding_law": "Linear viscous",              "initialization": "Data assimilation"},
    ("MUN",      "GSM1"):      {"ice_model": "GSM",       "sliding_law": "Coulomb and Weertman",        "initialization": "Spin-up"},
    ("MUN",      "GSM2"):      {"ice_model": "GSM",       "sliding_law": "Linear viscous and Weertman", "initialization": "Spin-up"},
    ("VUB",      "GISM"):      {"ice_model": "GISM",      "sliding_law": "Weertman",                    "initialization": "Data assimilation"},
    # Goelzer et al. 2020 Table A1/A6/A7/A12 -- previously missing exact
    # keys, silently mislabeled via get_ism_meta()'s substring fallback
    # (e.g. JPL/ISSM was borrowing JPL1's unrelated AIS entry).
    ("JPL",      "ISSM"):      {"ice_model": "ISSM",      "sliding_law": "Linear viscous",              "initialization": "Data assimilation"},
    ("JPL",      "ISSMPALEO"): {"ice_model": "ISSM",      "sliding_law": "Linear viscous",              "initialization": "Spin-up"},
    ("UCIJPL",   "ISSM1"):     {"ice_model": "ISSM",      "sliding_law": "Linear viscous",              "initialization": "Data assimilation"},
    ("UCIJPL",   "ISSM2"):     {"ice_model": "ISSM",      "sliding_law": "Linear viscous",              "initialization": "Data assimilation"},
    # AIS additional groups
    ("ARC",      "PISM1"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("ARC",      "PISM2"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("AWI",      "PISM1"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("DOE",      "MALI"):      {"ice_model": "MALI",      "sliding_law": "Coulomb",                  "initialization": "Data assimilation"},
    ("GRL",      "PISM"):      {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("GSFC",     "ISSM"):      {"ice_model": "ISSM",      "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("IGE",      "ElmerIce"):  {"ice_model": "Elmer/Ice", "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
    ("ILTS",     "SICOPOLIS"): {"ice_model": "SICOPOLIS", "sliding_law": "Weertman",                 "initialization": "Spin-up"},
    ("NEMO",     "fETISh"):    {"ice_model": "Kori-ULB/f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
    ("PSU",      "PSU3D1"):    {"ice_model": "PSU-ISM",   "sliding_law": "Coulomb / Weertman",       "initialization": "Spin-up"},
    ("PSU",      "PSU3D2"):    {"ice_model": "PSU-ISM",   "sliding_law": "Coulomb / Weertman",       "initialization": "Spin-up"},
    ("UTAS",     "ElmerIce"):  {"ice_model": "Elmer/Ice", "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
    ("VUB",      "AISMPALEO"): {"ice_model": "AISMPALEO", "sliding_law": "Weertman",                 "initialization": "Spin-up"},

    # Non-ISMIP6 published simulations (see utilities/external_sources.py for
    # provenance/derivation of each). "initialization" for Rahlves2025 covers
    # both its ERA5- and ESM-forced branches with one description, since Model
    # is "CISM" either way (unlike ISMIP6, this dict doesn't split them further).
    ("Rahlves2025", "CISM"):     {"ice_model": "CISM",      "sliding_law": "Weertman",                    "initialization": "Spin-up"},
    ("Coulon2024",  "Kori-ULB"): {"ice_model": "Kori-ULB/f.ETISh",  "sliding_law": "Weertman",                    "initialization": "Data assimilation"},
    ("Aschwanden2019", "PISM"):  {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",              "initialization": "Spin-up"},

    # Goelzer et al. (2025) PROTECT-Greenland ensemble (see
    # utilities/external_sources.py for provenance/scope decisions) -- one
    # entry per lab's ice sheet model (Group="Goelzer2025" throughout;
    # NORCE's many internal CISM resolution/tuning variants are collapsed
    # to Model="CISM", per the user's explicit choice), sliding
    # law/initialization transcribed from the paper's own Table 1/Sect. 2.
    ("Goelzer2025", "Elmer/Ice"): {"ice_model": "Elmer/Ice", "sliding_law": "Linear & Weertman (m=1/3)",   "initialization": "Inverse control method + 20-yr relaxation"},
    ("Goelzer2025", "IMAU-ICE"):  {"ice_model": "IMAU-ICE",  "sliding_law": "Basal inversion (variable)",  "initialization": "Hybrid: basal inversion + paleo spin-up"},
    ("Goelzer2025", "CISM"):      {"ice_model": "CISM",      "sliding_law": "Schoof (2005)",               "initialization": "Spin-up"},
    ("Goelzer2025", "GISM"):      {"ice_model": "GISM",      "sliding_law": "Optimized coefficients (variable)", "initialization": "Iterative assimilation + 2-cycle spin-up"},

    # Edwards et al. (2021) -- a Gaussian-process EMULATOR, not a physical
    # ice-sheet model, so it has no sliding law/spin-up to report; labeled
    # honestly rather than borrowed from an unrelated model (see
    # utilities/external_sources.py for provenance/scope decisions).
    ("Edwards2021", "emulandice"): {"ice_model": "emulandice",
                                     "sliding_law": "Not applicable (statistical emulator)",
                                     "initialization": "Not applicable (statistical emulator)"},

    # DeConto & Pollard (2016). Data Set S1 gives only the model and each
    # run's MICI parameters (shown in hover), not the sliding law or
    # initialization, so those defer to the paper rather than being guessed.
    ("DeConto2016", "PSU3D-ICE"): {"ice_model": "PSU3D-ICE", "sliding_law": "See paper",
                                   "initialization": "See paper"},
    # DeConto et al. (2021): same model (updated), same caveat.
    ("DeConto2021", "PSU3D-ICE"): {"ice_model": "PSU3D-ICE", "sliding_law": "See paper",
                                   "initialization": "See paper"},
}


def get_exp_meta(ice_sheet, exp):
    """Return experiment metadata dict for the given ice sheet and experiment code."""
    meta = ais_exp_meta if ice_sheet == "AIS" else gis_exp_meta
    return meta.get(exp, {"climate_model": "Unknown", "scenario": "Unknown", "protocol": "Unknown"})


def get_ism_meta(group, model):
    """Return ice sheet model metadata; falls back gracefully if unknown."""
    key = (group, model)
    if key in ism_meta:
        return ism_meta[key]
    # Try case-insensitive partial match on group
    for (g, m), v in ism_meta.items():
        if g.upper() in group.upper() or group.upper() in g.upper():
            return v
    return {"ice_model": model, "sliding_law": "See paper", "initialization": "See paper"}


# ─────────────────────────────────────────────────────────────────────────────
# Non-ISMIP6 published simulations -- toggleable overlays on the interactive
# rate-comparison figure. Self-contained like the rest of this file (reads
# the same bundled CSVs as utilities/external_sources.py's loaders, which
# this deliberately doesn't import -- see that module's docstrings for what
# each source is and how the bundled CSV was derived from the original
# archive; this app just reads the result from its own data/ directory so it
# stays deployable by copying dash_app/ alone).
# ─────────────────────────────────────────────────────────────────────────────

_TEXT_COLS = ["IS", "Group", "Model", "Exp", "climate_model", "scenario", "protocol", "RCP",
              "ocean_sensitivity", "basal_melt_param", "retreat_percentile"]


def _read_csv(name):
    """Bundled CSV with its repeated text columns parsed as categoricals --
    as plain object strings, the 255k-row Edwards files alone cost ~100 MB
    of short-lived Python strings at startup (Render caps the process at
    512 MB)."""
    path = os.path.join(DATA_DIR, name)
    cols = pd.read_csv(path, nrows=0).columns
    return pd.read_csv(path, dtype={c: "category" for c in _TEXT_COLS if c in cols})


def _exp_meta_from_df(df, extra_cols):
    cols = ["Exp", "climate_model", "scenario", "protocol"] + extra_cols
    unique = df[cols].drop_duplicates("Exp").set_index("Exp")
    return unique.to_dict(orient="index")


ismip6_ais = _read_csv("ismip6_ais.csv.gz")
ismip6_ais["IS"] = "AIS"
ismip6_gis = _read_csv("ismip6_gis_ctrl.csv.gz")
ismip6_gis["IS"] = "GIS"
rahlves2025_gis = _read_csv("external_sources_rahlves2025_gis.csv.gz")
coulon2024_ais = _read_csv("external_sources_coulon2024_ais.csv.gz")
aschwanden2022_gis = _read_csv("external_sources_aschwanden2022_gis.csv.gz")
goelzer2025_gis = _read_csv("external_sources_goelzer2025_gis.csv.gz")
edwards2021_ais = _read_csv("external_sources_edwards2021_ais.csv.gz")
edwards2021_gis = _read_csv("external_sources_edwards2021_gis.csv.gz")
# Risk-averse Antarctic variant (Edwards et al. 2021's S11_RISK sensitivity
# test). Its sample ids repeat MAIN's, so prefix them: exp metadata is keyed
# by (ice sheet, Exp), and the two variants must not overwrite each other.
edwards2021_ais_risk = _read_csv("external_sources_edwards2021_ais_risk.csv.gz")
edwards2021_ais_risk["Exp"] = "RISK_" + edwards2021_ais_risk["Exp"].astype(str)


def _load_deconto2016():
    """DeConto & Pollard (2016) Antarctic runs from Kopp et al. (2017,
    Earth's Future) Data Set S1, bundled byte-for-byte as downloaded (its
    CC BY-NC-ND licence allows sharing it unmodified, so the reshaping
    happens here at startup). Same transform as
    utilities/external_sources.py load_deconto2016_ais(), which documents
    provenance: raw (uncorrected) RCP2.6/4.5/8.5 runs, AIS total, decadal mm
    GMSL (zeroed at 2000) linearly interpolated to annual through 2100 and
    converted to Gt (x -362.5)."""
    raw = pd.read_csv(os.path.join(DATA_DIR, "eft2271-sup-0002-2017ef000663-ds01.tsv"), sep="\t")
    scen = {"RCP 2.6": "RCP2.6", "RCP 4.5": "RCP4.5", "RCP 8.5": "RCP8.5"}
    raw = raw[(raw["Ice Sheet"] == "AIS") & raw["Scenario"].isin(scen)]
    dec = np.array([int(c) for c in raw.columns if c.isdigit() and int(c) <= 2100])
    years = np.arange(dec.min(), dec.max() + 1)
    vals = raw[[str(y) for y in dec]].to_numpy(dtype=float)
    out = []
    for (_, r), v in zip(raw.iterrows(), vals):
        sc = scen[r["Scenario"]]
        out.append(pd.DataFrame({
            "Year": years, MASS_COL: np.interp(years, dec, v) / GT_TO_MM_SLE,  # mm GMSL -> Gt (loss negative)
            "Group": "DeConto2016", "Model": "PSU3D-ICE", "Exp": f"{sc.replace('.', '')}_m{int(r['Ensemble Member']):02d}",
            "climate_model": "Not reported in Data Set S1", "scenario": sc,
            "protocol": "DP16 members passing Pliocene (5-15 m) and LIG tests; uncorrected",
            "mici_params": (f"OCFAC {r['OCFAQ']:g}, CREVLIQ {r['CREVLIQ']:g}, VCLIF {r['VCLIF']:g} km/yr; "
                            f"Pliocene {r['Pliocene SL (m)']:.1f} m, LIG {r['LIG SL (m)']:.1f} m"),
        }))
    return pd.concat(out, ignore_index=True)


deconto2016_ais = _load_deconto2016()
deconto2021_ais = _read_csv("external_sources_deconto2021_ais.csv.gz")
gis_exp_meta.update(_exp_meta_from_df(rahlves2025_gis, ["ocean_sensitivity"]))
ais_exp_meta.update(_exp_meta_from_df(coulon2024_ais, ["basal_melt_param"]))
gis_exp_meta.update(_exp_meta_from_df(aschwanden2022_gis, []))
gis_exp_meta.update(_exp_meta_from_df(goelzer2025_gis, ["retreat_percentile"]))
ais_exp_meta.update(_exp_meta_from_df(edwards2021_ais, []))
ais_exp_meta.update(_exp_meta_from_df(edwards2021_ais_risk, []))
ais_exp_meta.update(_exp_meta_from_df(deconto2016_ais, ["mici_params"]))
ais_exp_meta.update(_exp_meta_from_df(deconto2021_ais, ["mici_params"]))
gis_exp_meta.update(_exp_meta_from_df(edwards2021_gis, []))

# Display order everywhere (sidebar, Study legends/boxes, README): ISMIP6
# first (SOURCE_LABELS below), then the other ice-sheet-model studies by
# publication date, then the Edwards et al. (2021) emulator studies last,
# since they aren't ice sheet model runs themselves (user decision, 2026-10-08).
EXTRA_SOURCES = [
    {"label": "DeConto & Pollard 2016", "df": deconto2016_ais, "color": "#7f3b08"},
    {"label": "Aschwanden 2019", "df": aschwanden2022_gis, "color": "#756bb1"},
    {"label": "DeConto 2021", "df": deconto2021_ais, "color": "#d6604d"},
    {"label": "Coulon 2024", "df": coulon2024_ais, "color": "#31a354"},
    {"label": "Rahlves 2025", "df": rahlves2025_gis, "color": "#e6550d"},  # TC 19, 1205 (Mar 2025)
    {"label": "Goelzer 2025 (PROTECT GIS)", "df": goelzer2025_gis, "color": "#3182bd"},  # TC 19, 6887 (later 2025)
    # Edwards et al. (2021) covers BOTH ice sheets, so its labels carry an
    # ice-sheet suffix to stay unambiguous. Same color for AIS Main and GIS --
    # they render on separate AIS/GIS panels, so there's no legend collision.
    {"label": "Edwards 2021 (AIS Main)", "df": edwards2021_ais, "color": "#e7298a"},
    {"label": "Edwards 2021 (AIS Risk Averse)", "df": edwards2021_ais_risk, "color": "#a6114f"},
    {"label": "Edwards 2021 (GIS)", "df": edwards2021_gis, "color": "#e7298a"},
]

for _src, _is in [
    ("Rahlves 2025", "GIS"), ("Coulon 2024", "AIS"), ("Aschwanden 2019", "GIS"),
    ("Goelzer 2025 (PROTECT GIS)", "GIS"), ("Edwards 2021 (AIS Main)", "AIS"),
    ("Edwards 2021 (AIS Risk Averse)", "AIS"), ("Edwards 2021 (GIS)", "GIS"), ("DeConto & Pollard 2016", "AIS"),
    ("DeConto 2021", "AIS"),
]:
    next(s for s in EXTRA_SOURCES if s["label"] == _src)["df"]["IS"] = _is

PUBLICATION_LABEL = {"AIS": "Seroussi 2020 (ISMIP6 AIS)", "GIS": "Goelzer 2020 (ISMIP6 GIS)"}
ISMIP6_COLOR = "#636363"
SOURCE_LABELS = [PUBLICATION_LABEL["AIS"], PUBLICATION_LABEL["GIS"]] + [s["label"] for s in EXTRA_SOURCES]
SOURCE_DEFAULT_CHECKED = [PUBLICATION_LABEL["AIS"], PUBLICATION_LABEL["GIS"]]
SOURCE_COLOR = {PUBLICATION_LABEL["AIS"]: ISMIP6_COLOR, PUBLICATION_LABEL["GIS"]: ISMIP6_COLOR,
                **{s["label"]: s["color"] for s in EXTRA_SOURCES}}

# (key, label) -- key is a run-table column; None = no grouping.
GROUP_DIMENSIONS = [
    (None, "All simulations"),
    ("publication", "Study"),
    ("initialization", "Initialization"),
    ("ice_model", "Ice sheet model"),
    ("sliding_law", "Sliding law"),
    ("scenario", "Climate scenario"),
    ("climate_model", "GCM"),
]
ANOVA_CHARACTERISTICS = ["ice_model", "sliding_law", "initialization", "scenario", "climate_model"]
DIM_LABEL = {k: v for k, v in GROUP_DIMENSIONS if k}

# "Collapse RCPs and SSPs": each RCP pooled with the SSP of the same nominal
# 2100 forcing (W/m^2). Not equivalent scenarios -- the app shows a warning
# whenever this is on. Scenarios without a counterpart keep their names.
SCENARIO_COMPOSITE = {
    "SSP1-1.9": "composite 1.9",
    "RCP2.6": "composite 2.6", "SSP1-2.6": "composite 2.6",
    "RCP4.5": "composite 4.5", "SSP2-4.5": "composite 4.5",
    "SSP3-7.0": "composite 7.0",
    "RCP8.5": "composite 8.5", "SSP5-8.5": "composite 8.5",
}
COMPOSITE_DIM = "scenario_composite"
DIM_LABEL[COMPOSITE_DIM] = "Climate scenario (RCP/SSP composite)"

_FIXED_COLORS = {
    "initialization": {"Data assimilation": "#1a7f5e", "Spin-up": "#003466", "See paper": "#6a6a6a"},
    "scenario": {"RCP2.6": "#003466", "RCP4.5": "#b8860b", "RCP8.5": "#990002", "SSP1-2.6": "#1a7f5e",
                 "SSP2-4.5": "#8a6d3a", "SSP5-8.5": "#8b1a00", "SSP1-1.9": "#4daf4a", "SSP3-7.0": "#d95f02",
                 "NDC (current pledges)": "#7570b3", "Control": "#555555", "Unknown": "#888888"},
    "publication": SOURCE_COLOR,
}
_FIXED_COLORS[COMPOSITE_DIM] = {
    **{k: v for k, v in _FIXED_COLORS["scenario"].items() if k not in SCENARIO_COMPOSITE},
    **{c: _FIXED_COLORS["scenario"][ssp] for ssp, c in SCENARIO_COMPOSITE.items() if ssp.startswith("SSP")},
}


def _build_hover(group, model, exp, ice_sheet, ism_m, exp_m, publication):
    extra = ""
    if "basal_melt_param" in exp_m:
        extra = f"<br>Basal melt param: {exp_m['basal_melt_param']}"
    elif "ocean_sensitivity" in exp_m:
        extra = f"<br>Ocean sensitivity: {exp_m['ocean_sensitivity']}"
    elif "mici_params" in exp_m:
        extra = f"<br>{exp_m['mici_params']}"
    return (
        f"<b>{group} / {model}</b> ({publication})<br>"
        f"Experiment: {exp}<br>"
        f"Ice model: {ism_m['ice_model']}<br>"
        f"Sliding law: {ism_m['sliding_law']}<br>"
        f"Initialization: {ism_m['initialization']}<br>"
        f"Climate model: {exp_m['climate_model']}<br>"
        f"Scenario: {exp_m['scenario']}<br>"
        f"Protocol: {exp_m['protocol']}{extra}"
    )


def _build_runs():
    """One row per simulation across every source, plus each run's raw
    (Year, cumulative Gt) arrays and an (n_runs x YEAR_GRID) matrix on
    integer years. Integer years only for the matrix: ISMIP6 AIS reports
    half-year steps where only cumulative (not rate) is populated, and
    grouping on mixed half/whole years is what broke the AIS lines in the
    notebook's projection plot."""
    sources = [(PUBLICATION_LABEL["AIS"], ismip6_ais), (PUBLICATION_LABEL["GIS"], ismip6_gis)]
    sources += [(s["label"], s["df"]) for s in EXTRA_SOURCES]
    keys = ["IS", "Group", "Model", "Exp"]
    meta, years, vals = [], [], []
    for publication, df in sources:
        df = df.loc[df[MASS_COL].notna(), keys + ["Year", MASS_COL]]
        gid = df.groupby(keys, sort=False, observed=True).ngroup().to_numpy()
        yr = df["Year"].to_numpy(dtype=float)
        mv = df[MASS_COL].to_numpy(dtype=float)
        order = np.lexsort((yr, gid))
        gid, yr, mv = gid[order], yr[order], mv[order]
        bounds = np.flatnonzero(np.diff(gid)) + 1
        firsts = np.r_[0, bounds]
        key_rows = df[keys].astype(str).to_numpy()[order][firsts]
        years += np.split(yr, bounds)
        vals += np.split(mv, bounds)
        for ice_sheet, group, model, exp in key_rows:
            ism_m = get_ism_meta(group, model)
            exp_m = get_exp_meta(ice_sheet, exp)
            meta.append({
                "ice_sheet": ice_sheet, "publication": publication,
                "group": group, "model": model, "exp": exp,
                "initialization": ism_m.get("initialization", "See paper"),
                "ice_model": ism_m.get("ice_model", "See paper"),
                "sliding_law": ism_m.get("sliding_law", "See paper"),
                "scenario": exp_m.get("scenario", "Unknown"),
                "climate_model": exp_m.get("climate_model", "Unknown"),
                "hover": _build_hover(group, model, exp, ice_sheet, ism_m, exp_m, publication),
            })
    runs = pd.DataFrame(meta)
    grid = np.arange(1950, 2102)
    mat = np.full((len(runs), len(grid)), np.nan, dtype=np.float32)
    for i, (y, v) in enumerate(zip(years, vals)):
        whole = y == np.round(y)
        idx = (y[whole] - grid[0]).astype(int)
        ok = (idx >= 0) & (idx < len(grid))
        mat[i, idx[ok]] = v[whole][ok]
    return runs, years, vals, grid, mat


RUNS, RUN_YEARS, RUN_VALS, YEAR_GRID, CUM = _build_runs()
RUNS[COMPOSITE_DIM] = RUNS["scenario"].map(lambda s: SCENARIO_COMPOSITE.get(s, s))


def _edwards_baseline_2015():
    """Edwards et al. (2021)'s emulandice series are cumulative change since
    2015 (the ISMIP6 projection start) with the first output in 2016 -- the
    archive's own results/proj_MAIN_TIMESERIES/summary_FAIR_*.csv gives a
    nonzero 2016 spread (GrIS median 0.032 cm, 5-95 % -0.021 to 0.088 cm),
    so 2016 is one year of change, not the zero point. Put each sample's
    implied 2015 value (0) on the integer-year matrix so time series and the
    2015->2100 change start from 2015 like IMBIE and ISMIP6. Rates still come
    from RUN_VALS, i.e. the published years only."""
    rows = RUNS["publication"].astype(str).str.startswith("Edwards 2021").to_numpy()
    CUM[rows, proj_start - YEAR_GRID[0]] = 0.0


_edwards_baseline_2015()

# The run table/arrays above hold everything downstream needs -- drop the
# source DataFrames so they don't sit in memory for the life of the worker
# (Render's free tier caps the whole process at 512 MB).
for _s in EXTRA_SOURCES:
    _s["df"] = None
del ismip6_ais, ismip6_gis, rahlves2025_gis, coulon2024_ais, aschwanden2022_gis, goelzer2025_gis
del edwards2021_ais, edwards2021_ais_risk, edwards2021_gis, deconto2016_ais, deconto2021_ais
import gc  # noqa: E402
gc.collect()
IS_ISMIP6 = RUNS["publication"].isin(PUBLICATION_LABEL.values()).to_numpy()


def _color_maps():
    palette = px.colors.qualitative.Dark24
    maps = {}
    for key in [k for k, _ in GROUP_DIMENSIONS if k] + [COMPOSITE_DIM]:
        fixed = _FIXED_COLORS.get(key, {})
        cmap, i = {}, 0
        for cat in sorted(RUNS[key].unique()):
            if cat in fixed:
                cmap[cat] = fixed[cat]
            else:
                cmap[cat] = palette[i % len(palette)]
                i += 1
        maps[key] = cmap
    return maps


COLOR_MAPS = _color_maps()
