"""
Standalone Dash app serving the "Interactive Figure - Rate Comparison" from
exploring_fits/analyze_slr_predictions_interactive.ipynb.

The figure is pure Plotly. Its "Distribution medians:" and "Group by:"
dropdowns are native Plotly updatemenus toggling/restyling pre-built traces,
so they run entirely client-side with no callbacks back to this server.
"Units:" is a small dcc.RadioItems + clientside_callback (see near the
bottom of this file) instead of a third native updatemenu -- a native
updatemenu can only apply a pre-baked data patch, which would mean shipping
a full second (Sea level rise) copy of every trace's x/text/hovertemplate
just for a unit conversion, roughly doubling the page's payload; the
clientside_callback does the same *2.759 conversion with plain JS on data
already in the browser instead. The "Years:" range slider and the
"Simulation studies:" checklist are different from all of the above -- the
KDEs, jitter, medians, and IMBIE slope all depend on the selected year
window and which sources are included, so those genuinely need a real Dash
callback that reruns plot_interactive_rate_comparison() on the server.

That server-side rebuild is the app's main cost: building weighted KDEs
across 7 "Group by:" dimensions x 2 panels, potentially over ~1700+ pooled
ISMIP6 + extra_sources simulations, plus JSON-serializing the ~13 MB result.
On a resource-constrained deploy target this can exceed a WSGI server's
default request timeout (gunicorn's is 30s) well before it exceeds what
actually feels "slow" locally, so raising --timeout matters (see render.yaml
alongside this file). This is specifically why the app moved off Plotly
Cloud: that platform's plotly-cloud.toml schema only accepts name/
description/app_id/app_url/team_id/team_name (checked directly against the
installed `plotly-cloud` CLI package's own AppDeploymentConfig/AppRequest
type definitions) -- no user-facing way to configure gunicorn's timeout at
all there, so a slow-but-legitimate request had no way to avoid getting
killed mid-response.

Run locally:
    pip install -r requirements.txt
    python app.py
    # -> http://127.0.0.1:8050

Deploy on Render: connect this directory's repo (see render.yaml alongside
this file, which sets `gunicorn app:server --timeout 120` as the start
command -- Render's own edge proxy allows responses up to 100 minutes, so
this app's own --timeout is genuinely the controlling limit, unlike Plotly
Cloud above).

Deploy elsewhere (any WSGI host you control, e.g. gunicorn behind nginx on
your own server):
    gunicorn app:server -b 0.0.0.0:8000 --timeout 120

Embed on an existing page once deployed:
    <iframe src="https://your-domain.example/" style="width:100%; height:900px; border:0;"></iframe>
"""

import os

import numpy as np
import pandas as pd
import scipy.stats
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.stats import linregress

import dash
from dash import dcc, html
from dash.dependencies import Input, Output

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

proj_start = 2015  # utilities/helper.py

# IMBIE community assessment mass balance time series (Otosaka et al., 2023, ESSD,
# https://doi.org/10.5194/essd-15-1597-2023), 1992-2020, from the BAS Polar Data Centre
# (https://doi.org/10.5285/77B64C55-7166-4A06-9DEF-2E400398E452).
IMBIE2023_GT_URLS = {
    "antarctica": (
        "https://ramadda.data.bas.ac.uk/repository/entry/get/imbie_antarctica_2021_Gt.csv"
        "?entryid=synth:77b64c55-7166-4a06-9def-2e400398e452:L2ltYmllX2FudGFyY3RpY2FfMjAyMV9HdC5jc3Y="
    ),
    "greenland": (
        "https://ramadda.data.bas.ac.uk/repository/entry/get/imbie_greenland_2021_Gt.csv"
        "?entryid=synth:77b64c55-7166-4a06-9def-2e400398e452:L2ltYmllX2dyZWVubGFuZF8yMDIxX0d0LmNzdg=="
    ),
}


def _load_imbie2023(region):
    """
    Loads the IMBIE 2023 community-assessment mass balance CSV for the given
    region ("antarctica" or "greenland"). Only the headline mass-balance
    columns are needed here (not the Greenland SMB/dynamics partitioning, which
    this figure doesn't use), so unlike utilities/imbie2023_loader.py this
    doesn't merge in the older dynamics workbook -- no openpyxl dependency.
    """
    df = pd.read_csv(IMBIE2023_GT_URLS[region])
    imbie = df.rename(
        columns={
            "Mass balance (Gt/yr)": "Rate of ice sheet mass change (Gt/yr)",
            "Cumulative mass balance (Gt)": "Cumulative ice sheet mass change (Gt)",
            "Cumulative mass balance uncertainty (Gt)": "Cumulative ice sheet mass change uncertainty (Gt)",
        }
    )[
        [
            "Year",
            "Cumulative ice sheet mass change (Gt)",
            "Cumulative ice sheet mass change uncertainty (Gt)",
            "Rate of ice sheet mass change (Gt/yr)",
        ]
    ].copy()

    imbie["Cumulative ice sheet mass change (Gt)"] -= imbie.loc[
        imbie["Year"] == proj_start, "Cumulative ice sheet mass change (Gt)"
    ].values

    imbie["Cumulative ice sheet mass change uncertainty (Gt)"] -= imbie[
        "Cumulative ice sheet mass change uncertainty (Gt)"
    ].values[-1]
    imbie["Cumulative ice sheet mass change uncertainty (Gt)"] *= -1

    return imbie


def load_ismip6_ais():
    """Reads the pre-generated ISMIP6 AIS scalar CSV bundled in data/."""
    return pd.read_csv(os.path.join(DATA_DIR, "ismip6_ais.csv.gz"))


def load_ismip6_gis():
    """Reads the pre-generated ISMIP6 GIS scalar CSV bundled in data/."""
    return pd.read_csv(os.path.join(DATA_DIR, "ismip6_gis_ctrl.csv.gz"))


# ─────────────────────────────────────────────────────────────────────────────
# Experiment metadata look-up tables
# Sources:
#   GIS – Goelzer et al. (2020) The Cryosphere 14, 3071-3096
#   AIS – Seroussi et al. (2020) The Cryosphere 14, 3033-3070
# ─────────────────────────────────────────────────────────────────────────────

# GIS experiment descriptions (Goelzer et al. 2020, Table 1 & Appendix)
gis_exp_meta = {
    # "exp01":  {"climate_model": "MIROC5",       "scenario": "RCP 8.5", "protocol": "Open",     "ocean_sensitivity": "Medium"},
    # "exp02":  {"climate_model": "MIROC5",       "scenario": "RCP 8.5", "protocol": "Open",     "ocean_sensitivity": "Low"},
    # "exp03":  {"climate_model": "MIROC5",       "scenario": "RCP 2.6", "protocol": "Open",     "ocean_sensitivity": "Medium"},
    # "exp04":  {"climate_model": "MIROC5",       "scenario": "RCP 2.6", "protocol": "Open",     "ocean_sensitivity": "Low"},
    "exp05":  {"climate_model": "MIROC5",       "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp06":  {"climate_model": "NorESM",       "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp07":  {"climate_model": "MIROC5",       "scenario": "RCP 2.6", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp08":  {"climate_model": "HadGEM2-ES",   "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "exp09":  {"climate_model": "MIROC5",       "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "High"},
    "exp10":  {"climate_model": "MIROC5",       "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Low"},
    # "exp11":  {"climate_model": "ACCESS1.3",    "scenario": "RCP 8.5", "protocol": "Open",     "ocean_sensitivity": "Medium"},
    # "exp12":  {"climate_model": "ACCESS1.3",    "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    # "exp13":  {"climate_model": "CESM2",        "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "High"},
    "expa01": {"climate_model": "IPSL-CM5A-MR", "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "expa02": {"climate_model": "CSIRO-Mk3.6",  "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
    "expa03": {"climate_model": "ACCESS1.3",    "scenario": "RCP 8.5", "protocol": "Standard", "ocean_sensitivity": "Medium"},
}

# AIS experiment descriptions (Seroussi et al. 2020, Table 1)
ais_exp_meta = {
    "exp01": {"climate_model": "NorESM",         "scenario": "RCP 8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp02": {"climate_model": "MIROC-ESM-CHEM", "scenario": "RCP 8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp03": {"climate_model": "NorESM",         "scenario": "RCP 2.6",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp04": {"climate_model": "CCSM4",          "scenario": "RCP 8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp05": {"climate_model": "NorESM",         "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp06": {"climate_model": "MIROC-ESM-CHEM", "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp07": {"climate_model": "NorESM",         "scenario": "RCP 2.6",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp08": {"climate_model": "CCSM4",          "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp09": {"climate_model": "NorESM",         "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "PIGL medium"},
    "exp10": {"climate_model": "NorESM",         "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "PIGL high"},
    "exp11": {"climate_model": "CCSM4",          "scenario": "RCP 8.5",  "protocol": "Open",     "basal_melt_param": "Standard"},
    "exp12": {"climate_model": "CCSM4",          "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "Standard"},
    "exp13": {"climate_model": "NorESM",         "scenario": "RCP 8.5",  "protocol": "Standard", "basal_melt_param": "PIGL very high"},
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
    # AIS additional groups
    ("ARC",      "PISM1"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("ARC",      "PISM2"):     {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("DOE",      "MALI"):      {"ice_model": "MALI",      "sliding_law": "Coulomb",                  "initialization": "Data assimilation"},
    ("GRL",      "PISM"):      {"ice_model": "PISM",      "sliding_law": "Pseudo-plastic",           "initialization": "Spin-up"},
    ("GSFC",     "ISSM"):      {"ice_model": "ISSM",      "sliding_law": "Weertman",                 "initialization": "Data assimilation"},
    ("IGE",      "ElmerIce"):  {"ice_model": "Elmer/Ice", "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
    ("ILTS",     "SICOPOLIS"): {"ice_model": "SICOPOLIS", "sliding_law": "Weertman",                 "initialization": "Spin-up"},
    ("NEMO",     "fETISh"):    {"ice_model": "Kori-ULB/f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
    ("PSU",      "PSU3D1"):    {"ice_model": "PSU-ISM",   "sliding_law": "Coulomb / Weertman",       "initialization": "Spin-up"},
    ("PSU",      "PSU3D2"):    {"ice_model": "PSU-ISM",   "sliding_law": "Coulomb / Weertman",       "initialization": "Spin-up"},
    ("UTAS",     "ElmerIce"):  {"ice_model": "Elmer/Ice", "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},

    # Non-ISMIP6 published simulations (see utilities/external_sources.py for
    # provenance/derivation of each).
    ("Rahlves2025", "CISM"):     {"ice_model": "CISM",      "sliding_law": "Weertman (power-law)",        "initialization": "ERA5- or ESM-forced spin-up"},
    ("Coulon2024",  "Kori-ULB"): {"ice_model": "Kori-ULB/f.ETISh",  "sliding_law": "See paper",                   "initialization": "Nudged to present-day geometry"},
    ("Aschwanden2022", "PISM"):  {"ice_model": "PISM",      "sliding_law": "See paper",                   "initialization": "Temperature-index SMB, present-day start"},
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

def _exp_meta_from_df(df, extra_cols):
    cols = ["Exp", "climate_model", "scenario", "protocol"] + extra_cols
    unique = df[cols].drop_duplicates("Exp").set_index("Exp")
    return unique.to_dict(orient="index")


rahlves2025_gis = pd.read_csv(os.path.join(DATA_DIR, "external_sources_rahlves2025_gis.csv.gz"))
coulon2024_ais = pd.read_csv(os.path.join(DATA_DIR, "external_sources_coulon2024_ais.csv.gz"))
aschwanden2022_gis = pd.read_csv(os.path.join(DATA_DIR, "external_sources_aschwanden2022_gis.csv.gz"))
gis_exp_meta.update(_exp_meta_from_df(rahlves2025_gis, ["ocean_sensitivity"]))
ais_exp_meta.update(_exp_meta_from_df(coulon2024_ais, ["basal_melt_param"]))
gis_exp_meta.update(_exp_meta_from_df(aschwanden2022_gis, []))

EXTRA_SOURCES = [
    {"label": "Rahlves 2025", "df": rahlves2025_gis, "color": "#e6550d"},
    {"label": "Coulon 2024", "df": coulon2024_ais, "color": "#31a354"},
    {"label": "Aschwanden 2022", "df": aschwanden2022_gis, "color": "#756bb1"},
]

# plot_interactive_rate_comparison's "Group by publication" category names/
# color for each panel's own ISMIP6 population (defined here, ahead of that
# function, so the checklist labels below can share the same strings).
PUBLICATION_LABEL = {"AIS": "Seroussi 2020 (ISMIP6 AIS)", "GIS": "Goelzer 2020 (ISMIP6 GIS)"}
ISMIP6_PUBLICATION_COLOR = "#636363"  # ISMIP6's own color under "Group by publication" -- distinct from unify-gray and from each extra_sources paper's own color

# The two core ISMIP6 panels (each *is* one whole precomputed_rows["AIS"/"GIS"]
# entry -- see _update_figure) get checkbox labels alongside the extra_sources
# papers above, named after each panel's own source paper (Seroussi et al.
# 2020 / Goelzer et al. 2020) rather than "ISMIP6 AIS"/"ISMIP6 GIS", to read
# as one flat list of papers rather than singling ISMIP6 out as different in
# kind from Rahlves2025/Coulon2024/Aschwanden2022. Sourced from
# PUBLICATION_LABEL (plot_interactive_rate_comparison's own "Group by
# publication" category names) so the checklist and the figure never drift
# apart on what to call each panel's ISMIP6 population.
ISMIP6_AIS_LABEL = PUBLICATION_LABEL["AIS"]
ISMIP6_GIS_LABEL = PUBLICATION_LABEL["GIS"]
DATA_SOURCE_LABELS = [src["label"] for src in EXTRA_SOURCES] + [ISMIP6_AIS_LABEL, ISMIP6_GIS_LABEL]


# ─────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

RCP_COLOR = {
    "RCP 2.6":  "#003466",
    "RCP 8.5":  "#990002",
    "SSP1-2.6": "#1a7f5e",
    "SSP5-8.5": "#8b1a00",
    "Unknown":  "#888888",
}
IMBIE_COLOR = "#08519c"

gt2cmSLE = 1.0 / 362.5 / 10.0  # Gt -> cm sea-level-equivalent
gt2mmSLE = gt2cmSLE * 10  # the "Units:" toggle uses mm, not cm -- typical rates are well under 1 cm/yr

INIT_COLOR = {
    "Data assimilation": "#1a7f5e",
    "Spin-up": "#003466",
    "See paper": "#6a6a6a",
}

rate_strip_y0, rate_strip_jitter, rate_kde_y0, rate_kde_height = 0.0, 0.15, 0.30, 0.55
rate_y_lim = (-0.35, 1.05)
rate_median_y = rate_kde_y0 + rate_kde_height + 0.08  # above the tallest possible KDE curve

# Each entry is (dropdown label, sim_df column to group by, fixed color map or
# None to assign colors dynamically). "All simulations" (dim=None) is the
# ungrouped gray default.
GROUP_DIMENSIONS = [
    ("All simulations", None, None),
    ("Group by initialization", "initialization", INIT_COLOR),
    ("Group by ice sheet model", "ice_model", None),
    ("Group by sliding law", "sliding_law", None),
    ("Group by climate scenario", "scenario", RCP_COLOR),
    ("Group by GCM", "climate_model", None),
]


def _build_hover(row_group, row_model, row_exp, ice_sheet):
    """Compose the hover-text shown for a single simulation trace."""
    exp_m = get_exp_meta(ice_sheet, row_exp)
    ism_m = get_ism_meta(row_group, row_model)

    extra = ""
    if ice_sheet == "AIS" and "basal_melt_param" in exp_m:
        extra = f"<br>Basal melt param : {exp_m['basal_melt_param']}"
    if ice_sheet == "GIS" and "ocean_sensitivity" in exp_m:
        extra = f"<br>Ocean sensitivity : {exp_m['ocean_sensitivity']}"

    return (
        f"<b>{row_group} / {row_model}</b><br>"
        f"Experiment : {row_exp}<br>"
        f"Ice model  : {ism_m['ice_model']}<br>"
        f"Sliding law: {ism_m['sliding_law']}<br>"
        f"Init. method: {ism_m['initialization']}<br>"
        f"Climate model: {exp_m['climate_model']}<br>"
        f"Scenario   : {exp_m['scenario']}<br>"
        f"Protocol   : {exp_m['protocol']}"
        f"{extra}"
    )


def imbie_mass_loss_slope(df, year_start=2000, year_end=2025):
    """Linear-regression slope (Gt/yr) of observed cumulative ice sheet mass change."""
    mask = (
        (df["Year"] >= year_start) & (df["Year"] <= year_end)
        & df["Cumulative ice sheet mass change (Gt)"].notna()
    )
    x = df.loc[mask, "Year"]
    y = df.loc[mask, "Cumulative ice sheet mass change (Gt)"]
    result = linregress(x, y)
    return result.slope, result.stderr


def _categorical_color_map(categories, fixed=None, palette=None):
    """Assigns a color to each category, reusing `fixed` entries where given
    and cycling through a qualitative palette for the rest."""
    palette = palette or px.colors.qualitative.Dark24
    color_map = dict(fixed or {})
    i = 0
    for cat in categories:
        if cat in color_map:
            continue
        color_map[cat] = palette[i % len(palette)]
        i += 1
    return color_map


def _rate_kde_raw(values, x_range, n=200, weights=None):
    """Raw KDE curve (integrates to 1 over its own support) -- NOT renormalized
    to its own peak. Callers weight and rescale this against a shared reference
    so peak heights stay comparable across subgroups of very different spread.

    `weights` lets each point in `values` contribute unequally to the
    estimate (passed straight through to gaussian_kde) -- used to blend
    ISMIP6 points (weight 1) with extra_sources points (weight
    typical_model_n/n_src, see plot_interactive_rate_comparison) into a
    single properly-proportioned curve, rather than computing two separate
    curves and summing them.

    Falls back to an all-zero (flat) curve if the covariance is singular --
    e.g. a narrow year window can leave a subgroup with near-but-not-exactly-
    zero variance, which passes a `std() > 0` guard but still isn't enough
    for gaussian_kde's Cholesky step."""
    xs = np.linspace(x_range[0], x_range[1], n)
    try:
        dens = scipy.stats.gaussian_kde(values, weights=weights)(xs)
    except np.linalg.LinAlgError:
        dens = np.zeros_like(xs)
    return xs, dens


def _rgba(color, alpha):
    """Bakes an alpha channel directly into an rgba(...) string. Plotly's
    trace-level `opacity` does not reliably apply to `fill` (a lone fill at
    opacity=0.14 still renders fully solid), so every translucent fill in
    this figure uses this instead of the `opacity` kwarg."""
    if color in ("gray", "grey"):
        r, g, b = 128, 128, 128
    else:
        c = color.lstrip("#")
        r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
    return f"rgba({r}, {g}, {b}, {alpha})"


def _dim_color_maps(simulated, extra_sources=None):
    """Category -> color, per grouping dimension, built from the FULL
    simulated dataframe AND every extra_sources entry (not windowed by
    year) so a category's color stays consistent across both panels, every
    year-range selection, and regardless of which papers are currently
    checked in the caller's own UI (pass the caller's *complete* list here,
    not a checkbox-filtered subset, so toggling a paper on/off never
    reassigns another paper's or ISMIP6's colors). Factored out of
    plot_interactive_rate_comparison so a caller (e.g. the Dash app) can
    compute this once and pass it back in via the
    `precomputed_dim_color_maps` kwarg, instead of paying its groupby cost
    on every request -- it's the same regardless of year_start/year_end."""
    dims = [d for _, d, _ in GROUP_DIMENSIONS if d is not None]
    dim_color_maps = {}
    dim_cats = {d: set() for d in dims}

    def _collect(df):
        if df is None:
            return
        for (ice_sheet_, group, model, exp), _ in df.groupby(["IS", "Group", "Model", "Exp"]):
            ism_m = get_ism_meta(group, model)
            exp_m = get_exp_meta(ice_sheet_, exp)
            row_meta = {**ism_m, **exp_m}
            for d in dims:
                dim_cats[d].add(row_meta.get(d, "See paper"))

    _collect(simulated)
    for src in (extra_sources or []):
        _collect(src["df"])

    for label, dim, fixed in GROUP_DIMENSIONS:
        if dim is None:
            continue
        dim_color_maps[dim] = _categorical_color_map(sorted(dim_cats[dim]), fixed=fixed)
    return dim_color_maps


def _sim_rows(simulated, ice_sheet, year_start, year_end):
    """Per-simulation rate + metadata rows for one ice sheet and year window
    -- the pandas-groupby-plus-linregress-per-simulation step that has to
    rerun for every distinct (year_start, year_end), factored out of
    plot_interactive_rate_comparison so a caller can precompute rows for
    many year windows up front (e.g. the Dash app precomputes every valid
    "Years:" slider position at startup) rather than redoing this specific
    step -- profiling showed it's roughly a third of a single call's total
    time -- on every request. The remaining trace/KDE-building work still
    has to happen per-request either way."""
    rows = []
    if simulated is None:
        return rows
    df = simulated[simulated["IS"] == ice_sheet]
    for (group, model, exp), g in df.groupby(["Group", "Model", "Exp"]):
        n_pts = (
            (g["Year"] >= year_start) & (g["Year"] <= year_end)
            & g["Cumulative ice sheet mass change (Gt)"].notna()
        ).sum()
        if n_pts < 2:
            continue
        slope, _ = imbie_mass_loss_slope(g, year_start=year_start, year_end=year_end)
        if slope is None or not np.isfinite(slope):
            continue
        ism_m = get_ism_meta(group, model)
        exp_m = get_exp_meta(ice_sheet, exp)
        base_hover = _build_hover(group, model, exp, ice_sheet)
        rows.append({
            "rate": slope, "group": group, "model": model, "exp": exp,
            "initialization": ism_m.get("initialization", "See paper"),
            "ice_model": ism_m.get("ice_model", "See paper"),
            "sliding_law": ism_m.get("sliding_law", "See paper"),
            "scenario": exp_m.get("scenario", "Unknown"),
            "climate_model": exp_m.get("climate_model", "Unknown"),
            "hover": base_hover + f"<br>Rate: {slope:.0f} Gt/yr",
            "hover_sle": base_hover + f"<br>Rate: {slope * gt2mmSLE:.2f} mm/yr",
        })
    return rows


GROUP_BY_LABELS = [label for label, _dim, _fixed in GROUP_DIMENSIONS] + ["Group by publication"]


def plot_interactive_rate_comparison(
    simulated=None, observed=None, year_start=2000, year_end=2020,
    title="Compare observed and simulated rates of ice sheet change", seed=42,
    show_title=True, show_subtitle=True,
    precomputed_dim_color_maps=None, precomputed_rows=None,
    extra_sources=None,
):
    """
    Interactive dual-panel (AIS / GIS) raincloud comparison of the mean
    ice-sheet mass-change rate, matching the static plot_rate_comparison
    figure's default look, plus two independent controls:

    - "Group by:" dropdown -- every option pools ISMIP6 (`simulated`) AND
      every entry in `extra_sources` into one population, then either
      leaves it ungrouped ("All simulations": gray jittered points + gray
      KDE) or recolors/splits it by a category (initialization method, ice
      sheet model, sliding law, climate scenario, GCM, or -- "Group by
      publication" -- which paper the simulation is from at all: ISMIP6
      counts as one category per panel, Seroussi et al. 2020 for AIS /
      Goelzer et al. 2020 for GIS, alongside one category per
      extra_sources entry). Every dimension shares the same machinery
      (see GROUP_DIMENSIONS/GROUP_BY_LABELS): each category's KDE is a
      properly weighted density estimate (scipy.stats.gaussian_kde's
      `weights`, via _rate_kde_raw) over however many ISMIP6 and/or
      extra_sources points fall into it, scaled by that category's share
      of the total weighted population and referenced against the SAME
      peak (the full pooled population's) as every other category -- not
      its own peak. This keeps the AREA under each colored curve
      proportional to its share of simulations, so a widely-scattered
      subgroup reads as low and wide rather than being inflated to look as
      prominent as any other category, and a tightly-clustered subgroup can
      legitimately show as a sharp, narrow spike even with relatively few
      points. Summing every category's curve approximately reconstructs
      the "All simulations" gray curve. Because an extra_sources entry can
      have far more points than a typical ISMIP6 institution (e.g. a
      100-member parameter ensemble), each of ITS points individually
      counts for only `typical_model_n / n_src` of an ISMIP6 point's
      weight (`typical_model_n` = ISMIP6's own average points-per-Group for
      this panel, `n_src` = that source's own total point count here) --
      so its ensemble collectively counts for about as much as one typical
      ISMIP6 institution, in every dimension, not just "Group by
      publication" -- otherwise one paper's internal ensemble size would
      visually dominate ISMIP6's ~30-40 independent models regardless of
      which dimension is selected. Point markers are similarly
      shrunk/faded for large-N sources.
    - "Medians:" dropdown -- independently shows/hides a downward-pointing
      triangle marker at each currently-visible category's median rate,
      drawn above the KDE curves. Implemented as a targeted `restyle` on
      trace `opacity` (not `visible`) for just the median-marker traces, so
      it never conflicts with the "Group by:" dropdown's own visibility
      updates -- the two controls are fully independent.

    Each simulation's vertical jitter offset is computed once and reused
    across every grouping dimension, so points don't jump up/down when the
    "Group by:" dropdown is toggled. IMBIE's mean-rate line and 2-sigma band
    are drawn as real traces (not shapes) so they get proper legend entries,
    and stay visible regardless of the dropdown state. A category's legend
    entry is shown the first time it appears (AIS is processed first), so a
    category that only exists in the GIS panel (e.g. an experiment code with
    no scenario/GCM metadata, falling back to "Unknown") still gets a legend
    entry instead of appearing as an unlabeled curve.

    show_title/show_subtitle -- set False to omit the in-figure title and/or
    explanatory-text annotation (e.g. the Dash app renders both itself as
    plain HTML instead, alongside rather than above the figure, so it passes
    both as False). Omitting either also shrinks the top margin back down,
    since that space exists only to hold whichever of the two is shown.

    precomputed_dim_color_maps/precomputed_rows -- optional outputs of
    _dim_color_maps()/_sim_rows(), letting a caller that already computed
    these (e.g. once at startup, or once per year-range ahead of time) skip
    redoing that work on every call. Both fall back to computing from
    `simulated`/`extra_sources` as usual when omitted. If passing a
    precomputed value, build it from the caller's *complete* extra_sources
    list (see _dim_color_maps), not a checkbox-filtered subset.

    extra_sources -- optional list of {"label", "df", "color"} dicts for
    non-ISMIP6 published simulations (see utilities/external_sources.py).
    Which papers are included is entirely up to the caller -- there is no
    separate on/off control inside the figure itself (the Dash app's
    "Simulation studies:" checkboxes filter this list before calling here;
    the notebook always passes all of them). Each source's `df` must have
    the same shape as `simulated` (Year, Cumulative ice sheet mass change
    (Gt), Group, Model, Exp, IS) -- rows for an ice sheet a source doesn't
    cover (e.g. Rahlves2025 is GIS-only) are simply absent, so that panel is
    skipped for it.
    """
    rng = np.random.default_rng(seed)
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=False,
        subplot_titles=["Antarctic Ice Sheet (AIS)", "Greenland Ice Sheet (GIS)"],
        vertical_spacing=0.20,
    )

    # Global color maps (shared between AIS/GIS panels) per grouping dimension,
    # built from the full simulated dataframe (+ extra_sources) so a
    # category's color is consistent across both panels.
    dim_color_maps = (
        precomputed_dim_color_maps if precomputed_dim_color_maps is not None
        else _dim_color_maps(simulated, extra_sources=extra_sources)
    )

    dim_only_idx = {label: [] for label, dim, _ in GROUP_DIMENSIONS if dim is not None}
    all_only_idx = []  # "All simulations": ISMIP6 + every extra_sources entry, pooled and gray
    publication_only_idx = []  # "Group by publication": one category per paper (ISMIP6 counts as one)
    median_trace_idx = []  # every median-marker trace, across all dimensions -- toggled by the "Medians:" dropdown
    legend_shown = set()  # (dim, cat) pairs already given a legend entry, across both panels
    x_range_by_panel = {}  # row (1=AIS, 2=GIS) -> (min, max), for the Units: dropdown's explicit axis ranges
    observed_annotation_x = {}  # annotation index -> its Gt/yr x-position (slope), for the Units: dropdown

    for k, ice_sheet in enumerate(["AIS", "GIS"], start=1):
        if precomputed_rows is not None and ice_sheet in precomputed_rows:
            rows = precomputed_rows[ice_sheet]
        else:
            rows = _sim_rows(simulated, ice_sheet, year_start, year_end)
        sim_df = pd.DataFrame(rows)
        n_total = len(sim_df)
        # Reference point count for de-weighting extra_sources below --
        # ISMIP6's own average points-per-institution for this panel. With
        # no ISMIP6 rows in this panel (e.g. that panel's checkbox is
        # unchecked in the Dash app), there's no such reference, so
        # extra_sources fall back to their natural (unweighted) contribution.
        typical_model_n = max(1, n_total / max(1, sim_df["group"].nunique())) if n_total else None

        # One merged population: every ISMIP6 row (weight 1) + every
        # extra_sources row (weight typical_model_n/n_src, see docstring),
        # each tagged with its own "publication" -- lets "Group by
        # publication" reuse the exact same generic per-dimension machinery
        # every other "Group by ..." option below uses, and lets every
        # OTHER dimension (ice model, scenario, ...) include extra_sources'
        # rows too rather than being ISMIP6-only.
        frames = []
        publication_color_map = {}
        if n_total:
            sim_df["jitter"] = rate_strip_y0 + rng.uniform(-rate_strip_jitter, rate_strip_jitter, size=n_total)
            sim_df["row_weight"] = 1.0
            sim_df["publication"] = PUBLICATION_LABEL[ice_sheet]
            frames.append(sim_df)
            publication_color_map[PUBLICATION_LABEL[ice_sheet]] = ISMIP6_PUBLICATION_COLOR
        for src in (extra_sources or []):
            src_rows = _sim_rows(src["df"], ice_sheet, year_start, year_end)
            if not src_rows:
                continue
            src_df = pd.DataFrame(src_rows)
            n_src = len(src_df)
            src_df["jitter"] = rate_strip_y0 + rng.uniform(-rate_strip_jitter, rate_strip_jitter, size=n_src)
            src_df["row_weight"] = (typical_model_n / n_src) if typical_model_n else 1.0
            src_df["publication"] = src["label"]
            frames.append(src_df)
            publication_color_map[src["label"]] = src.get("color", "#888888")

        if frames:
            merged_df = pd.concat(frames, ignore_index=True)

            pad = max(0.15 * (merged_df["rate"].max() - merged_df["rate"].min()), 1.0)
            x_range = (merged_df["rate"].min() - pad, merged_df["rate"].max() + pad)
            x_range_by_panel[k] = x_range

            total_weight = merged_df["row_weight"].sum()
            xs, dens_full_raw = _rate_kde_raw(
                merged_df["rate"].values, x_range, weights=merged_df["row_weight"].values,
            )
            full_max = dens_full_raw.max()
            if full_max <= 0:
                # Degenerate window (e.g. all simulations converged to the
                # same rate) -- avoid a division by zero; nothing meaningful
                # to scale against, so every curve stays flat.
                full_max = 1.0
            dens = dens_full_raw / full_max * rate_kde_height

            # "All simulations" -- ISMIP6 + every extra_sources entry, pooled
            # and gray. By construction this curve always peaks at exactly
            # rate_kde_height (it IS the full_max reference), matching how
            # the ISMIP6-only version worked before extra_sources existed.
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=xs, y=np.full_like(xs, rate_kde_y0), mode="lines", line=dict(width=0),
                hoverinfo="skip", showlegend=False,
            ), row=k, col=1)
            all_only_idx.append(idx)
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=xs, y=rate_kde_y0 + dens, mode="lines", line=dict(color="gray", width=1),
                fill="tonexty", fillcolor=_rgba("gray", 0.5), hoverinfo="skip",
                name="PDF of all simulations", legendgroup="kde_all", showlegend=(k == 1),
            ), row=k, col=1)
            all_only_idx.append(idx)

            # Gray jittered points -- every row in the merged population.
            # ISMIP6 rows render at the usual fixed size/opacity; each
            # extra_sources row is shrunk/faded by its OWN source's n_src
            # (same de-emphasis every other dimension below applies to that
            # source), so a large ensemble reads as a density cloud here
            # too, not a wall of dots outnumbering ISMIP6's own.
            size_by_pub, opacity_by_pub = {}, {}
            if n_total:
                size_by_pub[PUBLICATION_LABEL[ice_sheet]] = 4
                opacity_by_pub[PUBLICATION_LABEL[ice_sheet]] = 0.6
            for src in (extra_sources or []):
                n_src = (merged_df["publication"] == src["label"]).sum()
                if n_src == 0:
                    continue
                size_by_pub[src["label"]] = 3 if n_src > 30 else 4
                opacity_by_pub[src["label"]] = max(0.12, min(0.6, 15 / n_src))
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=merged_df["rate"], y=merged_df["jitter"], mode="markers",
                marker=dict(
                    color="gray", line_width=0,
                    size=merged_df["publication"].map(size_by_pub).tolist(),
                    opacity=merged_df["publication"].map(opacity_by_pub).tolist(),
                ),
                name="Simulation", legendgroup="all_pts", showlegend=(k == 1), visible=True,
                text=merged_df["hover"], customdata=merged_df["hover_sle"], hovertemplate="%{text}<extra></extra>",
            ), row=k, col=1)
            all_only_idx.append(idx)

            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=[merged_df["rate"].median()], y=[rate_median_y], mode="markers",
                marker=dict(symbol="triangle-down", size=11, color="gray", line=dict(width=1, color="black")),
                opacity=1, name="All simulations median", legendgroup="median:all", showlegend=(k == 1),
                visible=True, hovertemplate="All simulations median<br>" + "Rate: %{x:.0f} Gt/yr<extra></extra>",
            ), row=k, col=1)
            all_only_idx.append(idx)
            median_trace_idx.append(idx)

            # Colored KDE + colored jittered points + median marker per
            # category, for each grouping dimension (including the new
            # "publication" one), drawn on top of the gray KDE -- every
            # dimension shares this exact loop body, over the same merged
            # population, differing only in which column it groups by and
            # which color map/target trace-index list it uses.
            all_dims = [
                (label, dim, dim_color_maps[dim], dim_only_idx[label])
                for label, dim, _fixed in GROUP_DIMENSIONS if dim is not None
            ]
            all_dims.append(("Group by publication", "publication", publication_color_map, publication_only_idx))

            for label, dim, color_map, target_idx in all_dims:
                # Semi-transparent fills compound when many categories overlap in the
                # same x-range (n stacked layers at opacity a look like 1-(1-a)^n, so
                # even a=0.5 reads as ~90% opaque by the 4th overlapping layer). Scale
                # each category's fill opacity down for dimensions with more
                # categories (e.g. ice sheet model, ~14) so the *stacked* result stays
                # translucent, while dimensions with few categories (e.g.
                # initialization, ~4) keep the full 0.5.
                n_cats = max(len(color_map), 1)
                kde_fill_opacity = min(0.5, 2.0 / n_cats)
                for cat, g in merged_df.groupby(dim):
                    color = color_map.get(cat, "#888888")
                    cat_weight = g["row_weight"].sum()
                    if len(g) >= 2 and g["rate"].std() > 0:
                        _, dens2_raw = _rate_kde_raw(g["rate"].values, x_range, weights=g["row_weight"].values)
                        weight = cat_weight / total_weight
                        dens2 = (dens2_raw * weight) / full_max * rate_kde_height
                        idx = len(fig.data)
                        fig.add_trace(go.Scatter(
                            x=xs, y=np.full_like(xs, rate_kde_y0), mode="lines", line=dict(width=0),
                            hoverinfo="skip", showlegend=False, visible=False,
                        ), row=k, col=1)
                        target_idx.append(idx)
                        idx = len(fig.data)
                        fig.add_trace(go.Scatter(
                            x=xs, y=rate_kde_y0 + dens2, mode="lines", line=dict(color=color, width=1),
                            fill="tonexty", fillcolor=_rgba(color, kde_fill_opacity), hoverinfo="skip",
                            name=f"{cat} KDE", legendgroup=f"{dim}:{cat}", showlegend=False, visible=False,
                        ), row=k, col=1)
                        target_idx.append(idx)

                    show_this = (dim, cat) not in legend_shown
                    if show_this:
                        legend_shown.add((dim, cat))
                    idx = len(fig.data)
                    fig.add_trace(go.Scatter(
                        x=g["rate"], y=g["jitter"], mode="markers",
                        marker=dict(color=color, size=4, opacity=0.85, line_width=0),
                        name=str(cat), legendgroup=f"{dim}:{cat}", showlegend=show_this, visible=False,
                        text=g["hover"], customdata=g["hover_sle"], hovertemplate="%{text}<extra></extra>",
                    ), row=k, col=1)
                    target_idx.append(idx)

                    # Median marker: visibility tied to the grouping dropdown
                    # (like its KDE/points siblings), but opacity toggled
                    # independently by the "Medians:" dropdown.
                    idx = len(fig.data)
                    fig.add_trace(go.Scatter(
                        x=[g["rate"].median()], y=[rate_median_y], mode="markers",
                        marker=dict(symbol="triangle-down", size=11, color=color, line=dict(width=1, color="black")),
                        opacity=1, name=f"{cat} median", legendgroup=f"median:{dim}:{cat}", showlegend=False,
                        visible=False, hovertemplate=f"{cat} median<br>" + "Rate: %{x:.0f} Gt/yr<extra></extra>",
                    ), row=k, col=1)
                    target_idx.append(idx)
                    median_trace_idx.append(idx)

        if observed is not None:
            imbie_df = observed[observed["IS"] == ice_sheet]
            slope, stderr = imbie_mass_loss_slope(imbie_df, year_start=year_start, year_end=year_end)

            # Real traces (not add_vline/add_vrect shapes) so IMBIE gets legend
            # entries, matching the other interactive figures' convention.
            fig.add_trace(go.Scatter(
                x=[slope - 2 * stderr, slope - 2 * stderr, slope + 2 * stderr, slope + 2 * stderr],
                y=[rate_y_lim[0], rate_y_lim[1], rate_y_lim[1], rate_y_lim[0]],
                mode="lines", line=dict(width=0), fill="toself", fillcolor=_rgba(IMBIE_COLOR, 0.2),
                hoverinfo="skip",
                name="IMBIE ±2σ", legendgroup="imbie", showlegend=(k == 1), legendrank=9999,
            ), row=k, col=1)
            fig.add_trace(go.Scatter(
                x=[slope, slope], y=[rate_y_lim[0], rate_y_lim[1]],
                mode="lines", line=dict(color=IMBIE_COLOR, width=2), hoverinfo="skip",
                name="IMBIE (observed)", legendgroup="imbie", showlegend=(k == 1), legendrank=10000,
            ), row=k, col=1)
            fig.add_annotation(
                x=slope, y=0.55, xref=f"x{k}", yref=f"y{k}",
                text="<b>Observed mass change</b>",
                showarrow=True, arrowhead=2, arrowsize=1.1, arrowwidth=2.5,
                arrowcolor=IMBIE_COLOR,
                ax=60, ay=-50,
                font=dict(size=11, color=IMBIE_COLOR),
                bgcolor="rgba(255,255,255,0.85)",
            )
            observed_annotation_x[len(fig.layout.annotations) - 1] = slope

        fig.add_vline(x=0, line_dash="dot", line_color="black", line_width=0.8, row=k, col=1)
        fig.update_xaxes(
            title_text="Rate of ice sheet mass change (Gt/yr)",
            range=list(x_range_by_panel[k]) if k in x_range_by_panel else None,
            row=k, col=1,
        )
        fig.update_yaxes(showticklabels=False, range=list(rate_y_lim), row=k, col=1)

    n_traces = len(fig.data)
    # "Group by:" buttons manage every trace built above (dim_only_idx,
    # all_only_idx, publication_only_idx) -- exactly one of these groups is
    # visible at a time. Traces outside all of these (IMBIE) are left alone.
    managed_groups = dict(dim_only_idx)
    managed_groups["All simulations"] = all_only_idx
    managed_groups["Group by publication"] = publication_only_idx
    managed_idx = sum(managed_groups.values(), [])
    buttons = []
    for label in GROUP_BY_LABELS:
        visible = [i in managed_groups[label] for i in managed_idx]
        # method="restyle" (not "update"): this only ever touches trace data
        # (visible), never layout, and restyle's args signature is exactly
        # [dataUpdate, traceIndices] -- "update"'s signature is instead
        # [dataUpdate, layoutUpdate, traceIndices], so passing managed_idx as
        # the 2nd positional arg to "update" put a list of integers where
        # Plotly expects a layout-patch object, silently breaking which
        # traces the visibility patch actually applied to (this is what
        # caused a real bug: empty panels and stray legend entries from
        # non-selected dimensions after clicking "Group by:").
        buttons.append(dict(label=label, method="restyle", args=[{"visible": visible}, managed_idx]))

    # "Units:" is NOT an in-figure Plotly updatemenu here (unlike the
    # notebook copy of this function) -- a native Plotly updatemenu can only
    # apply a static, pre-baked data patch, which would mean shipping a full
    # second (Sea level rise) copy of every trace's x/text/hovertemplate to
    # the browser just for a *2.759 unit conversion, roughly doubling the
    # page's data payload. Instead, the Dash app wires "Units:" up as its own
    # dcc.RadioItems (see app.layout) plus a clientside_callback (see below
    # this function) that does the conversion with plain JS math on data
    # already sitting in the browser -- no server round trip, and no
    # duplicated payload. x_range_by_panel/observed_annotation_x below feed
    # that clientside callback's initial-state expectations (it re-derives
    # everything from the CURRENTLY-rendered mass-Gt figure, so it works
    # whether this FIG came from initial load or `_update_figure`).
    layout_kwargs = dict(
        updatemenus=[
            dict(
                type="dropdown", direction="down",
                x=1.0, y=1.13, xanchor="right", yanchor="bottom",
                buttons=buttons,
            ),
            dict(
                type="dropdown", direction="down",
                x=0.48, y=1.13, xanchor="left", yanchor="bottom",
                active=1,
                buttons=[
                    dict(label="Off", method="restyle", args=[{"opacity": 0}, median_trace_idx]),
                    dict(label="On", method="restyle", args=[{"opacity": 1}, median_trace_idx]),
                ],
            ),
        ],
        margin=dict(t=300 if show_subtitle else 130),
        height=750, template="plotly_white",
        legend=dict(bgcolor="rgba(255,255,255,0.85)", bordercolor="#ccc", borderwidth=1),
    )
    if show_title:
        layout_kwargs["title"] = dict(text=title, x=0.02, xanchor="left", y=0.99, yanchor="top")
    fig.update_layout(**layout_kwargs)
    # add_annotation appends to the existing (subplot-title) annotations list,
    # rather than replacing it the way update_layout(annotations=[...]) would.
    if show_subtitle:
        fig.add_annotation(
            text=(
                "PDFs show the distribution of simulated mass changes over the selected time span.<br>"
                "PDFs are calculated collectively for all simulations, or (optionally) grouped by a<br>"
                "variety of different simulation characteristics.<br>"
                "Beneath the PDFs, each dot represents a single simulation, plotted at its simulated<br>"
                "rate of mass change.<br>"
                "Hover over each dot to learn more about its characteristics, or click on legend<br>"
                "elements to hide or show plot components."
            ),
            x=0.0, y=1.55, xref="paper", yref="paper",
            xanchor="left", yanchor="top", showarrow=False,
            align="left", font=dict(size=11, color="#555555"),
        )
    fig.add_annotation(
        text="Group by:", x=0.73, y=1.13, xref="paper", yref="paper",
        xanchor="right", yanchor="bottom", showarrow=False, font=dict(size=13),
    )
    fig.add_annotation(
        text="Distribution medians:", x=0.28, y=1.13, xref="paper", yref="paper",
        xanchor="left", yanchor="bottom", showarrow=False, font=dict(size=13),
    )
    # No in-figure "Units:" annotation here -- that control is a Dash-level
    # dcc.RadioItems + clientside_callback now (see the module docstring and
    # the note above), not an in-figure dropdown, so its label lives in
    # app.layout instead of as a Plotly annotation.
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Precomputed rate cache -- the "Years:" slider is the only control that can't
# be a pure client-side Plotly updatemenu (the KDEs/points/IMBIE slope all
# depend on year_start/year_end), so every slider move reruns
# plot_interactive_rate_comparison() on the server. Its costliest inner step
# is _sim_rows()'s per-simulation rate lookup (scipy.stats.linregress through
# pandas .groupby()/.loc slicing) -- profiling showed that's roughly a third
# of a single build's ~0.39s. Precomputing _sim_rows()'s output for all 351
# valid (year_start, year_end) combos up front via the *same* pandas/scipy
# path took ~61s (too slow to redo on every process start); _fast_slope
# below reproduces the identical rate using a closed-form OLS slope over raw
# numpy arrays instead, cutting that to ~1-2s, verified to match _sim_rows()
# exactly. The callback then becomes a dict lookup, not a recompute.
# ─────────────────────────────────────────────────────────────────────────────

def _fast_slope(x, y):
    """Closed-form OLS slope, numerically equivalent to
    scipy.stats.linregress(x, y).slope but without linregress's per-call
    nan-policy wrapper or the pandas-Series overhead of feeding it g["Year"]/
    g["..."] slices -- see the module-level comment above. Only used for
    bulk-precomputing _sim_rows()-equivalent rates; _sim_rows() itself (used
    for any single, non-precomputed call) still goes through
    imbie_mass_loss_slope()/scipy, matching the notebook's implementation."""
    n = len(x)
    if n < 2:
        return None
    xd = x - x.mean()
    ssxx = np.dot(xd, xd)
    if ssxx == 0:
        return None
    yd = y - y.mean()
    return np.dot(xd, yd) / ssxx


def _valid_year_ranges(year_min, year_max, min_span):
    """Every (lo, hi) the "Years:" RangeSlider can land on (post min-span
    self-correction): integers with year_min <= lo, hi <= year_max, and
    hi - lo >= min_span."""
    for lo in range(year_min, year_max - min_span + 1):
        for hi in range(lo + min_span, year_max + 1):
            yield (lo, hi)


def _precompute_rows_cache(simulated, year_min, year_max, min_span):
    """Builds {(year_start, year_end): {"AIS": [...], "GIS": [...]}}, with
    each inner list in exactly the row-dict shape _sim_rows() returns, for
    every valid (year_start, year_end). Extracts each (ice_sheet, group,
    model, exp) combo's NaN-dropped Year/mass-change arrays and its
    year-independent metadata (ice model, sliding law, GCM, scenario, base
    hover string) once, then reuses them across all 351 year windows --
    only the boolean mask and _fast_slope call vary per window."""
    combo_arrays = {}
    combo_meta = {}
    for ice_sheet in ["AIS", "GIS"]:
        df = simulated[simulated["IS"] == ice_sheet]
        for (group, model, exp), g in df.groupby(["Group", "Model", "Exp"]):
            mask = g["Cumulative ice sheet mass change (Gt)"].notna()
            key = (ice_sheet, group, model, exp)
            combo_arrays[key] = (
                g.loc[mask, "Year"].to_numpy(),
                g.loc[mask, "Cumulative ice sheet mass change (Gt)"].to_numpy(),
            )
            ism_m = get_ism_meta(group, model)
            exp_m = get_exp_meta(ice_sheet, exp)
            combo_meta[key] = {
                "initialization": ism_m.get("initialization", "See paper"),
                "ice_model": ism_m.get("ice_model", "See paper"),
                "sliding_law": ism_m.get("sliding_law", "See paper"),
                "scenario": exp_m.get("scenario", "Unknown"),
                "climate_model": exp_m.get("climate_model", "Unknown"),
                "base_hover": _build_hover(group, model, exp, ice_sheet),
            }

    cache = {}
    for lo, hi in _valid_year_ranges(year_min, year_max, min_span):
        per_ice_sheet = {"AIS": [], "GIS": []}
        for (ice_sheet, group, model, exp), (years, vals) in combo_arrays.items():
            m = (years >= lo) & (years <= hi)
            if m.sum() < 2:
                continue
            slope = _fast_slope(years[m], vals[m])
            if slope is None or not np.isfinite(slope):
                continue
            meta = combo_meta[(ice_sheet, group, model, exp)]
            per_ice_sheet[ice_sheet].append({
                "rate": slope, "group": group, "model": model, "exp": exp,
                "initialization": meta["initialization"], "ice_model": meta["ice_model"],
                "sliding_law": meta["sliding_law"], "scenario": meta["scenario"],
                "climate_model": meta["climate_model"],
                "hover": meta["base_hover"] + f"<br>Rate: {slope:.0f} Gt/yr",
                "hover_sle": meta["base_hover"] + f"<br>Rate: {slope * gt2mmSLE:.2f} mm/yr",
            })
        cache[(lo, hi)] = per_ice_sheet
    return cache


# ─────────────────────────────────────────────────────────────────────────────
# Build the figure once at startup (data doesn't change at runtime -- no need
# to recompute per-request), then serve it as a static Dash app.
# ─────────────────────────────────────────────────────────────────────────────

imbie_ais = _load_imbie2023("antarctica")
imbie_gis = _load_imbie2023("greenland")
imbie_ais["IS"] = "AIS"
imbie_gis["IS"] = "GIS"
imbie = pd.concat([imbie_ais, imbie_gis])

ismip6_ais = load_ismip6_ais()
ismip6_ais["IS"] = "AIS"
ismip6_gis = load_ismip6_gis()
ismip6_gis["IS"] = "GIS"
ismip6 = pd.concat([ismip6_ais, ismip6_gis])

YEAR_MIN, YEAR_MAX, MIN_YEAR_SPAN = 1990, 2020, 5
YEAR_DEFAULT = [2010, 2020]

TITLE_TEXT = "Compare observed and simulated rates of ice sheet change"

# Both are independent of year_start/year_end (colors are assigned from the
# full population; rates are precomputed for every valid window) so this
# happens once at startup, not per slider move -- see the precomputed rate
# cache comment above for why _fast_slope replaces _sim_rows()'s own
# pandas/scipy path here specifically. DIM_COLOR_MAPS is built from the
# COMPLETE EXTRA_SOURCES list (not whatever the "Simulation studies:"
# checklist currently has checked), so a category's color never shifts
# depending on which papers happen to be checked right now.
DIM_COLOR_MAPS = _dim_color_maps(ismip6, extra_sources=EXTRA_SOURCES)
ROWS_CACHE = _precompute_rows_cache(ismip6, YEAR_MIN, YEAR_MAX, MIN_YEAR_SPAN)

FIG = plot_interactive_rate_comparison(
    simulated=ismip6, observed=imbie, year_start=YEAR_DEFAULT[0], year_end=YEAR_DEFAULT[1],
    show_title=False, show_subtitle=False, extra_sources=EXTRA_SOURCES,
    precomputed_dim_color_maps=DIM_COLOR_MAPS, precomputed_rows=ROWS_CACHE[tuple(YEAR_DEFAULT)],
)

app = dash.Dash(__name__)
app.title = "ISMIP6 Rate of Ice Sheet Mass Change"
app.layout = html.Div(
    [
        # The figure's own title/subtitle are turned off above (show_title=
        # False, show_subtitle=False) and rendered here instead as plain
        # Dash HTML -- the title above the "Years:" slider, and the
        # explanatory text in a column alongside the graph rather than
        # crammed into the figure's top margin.
        html.Div(
            [
                html.H2(
                    TITLE_TEXT,
                    style={
                        "color": "#2a3f5f", "fontFamily": "Arial, sans-serif",
                        "fontWeight": "normal", "fontSize": "26px",
                        "margin": "0",
                    },
                ),
                # Small "app is updating" indicator -- shown whenever the
                # server is recomputing the graph (Years: slider drag or a
                # Simulation studies: checkbox change), via target_components
                # rather than wrapping the graph itself, so nothing dims or
                # blocks interaction while it spins: the slider/checklist/
                # graph all stay fully usable, this is purely informational.
                dcc.Loading(
                    target_components={"rate-comparison-graph": "figure"},
                    custom_spinner=html.Div(className="custom-slow-spinner"),
                    display="auto",
                    # Keeps the spinner visible at least this long once shown,
                    # so a fast callback (e.g. a Years: drag that lands back
                    # on an already-cached window) doesn't flash it on and
                    # off too quickly to actually notice.
                    delay_hide=400,
                ),
            ],
            style={
                "display": "flex", "alignItems": "center", "gap": "14px",
                "margin": "24px 0 0 40px",
            },
        ),
        html.Div(
            [
                html.Div(
                    "Select the time span over which to calculate the average rate of "
                    "change. Adjustments of this time span can take ~10 sec. to load.",
                    style={
                        "color": "#555555", "fontFamily": "Arial, sans-serif",
                        "fontSize": "13px", "lineHeight": "1.4", "maxWidth": "300px",
                        "marginRight": "24px",
                    },
                ),
                html.Label(
                    "Years:",
                    style={"fontWeight": "bold", "marginRight": "16px", "whiteSpace": "nowrap"},
                ),
                html.Div(
                    dcc.RangeSlider(
                        id="year-range-slider",
                        min=YEAR_MIN, max=YEAR_MAX, step=1,
                        value=YEAR_DEFAULT,
                        allowCross=False,
                        # Just the two end years -- at 420px wide, every-5-years marks
                        # (7 labels) overlapped/got clipped.
                        marks={YEAR_MIN: str(YEAR_MIN), YEAR_MAX: str(YEAR_MAX)},
                        # Only while dragging/hovering a handle, not persistently --
                        # the boundary-year marks are enough the rest of the time.
                        tooltip={"placement": "bottom"},
                    ),
                    # ~40% of the slider's previous width (it used to fill a flex:1
                    # slot spanning most of a 70%-wide row, roughly 1000-1100px at
                    # typical desktop widths). Fixed px, not %, since this div's
                    # immediate parent has no explicit width of its own for a
                    # percentage to resolve against.
                    style={"width": "420px", "minWidth": "220px"},
                ),
            ],
            style={
                "display": "flex", "alignItems": "center", "justifyContent": "center",
                "width": "100%", "margin": "20px auto 0 auto",
            },
        ),
        html.Div(
            [
                html.Label(
                    "Simulation studies:",
                    style={"fontWeight": "bold", "marginRight": "16px", "whiteSpace": "nowrap"},
                ),
                dcc.Checklist(
                    id="data-sources-checklist",
                    # Every extra_sources paper plus both core ISMIP6 panels,
                    # so unchecking Seroussi/Goelzer removes that whole panel's
                    # ISMIP6 overlay the same way unchecking a paper removes
                    # its own -- see _update_figure, which turns each checked
                    # label into either an EXTRA_SOURCES filter or an empty
                    # precomputed_rows["AIS"/"GIS"] list.
                    options=[{"label": f" {label}", "value": label} for label in DATA_SOURCE_LABELS],
                    value=list(DATA_SOURCE_LABELS),
                    inline=True,
                    style={"display": "flex", "flexWrap": "wrap", "gap": "4px 20px"},
                    inputStyle={"marginRight": "4px"},
                ),
            ],
            style={
                "display": "flex", "alignItems": "center", "justifyContent": "center",
                "width": "100%", "margin": "12px auto 0 auto",
                "fontFamily": "Arial, sans-serif", "fontSize": "14px", "color": "#2a3f5f",
            },
        ),
        html.Div(
            [
                html.Label(
                    "Units:",
                    style={"fontWeight": "bold", "marginRight": "16px", "whiteSpace": "nowrap"},
                ),
                dcc.RadioItems(
                    id="units-radio",
                    # A Dash-level control (not an in-figure Plotly updatemenu,
                    # unlike Group by/Distribution medians) -- see the
                    # clientside_callback below plot_interactive_rate_comparison
                    # for why: it lets the Mass<->Sea level rise conversion run
                    # as plain JS math on data already in the browser, instead
                    # of the server having to ship a full second copy of every
                    # trace's x/text/hovertemplate.
                    options=[
                        {"label": " Mass change", "value": "mass"},
                        {"label": " Sea level rise", "value": "sle"},
                    ],
                    value="mass",
                    inline=True,
                    style={"display": "flex", "gap": "4px 20px"},
                    inputStyle={"marginRight": "4px"},
                ),
                # Required Output target for the clientside_callback below --
                # it performs its work as a Plotly.restyle/relayout side
                # effect on the graph div directly, so this never needs to
                # hold anything meaningful, just satisfy Dash's requirement
                # that every callback have an Output.
                html.Div(id="units-clientside-dummy", style={"display": "none"}),
            ],
            style={
                "display": "flex", "alignItems": "center", "justifyContent": "center",
                "width": "100%", "margin": "12px auto 0 auto",
                "fontFamily": "Arial, sans-serif", "fontSize": "14px", "color": "#2a3f5f",
            },
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.P(
                            "This app allows for comparison between ISMIP6 simulations of "
                            "recent ice sheet change (Goelzer et al., 2020; Seroussi et al., "
                            "2020) and observations of ice sheet change (Otosaka et al., "
                            "2023). The goal of this app is to facilitate exploration of how "
                            "different modeling decisions affect simulated mass change."
                        ),
                        html.P(
                            "PDFs show the distribution of simulated mass changes over the "
                            "selected time span. PDFs are calculated collectively for all "
                            "simulations, or (optionally) grouped by a variety of different "
                            "simulation characteristics."
                        ),
                        html.P(
                            "Beneath the PDFs, each dot represents a single simulation, "
                            "plotted at its simulated rate of mass change."
                        ),
                        html.P(
                            "Hover over each dot to learn more about its characteristics, or "
                            "click on legend elements to hide or show plot components."
                        ),
                    ],
                    style={
                        "width": "15%", "minWidth": "220px", "padding": "24px 16px 0 40px",
                        "boxSizing": "border-box",
                        "color": "#555555", "fontFamily": "Arial, sans-serif",
                        "fontSize": "14px", "lineHeight": "1.5",
                    },
                ),
                dcc.Graph(
                    id="rate-comparison-graph",
                    figure=FIG,
                    # minWidth guards the in-figure "Distribution medians:"/"Group by:"
                    # dropdown row -- those controls have roughly fixed pixel footprints
                    # regardless of the figure's paper-coordinate width, so below ~1000px
                    # they start to visually collide. Letting the page scroll horizontally
                    # on narrow windows beats a broken control row.
                    style={"flex": "1", "minWidth": "1050px", "height": "85vh"},
                    config={"responsive": True, "displaylogo": False},
                ),
            ],
            style={"display": "flex", "alignItems": "flex-start"},
        ),
    ],
    style={"margin": 0, "padding": 0},
)


@app.callback(
    Output("rate-comparison-graph", "figure"),
    Output("year-range-slider", "value"),
    Input("year-range-slider", "value"),
    Input("data-sources-checklist", "value"),
)
def _update_figure(year_range, checked_sources):
    """Rebuilds the figure for the selected (year_start, year_end) and the
    currently-checked "Simulation studies:" checkboxes, enforcing a minimum
    MIN_YEAR_SPAN-year window.

    Unlike the "Units:"/"Group by:"/"Distribution medians:" dropdowns (pure
    client-side Plotly updatemenus toggling pre-built traces), the KDEs,
    jitter, medians, and IMBIE slope all depend on year_start/year_end, so
    this still reruns plot_interactive_rate_comparison on the server on
    every change -- mirroring the notebook's ipywidgets.IntRangeSlider
    callback -- but the per-simulation rates it needs are a ROWS_CACHE
    lookup (precomputed for every valid window at startup) rather than a
    fresh scipy/pandas computation, so only the Plotly trace-building itself
    still happens per request.

    An unchecked checkbox drops that source before the figure is even built
    (the checklist is the ONLY on/off control for a source -- there is no
    longer an in-figure dropdown duplicating this), so an unchecked paper's
    traces and legend entry disappear entirely rather than just being
    hidden. Each of EXTRA_SOURCES' own labels filters that list;
    ISMIP6_AIS_LABEL/ISMIP6_GIS_LABEL each gate one whole
    precomputed_rows["AIS"/"GIS"] entry -- since that's a plain list lookup
    (already computed in ROWS_CACHE at startup), swapping it for [] is free
    and skips that panel's entire ISMIP6 contribution inside
    plot_interactive_rate_comparison (its merged population is simply
    whatever extra_sources rows remain, same as an ISMIP6-only panel with
    no extra_sources checked at all).

    If the span is too narrow, the window is widened (growing from its
    center, then clamped to [YEAR_MIN, YEAR_MAX]) and the corrected value is
    written back to the slider via the second Output, which re-triggers this
    same callback once more with a valid span -- so an invalid, too-narrow
    window is never rendered.
    """
    lo, hi = year_range
    if hi - lo < MIN_YEAR_SPAN:
        center = (lo + hi) / 2
        lo, hi = center - MIN_YEAR_SPAN / 2, center + MIN_YEAR_SPAN / 2
        lo = max(YEAR_MIN, lo)
        hi = min(YEAR_MAX, hi)
        if hi - lo < MIN_YEAR_SPAN:
            lo = max(YEAR_MIN, hi - MIN_YEAR_SPAN)
            hi = min(YEAR_MAX, lo + MIN_YEAR_SPAN)
        lo, hi = round(lo), round(hi)

    checked = set(checked_sources or [])
    active_extra_sources = [src for src in EXTRA_SOURCES if src["label"] in checked]
    rows_for_window = ROWS_CACHE[(lo, hi)]
    active_rows = {
        "AIS": rows_for_window["AIS"] if ISMIP6_AIS_LABEL in checked else [],
        "GIS": rows_for_window["GIS"] if ISMIP6_GIS_LABEL in checked else [],
    }

    fig = plot_interactive_rate_comparison(
        simulated=ismip6, observed=imbie, year_start=lo, year_end=hi,
        show_title=False, show_subtitle=False, extra_sources=active_extra_sources,
        precomputed_dim_color_maps=DIM_COLOR_MAPS, precomputed_rows=active_rows,
    )
    return fig, [lo, hi]


# ─────────────────────────────────────────────────────────────────────────────
# "Units:" -- a clientside_callback (plain JS, runs entirely in the browser,
# no request back to this server) instead of a third native Plotly
# updatemenu. A native updatemenu can only apply a pre-baked data patch, so
# giving it a "Sea level rise" option would mean the server precomputing and
# shipping a full second copy of every trace's x/text/hovertemplate on every
# request -- previously ~10 MB of a ~27 MB payload (see the module docstring)
# for nothing but a *2.759 unit conversion. This does that multiplication in
# JS instead, always deriving from the `figure` Input (the FRESH, mass-Gt
# figure plot_interactive_rate_comparison() just built -- it never builds an
# SLE version itself), so it's correct both when the units radio changes and
# when the figure itself is rebuilt (Years:/Simulation studies: change)
# while "Sea level rise" is already selected -- both are listed as Inputs so
# either one re-applies the current unit choice.
#
# `text`/`customdata` on each point trace already carry BOTH hover-string
# variants (baked in server-side by _sim_rows, see plot_interactive_rate_
# comparison) -- this only ever swaps which one is assigned to `text`, it
# never computes new hover text itself.
app.clientside_callback(
    """
    function(units_value, figure) {
        if (!figure || !figure.data) {
            return "";
        }
        window.requestAnimationFrame(function() {
            var gd = document.getElementById("rate-comparison-graph");
            if (!gd || !gd.data) { return; }
            var GT2MMSLE = 1.0 / 362.5;  // Gt -> mm sea-level-equivalent, matches gt2mmSLE in app.py
            var factor = (units_value === "sle") ? GT2MMSLE : 1.0;

            var traceIdx = [], newX = [], newText = [], newHovertemplate = [];
            for (var i = 0; i < figure.data.length; i++) {
                var tr = figure.data[i];
                if (!tr.x) { continue; }
                traceIdx.push(i);
                newX.push(tr.x.map(function(v) { return v * factor; }));
                if (tr.customdata) {
                    newText.push(units_value === "sle" ? tr.customdata : tr.text);
                } else {
                    newText.push(tr.text || null);
                }
                var ht = tr.hovertemplate;
                if (ht && units_value === "sle") {
                    ht = ht.replace("%{x:.0f} Gt/yr", "%{x:.2f} mm/yr");
                }
                newHovertemplate.push(ht || null);
            }
            Plotly.restyle(gd, {x: newX, text: newText, hovertemplate: newHovertemplate}, traceIdx);

            var xTitle = (units_value === "sle")
                ? "Contribution to sea level rise (mm/yr)"
                : "Rate of ice sheet mass change (Gt/yr)";
            var relayoutUpdate = {"xaxis.title.text": xTitle, "xaxis2.title.text": xTitle};
            var L = figure.layout || {};
            if (L.xaxis && L.xaxis.range) {
                relayoutUpdate["xaxis.range"] = L.xaxis.range.map(function(v) { return v * factor; });
            }
            if (L.xaxis2 && L.xaxis2.range) {
                relayoutUpdate["xaxis2.range"] = L.xaxis2.range.map(function(v) { return v * factor; });
            }
            var annotations = L.annotations || [];
            for (var j = 0; j < annotations.length; j++) {
                var ann = annotations[j];
                if (ann.text && ann.text.indexOf("Observed mass change") !== -1) {
                    relayoutUpdate["annotations[" + j + "].text"] = (units_value === "sle")
                        ? "<b>Observed sea level contribution</b>" : "<b>Observed mass change</b>";
                    relayoutUpdate["annotations[" + j + "].x"] = ann.x * factor;
                }
            }
            Plotly.relayout(gd, relayoutUpdate);
        });
        return "";
    }
    """,
    Output("units-clientside-dummy", "children"),
    Input("units-radio", "value"),
    Input("rate-comparison-graph", "figure"),
)


server = app.server  # WSGI entry point, e.g. `gunicorn app:server`

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
