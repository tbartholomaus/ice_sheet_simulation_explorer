"""
Standalone Dash app serving the "Interactive Figure - Rate Comparison" from
exploring_fits/analyze_slr_predictions_interactive.ipynb.

The figure is pure Plotly. Its "Units:", "Distribution medians:", and
"Group by:" dropdowns are all native Plotly updatemenus toggling/restyling
pre-built traces, so they run entirely client-side with no callbacks back
to this server. (An earlier version tried moving "Units:" to a Dash-level
dcc.RadioItems + clientside_callback specifically to avoid precomputing a
second Sea-level-rise copy of every trace's x/text/hovertemplate
server-side -- ~10 MB of extra payload. That never actually worked in a
real browser despite checking out in every way testable without one, so it
reverted to the plain in-figure dropdown below, matching the notebook.) The
"Years:" range slider and the "Simulation studies:" checklist are different
from all of the above -- the KDEs, jitter, medians, and IMBIE slope all
depend on the selected year window and which sources are included, so
those genuinely need a real Dash callback that reruns
plot_interactive_rate_comparison() on the server.

That server-side rebuild is the app's main cost: building weighted KDEs
across 7 "Group by:" dimensions x 2 panels, potentially over ~1700+ pooled
ISMIP6 + extra_sources simulations, plus JSON-serializing the ~23 MB result
(back up from ~13 MB after reverting "Units:" -- see above). On a
resource-constrained deploy target this can exceed a WSGI server's
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

import collections
import math
import os
import time

# Pin every BLAS/OpenMP-based library this app touches (numpy, scipy,
# pandas' numexpr backend) to a single thread each, BEFORE they're
# imported -- these libraries size their internal thread pools once, the
# first time they're used, by reading these env vars (or falling back to
# "one thread per detected CPU core" if unset). A container's CPU CORE
# COUNT (what /proc/cpuinfo reports) and its actual CPU QUOTA (a thin
# fraction of a core, on Render's free tier) are different numbers --
# spawning threads sized to the former onto the latter causes thread
# contention/context-switching overhead that can make things far slower
# than the raw throttling alone would (observed: ~100x, not the ~10x a
# straightforward CPU-quota cut would predict). None of this app's own
# code benefits from BLAS/OpenMP parallelism -- every gaussian_kde/
# linregress call here operates on at most a few hundred points, well
# below where multi-threading would pay for its own overhead -- so pinning
# to 1 has no downside even on an unconstrained machine.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import psutil
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
    {"label": "Aschwanden 2019", "df": aschwanden2022_gis, "color": "#756bb1"},
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
# kind from Rahlves2025/Coulon2024/Aschwanden2019. Sourced from
# PUBLICATION_LABEL (plot_interactive_rate_comparison's own "Group by
# publication" category names) so the checklist and the figure never drift
# apart on what to call each panel's ISMIP6 population.
ISMIP6_AIS_LABEL = PUBLICATION_LABEL["AIS"]
ISMIP6_GIS_LABEL = PUBLICATION_LABEL["GIS"]
# ISMIP6 first (it's the app's core dataset, per its own title/intro text),
# then the three extra_sources papers.
DATA_SOURCE_LABELS = [ISMIP6_AIS_LABEL, ISMIP6_GIS_LABEL] + [src["label"] for src in EXTRA_SOURCES]
# Only the two ISMIP6 panels checked by default -- the extra_sources papers
# are opt-in, not on-by-default, so a first-time visitor sees the app's core
# ISMIP6-vs-IMBIE comparison before discovering the added papers.
DATA_SOURCE_DEFAULT_CHECKED = [ISMIP6_AIS_LABEL, ISMIP6_GIS_LABEL]


# ─────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

RCP_COLOR = {
    "RCP2.6":  "#003466",
    "RCP4.5":  "#b8860b",
    "RCP8.5":  "#990002",
    "SSP1-2.6": "#1a7f5e",
    "SSP2-4.5": "#8a6d3a",
    "SSP5-8.5": "#8b1a00",
    "Control":  "#555555",
    "Unknown":  "#888888",
}
IMBIE_COLOR = "#08519c"

gt2cmSLE = 1.0 / 362.5 / 10.0  # Gt -> cm sea-level-equivalent (magnitude only, see the sign flip below)
# Negative: mass change and sea-level rise are physically opposite in sign
# (losing ice -- negative Gt/yr -- raises sea level -- positive mm/yr -- and
# vice versa), so the "Sea level rise" unit view has to flip sign, not just
# rescale. The "Units:" toggle uses mm, not cm -- typical rates are well
# under 1 cm/yr.
gt2mmSLE = -gt2cmSLE * 10

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


def _rate_kde_raw(values, x_range, n=100, weights=None):
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
    this figure uses this instead of the `opacity` kwarg.

    Every color constant in this file (RCP_COLOR/INIT_COLOR/extra_sources'
    "color"/px.colors.qualitative.Dark24) is a hex string or literal "gray",
    so this deliberately doesn't pull in matplotlib (unlike the notebook's
    own _rgba, which uses colors.to_rgb and so accepts any named CSS color)
    -- that's real import/memory weight for a case that never happens today.
    It does handle 3-digit hex shorthand and passthrough rgb(...)/rgba(...)
    strings, and raises a clear error for anything else instead of a
    confusing IndexError/ValueError deep in string slicing, in case a
    future color constant is added in one of those other forms."""
    if color in ("gray", "grey"):
        r, g, b = 128, 128, 128
    elif color.startswith("#"):
        c = color.lstrip("#")
        if len(c) == 3:
            c = "".join(ch * 2 for ch in c)
        if len(c) != 6:
            raise ValueError(f"_rgba: unsupported hex color {color!r} (expected #rgb or #rrggbb)")
        r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
    elif color.startswith(("rgb(", "rgba(")):
        r, g, b = (int(v) for v in color[color.index("(") + 1:color.index(")")].split(",")[:3])
    else:
        raise ValueError(
            f"_rgba: unsupported color format {color!r} -- only hex (#rgb/#rrggbb), "
            f"'gray'/'grey', or rgb(...)/rgba(...) strings are supported"
        )
    return f"rgba({r}, {g}, {b}, {alpha})"


def _nice_sle_ticks(gt_range, factor, target_ticks=6):
    """Picks ~target_ticks "nice" (1/2/5 x 10^n) round sea-level-rise
    (mm/yr) values spanning gt_range (a (min, max) tuple in Gt/yr, i.e. a
    panel's x_range_by_panel entry), and returns (gt_positions, sle_labels):
    gt_positions are where those nice SLE values actually fall on the
    (never-moved) Gt/yr axis (sle_value / factor), and sle_labels are their
    display text. Used so the "Units:" dropdown's "Sea level rise" option
    can relabel the x-axis with sign-flipped, rescaled tick text WITHOUT
    moving a single point -- see that dropdown's construction below."""
    lo, hi = sorted(v * factor for v in gt_range)
    span = hi - lo
    if span <= 0:
        return [], []
    raw_step = span / target_ticks
    magnitude = 10 ** math.floor(math.log10(raw_step))
    step = 10 * magnitude
    for m in (1, 2, 5, 10):
        if raw_step <= m * magnitude:
            step = m * magnitude
            break
    start = math.ceil(lo / step) * step
    sle_values = []
    v = start
    while v <= hi + step * 1e-6:
        sle_values.append(round(v, 10))
        v += step
    gt_positions = [v / factor for v in sle_values]
    labels = [f"{v:g}" for v in sle_values]
    return gt_positions, labels


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


# "Group by publication" is a special-cased dimension (see the weighting note
# in the merged_df.groupby(dim) loop below), not one of GROUP_DIMENSIONS's
# generic entries -- inserted right after "All simulations" (GROUP_DIMENSIONS[0])
# rather than appended at the end, so it reads as the second-most-fundamental
# way to view the data, ahead of the more granular modeling-choice dimensions.
GROUP_BY_LABELS = (
    [GROUP_DIMENSIONS[0][0], "Group by publication"]
    + [label for label, _dim, _fixed in GROUP_DIMENSIONS[1:]]
)


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
    # legendgroups for the fixed "All simulations" traces (kde_all/all_pts/
    # median:all) already given a legend entry. NOT hardcoded to "whichever
    # trace belongs to k==1 (AIS)" -- a panel can have zero data (e.g. its
    # ISMIP6 checkbox is unchecked and no extra source covers that ice
    # sheet), skipping these traces entirely for that k, which silently
    # dropped the legend entry for BOTH panels when k==1 was assumed to
    # always be the one bearing it. Tracked dynamically instead, same
    # principle as legend_shown above.
    fixed_legend_shown = set()
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
            # Use the lazily-memoized cache when this source has one (see
            # its construction above, near ROWS_CACHE) -- indexing with
            # [key], not .get(key), so a first-time window actually
            # triggers LazyRowsCache's __missing__ compute-and-memoize
            # rather than silently falling through every time. Falls back
            # to _sim_rows() itself only when there's no cache at all
            # (e.g. the notebook's own extra_sources, which never sets
            # "rows_cache").
            rows_cache = src.get("rows_cache")
            src_rows = rows_cache[(year_start, year_end)][ice_sheet] if rows_cache is not None else _sim_rows(src["df"], ice_sheet, year_start, year_end)
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
            show_kde_all = "kde_all" not in fixed_legend_shown
            fixed_legend_shown.add("kde_all")
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=xs, y=rate_kde_y0 + dens, mode="lines", line=dict(color="gray", width=1),
                fill="tonexty", fillcolor=_rgba("gray", 0.5), hoverinfo="skip",
                name="PDF of all simulations", legendgroup="kde_all", showlegend=show_kde_all,
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
            show_all_pts = "all_pts" not in fixed_legend_shown
            fixed_legend_shown.add("all_pts")
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=merged_df["rate"], y=merged_df["jitter"], mode="markers",
                marker=dict(
                    color="gray", line_width=0,
                    size=merged_df["publication"].map(size_by_pub).tolist(),
                    opacity=merged_df["publication"].map(opacity_by_pub).tolist(),
                ),
                name="Simulation", legendgroup="all_pts", showlegend=show_all_pts, visible=True,
                text=merged_df["hover"], customdata=merged_df["hover_sle"], hovertemplate="%{text}<extra></extra>",
            ), row=k, col=1)
            all_only_idx.append(idx)

            show_median_all = "median:all" not in fixed_legend_shown
            fixed_legend_shown.add("median:all")
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=[merged_df["rate"].median()], y=[rate_median_y], mode="markers",
                marker=dict(symbol="triangle-down", size=11, color="gray", line=dict(width=1, color="black")),
                opacity=1, name="All simulations median", legendgroup="median:all", showlegend=show_median_all,
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
                        # "Group by publication" compares whole sources/studies
                        # to each other, not simulation counts -- ISMIP6 (summed
                        # across ~15-20 institutions) has vastly more total
                        # row_weight than any one paper's own (already
                        # de-weighted-to-typical_model_n) equivalent, so
                        # weighting by row_weight share here would make every
                        # paper's curve nearly invisible next to ISMIP6's.
                        # Every OTHER dimension keeps proportional-by-weight
                        # scaling (a category with more simulations legitimately
                        # gets a bigger curve there), but publication instead
                        # gives every category (ISMIP6 counts as one) equal
                        # area, so it's each source's distribution SHAPE being
                        # compared, not how many institutions/ensemble members
                        # went into computing it.
                        weight = (1.0 / n_cats) if dim == "publication" else (cat_weight / total_weight)
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

    # "Units:" dropdown -- switches between showing rates as ice-sheet mass
    # change (Gt/yr) and as their sea-level-rise equivalent (mm/yr). Every
    # trace's x data stays in Gt/yr NO MATTER which unit is selected -- an
    # explicit user request that points never move or get replotted when
    # this is toggled. "Sea level rise" is achieved purely by RELABELING
    # the same physical x-axis positions: explicit tickvals/ticktext (see
    # _nice_sle_ticks) showing each position's sign-flipped, rescaled
    # sea-level-equivalent value, instead of Plotly's own auto-ticking
    # (which can only label an axis by the data's actual values, not some
    # other unit's equivalent). Not moving x also means neither button
    # needs to ship its own copy of every trace's x array (previously the
    # single largest contributor to this dropdown's payload cost).
    mass_text_all = [tr.text for tr in fig.data]
    # Point traces stash their alt-unit hover string in customdata (set
    # alongside "text" at trace-creation); every other trace type either
    # has no text at all or doesn't depend on units (KDE fills, IMBIE),
    # so it falls back to its own (untouched) text.
    sle_text_all = [
        tr.customdata if tr.customdata is not None else tr.text for tr in fig.data
    ]
    mass_hovertemplate_all = [tr.hovertemplate for tr in fig.data]
    # Only the median markers' hovertemplate hardcodes a unit (the point
    # traces' hovertemplate is just "%{text}...", already unit-agnostic).
    sle_hovertemplate_all = [
        ht.replace("%{x:.0f} Gt/yr", "%{x:.2f} mm/yr") if ht is not None else None
        for ht in mass_hovertemplate_all
    ]

    sle_ticks_1 = _nice_sle_ticks(x_range_by_panel[1], gt2mmSLE) if 1 in x_range_by_panel else ([], [])
    sle_ticks_2 = _nice_sle_ticks(x_range_by_panel[2], gt2mmSLE) if 2 in x_range_by_panel else ([], [])

    # The "Observed ..." arrow annotation is pinned to the IMBIE slope in
    # data coordinates (annotation.x) -- left untouched here (same "nothing
    # moves" rule as the traces above), only its text changes.
    mass_annotation_updates, sle_annotation_updates = {}, {}
    for _idx, _slope in observed_annotation_x.items():
        mass_annotation_updates[f"annotations[{_idx}].text"] = "<b>Observed mass change</b>"
        sle_annotation_updates[f"annotations[{_idx}].text"] = "<b>Observed sea level contribution</b>"

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
            dict(
                type="dropdown", direction="down",
                x=0.06, y=1.13, xanchor="left", yanchor="bottom",
                active=0,
                buttons=[
                    dict(label="Mass change", method="update", args=[
                        {"text": mass_text_all, "hovertemplate": mass_hovertemplate_all},
                        {"xaxis.title.text": "Rate of ice sheet mass change (Gt/yr)",
                         "xaxis2.title.text": "Rate of ice sheet mass change (Gt/yr)",
                         "xaxis.tickmode": "auto", "xaxis2.tickmode": "auto",
                         "xaxis.tickvals": None, "xaxis2.tickvals": None,
                         "xaxis.ticktext": None, "xaxis2.ticktext": None,
                         **mass_annotation_updates},
                    ]),
                    dict(label="Sea level rise", method="update", args=[
                        {"text": sle_text_all, "hovertemplate": sle_hovertemplate_all},
                        {"xaxis.title.text": "Contribution to sea level rise (mm/yr)",
                         "xaxis2.title.text": "Contribution to sea level rise (mm/yr)",
                         "xaxis.tickmode": "array", "xaxis.tickvals": sle_ticks_1[0], "xaxis.ticktext": sle_ticks_1[1],
                         "xaxis2.tickmode": "array", "xaxis2.tickvals": sle_ticks_2[0], "xaxis2.ticktext": sle_ticks_2[1],
                         **sle_annotation_updates},
                    ]),
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
    fig.add_annotation(
        text="Units:", x=0.0, y=1.13, xref="paper", yref="paper",
        xanchor="left", yanchor="bottom", showarrow=False, font=dict(size=13),
    )
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


def _build_combo_arrays(simulated):
    """The expensive, one-time (not per-window) prep step: for every
    (ice_sheet, group, model, exp) combo, its NaN-dropped Year/mass-change
    arrays and its year-independent metadata (ice model, sliding law, GCM,
    scenario, base hover string) -- the groupby/get_ism_meta/get_exp_meta/
    _build_hover work, paid once regardless of how many year-windows are
    ever actually requested. Cost scales with the number of (group, model,
    exp) combos in `simulated`, NOT with the number of valid year windows."""
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
    return combo_arrays, combo_meta


def _rows_for_window(combo_arrays, combo_meta, lo, hi):
    """The cheap per-window step (a boolean mask + _fast_slope call over
    already-prepared arrays) -- what used to run for all 351 valid windows
    eagerly at startup. Called lazily instead, see LazyRowsCache below."""
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
    return per_ice_sheet


class LazyRowsCache(dict):
    """{(year_start, year_end): {"AIS": [...], "GIS": [...]}}, computed and
    memoized lazily on first access (via dict's own __missing__ hook, so
    ordinary `cache[(lo, hi)]` indexing everywhere else needs no changes)
    instead of eagerly building all 351 valid windows at import time.

    This replaces an earlier version (_precompute_rows_cache) that built
    every window upfront for ISMIP6, which was then naively extended to
    each of the 3 extra_sources dataframes too -- for Aschwanden2019's
    ~500-member ensemble alone that meant ~350,000 pre-built row-dicts held
    in memory before a single request ever arrived, nearly doubling this
    app's import-time RSS (measured locally: ~513 MB -> ~960 MB) and
    causing Render's free-tier instance to OOM during startup, before
    gunicorn could even bind a port. A real deployed session only ever
    touches a handful of the 351 possible windows (wherever the Years:
    slider actually gets dragged to), so eagerly building all of them was
    pure waste. The one-time _build_combo_arrays() setup this still does at
    startup is cheap (proportional to the number of (group, model, exp)
    combos, not the number of windows) -- only the expensive "build every
    window" step is now deferred and memoized per-window instead."""

    def __init__(self, simulated):
        super().__init__()
        self._combo_arrays, self._combo_meta = _build_combo_arrays(simulated)

    def __missing__(self, key):
        lo, hi = key
        value = _rows_for_window(self._combo_arrays, self._combo_meta, lo, hi)
        self[key] = value
        return value


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

YEAR_MIN, YEAR_MAX, MIN_YEAR_SPAN = 2000, 2020, 5
YEAR_DEFAULT = [2015, 2020]

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
ROWS_CACHE = LazyRowsCache(ismip6)
# Same lazily-memoized per-window-rate cache ISMIP6 gets above, one per
# extra_sources entry -- _sim_rows() on a paper's own dataframe (e.g.
# Aschwanden2019's ~500 LHS ensemble members) was previously recomputed via
# the slow scipy.stats.linregress path on every single request regardless
# of whether the year window had changed, which was real, avoidable
# per-request CPU/memory pressure. LazyRowsCache is already generic over
# any (Group, Model, Exp, IS)-shaped dataframe, so this is a direct reuse,
# not a new code path. plot_interactive_rate_comparison() below checks
# each source dict for this key and falls back to the original _sim_rows()
# call when it's absent -- keeping the notebook (which passes plain
# extra_sources dicts with no such cache) unaffected.
for _src in EXTRA_SOURCES:
    _src["rows_cache"] = LazyRowsCache(_src["df"])

# Matches DATA_SOURCE_DEFAULT_CHECKED (only the two ISMIP6 panels) so the
# page's very first render is consistent with what the "Simulation studies:"
# checklist shows checked -- otherwise a user would see all 5 sources on
# load despite only 2 checkboxes being ticked, until their first interaction.
_default_checked = set(DATA_SOURCE_DEFAULT_CHECKED)
FIG = plot_interactive_rate_comparison(
    simulated=ismip6, observed=imbie, year_start=YEAR_DEFAULT[0], year_end=YEAR_DEFAULT[1],
    show_title=False, show_subtitle=False,
    extra_sources=[src for src in EXTRA_SOURCES if src["label"] in _default_checked],
    precomputed_dim_color_maps=DIM_COLOR_MAPS,
    precomputed_rows={
        "AIS": ROWS_CACHE[tuple(YEAR_DEFAULT)]["AIS"] if ISMIP6_AIS_LABEL in _default_checked else [],
        "GIS": ROWS_CACHE[tuple(YEAR_DEFAULT)]["GIS"] if ISMIP6_GIS_LABEL in _default_checked else [],
    },
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
                # Simulation studies: checkbox change). Uses the classic
                # dcc.Loading(children=...) pattern (wrapping a hidden div
                # that _update_figure also targets as an Output), not the
                # newer target_components prop -- target_components never
                # visibly showed anything in practice, so this switches to
                # the long-established, heavily-used mechanism instead. The
                # wrapped div itself is empty/invisible either way; only the
                # spinner (positioned here, next to the title, not over the
                # graph) is meant to be seen, so nothing dims or blocks
                # interaction while it spins.
                dcc.Loading(
                    children=html.Div(id="loading-trigger", style={"display": "none"}),
                    custom_spinner=html.Div(className="custom-slow-spinner"),
                    display="auto",
                    # Keeps the spinner visible at least this long once shown.
                    # Both Years: and Simulation studies: target the SAME
                    # _update_figure callback -- interacting with one while
                    # the other's request is still in flight (or landing back
                    # in an already-cached window quickly) can end one
                    # loading event and start another in close succession;
                    # without a large enough bridge here, that shows up as
                    # the spinner flickering off between them instead of
                    # staying continuously visible for the whole time the
                    # user is actively interacting with either control.
                    delay_hide=1000,
                ),
            ],
            style={
                "display": "flex", "alignItems": "center", "gap": "14px",
                "margin": "24px 0 0 40px",
            },
        ),
        # Three text boxes around the Years:/Simulation studies: controls: a
        # "this can take ~10 sec" note on the far left, vertically centered
        # against BOTH rows combined (not just one), and a "Select the time
        # window..." label + the "Simulation studies:" label left-aligned
        # with each other directly above the controls each one describes --
        # achieved by giving both labels the same fixed width and stacking
        # the two rows in their own flex column, so the outer row's
        # alignItems:center centers the note against that whole column's
        # height, and both inner rows share one consistent left edge rather
        # than each independently self-centering (which is what the
        # previous version -- two separate top-level, individually
        # justifyContent:center'd rows -- could not express: two rows
        # centered as a group don't share a left edge unless their content
        # widths happen to match by coincidence).
        html.Div(
            [
                html.Div(
                    "Adjustments to the averaging window or the simulation "
                    "studies can take ~10 seconds to load.",
                    style={
                        "color": "#555555", "fontFamily": "Arial, sans-serif",
                        "fontSize": "13px", "lineHeight": "1.4", "maxWidth": "220px",
                        "marginRight": "32px",
                    },
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    "Select the time window over which to "
                                    "calculate the average rate of ice sheet "
                                    "change.",
                                    style={
                                        "width": "220px", "minWidth": "220px",
                                        "color": "#555555", "fontFamily": "Arial, sans-serif",
                                        "fontSize": "13px", "lineHeight": "1.4",
                                        "marginRight": "16px",
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
                            style={"display": "flex", "alignItems": "center"},
                        ),
                        html.Div(
                            [
                                html.Label(
                                    "Simulation studies:",
                                    style={
                                        "width": "220px", "minWidth": "220px",
                                        "fontWeight": "bold", "marginRight": "16px",
                                    },
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
                                    value=list(DATA_SOURCE_DEFAULT_CHECKED),
                                    inline=True,
                                    style={"display": "flex", "flexWrap": "wrap", "gap": "4px 20px"},
                                    inputStyle={"marginRight": "4px"},
                                ),
                            ],
                            style={
                                "display": "flex", "alignItems": "center",
                                "marginTop": "12px",
                                "fontFamily": "Arial, sans-serif", "fontSize": "14px", "color": "#2a3f5f",
                            },
                        ),
                    ],
                    style={"display": "flex", "flexDirection": "column"},
                ),
            ],
            style={
                "display": "flex", "alignItems": "center", "justifyContent": "center",
                "width": "100%", "margin": "20px auto 0 auto",
            },
        ),
        # "Units:" is a native in-figure Plotly dropdown now (see the
        # "Mass change"/"Sea level rise" buttons inside
        # plot_interactive_rate_comparison), not a Dash-level control -- an
        # earlier RadioItems + clientside_callback version lived here, but
        # never worked in practice (see that function's docstring), so this
        # was reverted to the same in-figure mechanism the notebook uses.
        html.Div(
            [
                html.Div(
                    [
                        html.P(
                            "This app allows comparison between simulations of recent ice "
                            "sheet change, including from ISMIP6, and observations of ice "
                            "sheet change (Otosaka et al., 2023). The goal of this app is to "
                            "facilitate exploration of how different modeling decisions "
                            "affect simulated mass change."
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
                            "Hover over each dot to learn more about its characteristics, "
                            "click on legend elements to hide or show plot components, or "
                            "use the tools at top right of the plots to zoom or pan through "
                            "the plots."
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
        # Small live CPU/memory readout -- position:fixed takes it entirely
        # out of the normal document flow (pinned to the viewport corner),
        # so adding it here can't shift or resize anything else in the
        # layout above, regardless of where in this children list it's
        # appended. Polls this process's own usage every 2s via
        # dcc.Interval -- added after chasing a real Render OOM (see
        # LazyRowsCache's docstring and render.yaml's --max-requests
        # comment) to make memory/CPU behavior visible during local
        # testing without needing to shell in or read server logs.
        dcc.Interval(id="perf-monitor-interval", interval=2000, n_intervals=0),
        html.Div(
            id="perf-monitor",
            style={
                "position": "fixed", "bottom": "8px", "right": "8px",
                "backgroundColor": "rgba(255,255,255,0.9)",
                "border": "1px solid #ccc", "borderRadius": "4px",
                "padding": "3px 8px", "fontFamily": "monospace",
                "fontSize": "11px", "color": "#555555", "zIndex": 1000,
                "pointerEvents": "none",
            },
        ),
    ],
    style={"margin": 0, "padding": 0},
)


# A single, reused Process handle -- cpu_percent(interval=None) reports
# usage SINCE THAT SAME OBJECT'S LAST CALL (0.0 on an object's first-ever
# call, with no prior reading to diff against), so calling
# psutil.Process().cpu_percent(...) fresh inside the callback below would
# silently read 0.0% forever. The throwaway priming call establishes a
# baseline immediately so the very first tick already reports something
# meaningful instead of 0.0.
_PERF_PROC = psutil.Process()
_PERF_PROC.cpu_percent(interval=None)

# (monotonic timestamp, cpu_pct, mem_mb) readings from the last
# _PERF_WINDOW_SECONDS, for the "Max of last 5 min." line -- time.monotonic()
# rather than time.time(), since this is only ever compared against other
# readings from this same process's uptime, and monotonic is immune to
# wall-clock adjustments (NTP sync, DST, etc.) that could otherwise make the
# trailing-window cutoff jump backward or forward.
_PERF_HISTORY = collections.deque()
_PERF_WINDOW_SECONDS = 5 * 60


@app.callback(
    Output("perf-monitor", "children"),
    Input("perf-monitor-interval", "n_intervals"),
)
def _update_perf_monitor(_n_intervals):
    """Live CPU%/memory readout for this worker process specifically (not
    the whole machine), plus the running max of each over the trailing 5
    minutes -- cpu_percent(interval=None) is non-blocking and reports usage
    since _PERF_PROC's last call, which is exactly what a per-tick poll on
    a persistent handle wants; interval=<float> would instead block this
    callback for that long on every single tick."""
    now = time.monotonic()
    cpu_pct = _PERF_PROC.cpu_percent(interval=None)
    mem_mb = _PERF_PROC.memory_info().rss / 1e6

    _PERF_HISTORY.append((now, cpu_pct, mem_mb))
    cutoff = now - _PERF_WINDOW_SECONDS
    while _PERF_HISTORY and _PERF_HISTORY[0][0] < cutoff:
        _PERF_HISTORY.popleft()
    max_cpu = max(r[1] for r in _PERF_HISTORY)
    max_mem = max(r[2] for r in _PERF_HISTORY)

    return [
        html.Div(f"Live: CPU {cpu_pct:5.1f}%  |  Mem {mem_mb:6.0f} MB"),
        html.Div(f"Max of last 5 min.: CPU {max_cpu:5.1f}%  |  Mem {max_mem:6.0f} MB"),
    ]


@app.callback(
    Output("rate-comparison-graph", "figure"),
    Output("year-range-slider", "value"),
    Output("loading-trigger", "children"),
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

    The third Output (loading-trigger.children) carries no information of
    its own -- it exists purely so the dcc.Loading wrapping that div (see
    app.layout) shows its spinner for exactly the duration of this callback,
    the standard dcc.Loading(children=...) pattern.
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
    # The third return value is meaningless on its own -- its only purpose is
    # being this callback's Output("loading-trigger", "children"), which is
    # what makes the dcc.Loading wrapping that div show its spinner for the
    # duration of this callback (see app.layout).
    return fig, [lo, hi], ""


server = app.server  # WSGI entry point, e.g. `gunicorn app:server`

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
