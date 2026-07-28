"""
Standalone Dash app serving the "Interactive Figure - Rate Comparison" from
exploring_fits/analyze_slr_predictions_interactive.ipynb.

The figure is pure Plotly. Its "Units:", "Distribution medians:", and
"Group by:" dropdowns are native Plotly updatemenus toggling/restyling
pre-built traces, so they run entirely client-side with no callbacks back to
this server. The "Years:" range slider is different -- the KDEs, jitter,
medians, and IMBIE slope all depend on the selected year window, so it's a
real Dash callback that reruns plot_interactive_rate_comparison() on the
server for each new (year_start, year_end).

Run locally:
    pip install -r requirements.txt
    python app.py
    # -> http://127.0.0.1:8050

Deploy (any WSGI host, e.g. gunicorn behind nginx on your own server):
    gunicorn app:server -b 0.0.0.0:8000

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
    ("ULB",      "FETISH1"):   {"ice_model": "f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
    ("ULB",      "FETISH2"):   {"ice_model": "f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
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
    ("NEMO",     "fETISh"):    {"ice_model": "f.ETISh",   "sliding_law": "Weertman / Coulomb",       "initialization": "Data assimilation"},
    ("PSU",      "PSU3D1"):    {"ice_model": "PSU-ISM",   "sliding_law": "Coulomb / Weertman",       "initialization": "Spin-up"},
    ("PSU",      "PSU3D2"):    {"ice_model": "PSU-ISM",   "sliding_law": "Coulomb / Weertman",       "initialization": "Spin-up"},
    ("UTAS",     "ElmerIce"):  {"ice_model": "Elmer/Ice", "sliding_law": "Regularized Coulomb",      "initialization": "Data assimilation"},
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


def _rate_kde_raw(values, x_range, n=200):
    """Raw KDE curve (integrates to 1 over its own support) -- NOT renormalized
    to its own peak. Callers weight and rescale this against a shared reference
    so peak heights stay comparable across subgroups of very different spread.

    Falls back to an all-zero (flat) curve if the covariance is singular --
    e.g. a narrow year window can leave a subgroup with near-but-not-exactly-
    zero variance, which passes a `std() > 0` guard but still isn't enough
    for gaussian_kde's Cholesky step."""
    xs = np.linspace(x_range[0], x_range[1], n)
    try:
        dens = scipy.stats.gaussian_kde(values)(xs)
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


def _dim_color_maps(simulated):
    """Category -> color, per grouping dimension, built from the FULL
    simulated dataframe (not windowed by year) so a category's color stays
    consistent across both panels and every year-range selection. Factored
    out of plot_interactive_rate_comparison so a caller (e.g. the Dash app)
    can compute this once and pass it back in via the `precomputed_dim_color_maps`
    kwarg, instead of paying its groupby cost on every request -- it's the
    same regardless of year_start/year_end."""
    dims = [d for _, d, _ in GROUP_DIMENSIONS if d is not None]
    dim_color_maps = {}
    if simulated is None:
        return dim_color_maps
    dim_cats = {d: set() for d in dims}
    for (ice_sheet_, group, model, exp), _ in simulated.groupby(["IS", "Group", "Model", "Exp"]):
        ism_m = get_ism_meta(group, model)
        exp_m = get_exp_meta(ice_sheet_, exp)
        row_meta = {**ism_m, **exp_m}
        for d in dims:
            dim_cats[d].add(row_meta.get(d, "See paper"))
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


def plot_interactive_rate_comparison(
    simulated=None, observed=None, year_start=2000, year_end=2020,
    title="Compare observed and simulated rates of ice sheet change", seed=42,
    show_title=True, show_subtitle=True,
    precomputed_dim_color_maps=None, precomputed_rows=None,
):
    """
    Interactive dual-panel (AIS / GIS) raincloud comparison of the mean
    ice-sheet mass-change rate, matching the static plot_rate_comparison
    figure's default look, plus two independent controls:

    - "Group by:" dropdown -- recolor by initialization method, ice sheet
      model, sliding law, climate scenario, or GCM (see GROUP_DIMENSIONS).
      "All simulations": gray jittered points + gray KDE, exactly like the
      static IS_rate_comparison.pdf figure. Any "Group by ..." option: colored
      jittered points and colored KDEs (one per category of the chosen
      dimension, drawn on top of the ever-present gray KDE). Each category's
      raw (un-renormalized) KDE is weighted by its share of the total
      simulation count and then scaled using the SAME reference (the full
      population's peak density) as the gray curve -- not its own peak. This
      keeps the AREA under each colored curve proportional to its share of
      simulations, so a widely-scattered subgroup reads as low and wide
      rather than being inflated to look as prominent as any other category,
      and a tightly-clustered subgroup can legitimately show as a sharp,
      narrow spike even with relatively few points. Summing every category's
      curve approximately reconstructs the gray curve.
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
    `simulated` as usual when omitted.
    """
    rng = np.random.default_rng(seed)
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=False,
        subplot_titles=["Antarctic Ice Sheet (AIS)", "Greenland Ice Sheet (GIS)"],
        vertical_spacing=0.20,
    )

    # Global color maps (shared between AIS/GIS panels) per grouping dimension,
    # built from the full simulated dataframe so a category's color is
    # consistent across both panels.
    dim_color_maps = (
        precomputed_dim_color_maps if precomputed_dim_color_maps is not None
        else _dim_color_maps(simulated)
    )

    dim_only_idx = {label: [] for label, dim, _ in GROUP_DIMENSIONS if dim is not None}
    all_only_idx = []
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

        if len(sim_df):
            n_total = len(sim_df)
            pad = max(0.15 * (sim_df["rate"].max() - sim_df["rate"].min()), 1.0)
            x_range = (sim_df["rate"].min() - pad, sim_df["rate"].max() + pad)
            x_range_by_panel[k] = x_range

            # Fixed per-simulation jitter, computed once and reused by both the
            # "all simulations" trace and every grouped-by-dimension trace below,
            # so a given simulation's marker stays at the same y position no
            # matter which dropdown option is selected.
            sim_df["jitter"] = rate_strip_y0 + rng.uniform(
                -rate_strip_jitter, rate_strip_jitter, size=len(sim_df)
            )

            # Gray "all simulations" KDE -- always visible, matches the static
            # figure. full_max is the shared reference every colored subgroup
            # below is scaled against (not each subgroup's own peak).
            xs, dens_full_raw = _rate_kde_raw(sim_df["rate"].values, x_range)
            full_max = dens_full_raw.max()
            if full_max <= 0:
                # Degenerate window (e.g. all simulations converged to the
                # same rate) -- avoid a division by zero; nothing meaningful
                # to scale against, so every curve stays flat.
                full_max = 1.0
            dens = dens_full_raw / full_max * rate_kde_height
            fig.add_trace(go.Scatter(
                x=xs, y=np.full_like(xs, rate_kde_y0), mode="lines", line=dict(width=0),
                hoverinfo="skip", showlegend=False,
            ), row=k, col=1)
            fig.add_trace(go.Scatter(
                x=xs, y=rate_kde_y0 + dens, mode="lines", line=dict(color="gray", width=1),
                fill="tonexty", fillcolor=_rgba("gray", 0.5), hoverinfo="skip",
                name="PDF of all simulations", legendgroup="kde_all", showlegend=(k == 1),
            ), row=k, col=1)

            # Colored KDE + colored jittered points + median marker per
            # category, for each grouping dimension, drawn on top of the gray
            # KDE.
            for label, dim, _fixed in GROUP_DIMENSIONS:
                if dim is None:
                    continue
                color_map = dim_color_maps[dim]
                # Semi-transparent fills compound when many categories overlap in the
                # same x-range (n stacked layers at opacity a look like 1-(1-a)^n, so
                # even a=0.5 reads as ~90% opaque by the 4th overlapping layer). Scale
                # each category's fill opacity down for dimensions with more
                # categories (e.g. ice sheet model, ~14) so the *stacked* result stays
                # translucent, while dimensions with few categories (e.g.
                # initialization, ~4) keep the full 0.5.
                n_cats = max(len(color_map), 1)
                kde_fill_opacity = min(0.5, 2.0 / n_cats)
                for cat, g in sim_df.groupby(dim):
                    color = color_map.get(cat, "#888888")
                    if len(g) >= 2 and g["rate"].std() > 0:
                        _, dens2_raw = _rate_kde_raw(g["rate"].values, x_range)
                        weight = len(g) / n_total
                        dens2 = (dens2_raw * weight) / full_max * rate_kde_height
                        idx = len(fig.data)
                        fig.add_trace(go.Scatter(
                            x=xs, y=np.full_like(xs, rate_kde_y0), mode="lines", line=dict(width=0),
                            hoverinfo="skip", showlegend=False, visible=False,
                        ), row=k, col=1)
                        dim_only_idx[label].append(idx)
                        idx = len(fig.data)
                        fig.add_trace(go.Scatter(
                            x=xs, y=rate_kde_y0 + dens2, mode="lines", line=dict(color=color, width=1),
                            fill="tonexty", fillcolor=_rgba(color, kde_fill_opacity), hoverinfo="skip",
                            name=f"{cat} KDE", legendgroup=f"{dim}:{cat}", showlegend=False, visible=False,
                        ), row=k, col=1)
                        dim_only_idx[label].append(idx)

                    idx = len(fig.data)
                    show_this = (dim, cat) not in legend_shown
                    if show_this:
                        legend_shown.add((dim, cat))
                    fig.add_trace(go.Scatter(
                        x=g["rate"], y=g["jitter"], mode="markers",
                        marker=dict(color=color, size=4, opacity=0.85, line_width=0),
                        name=str(cat), legendgroup=f"{dim}:{cat}", showlegend=show_this, visible=False,
                        text=g["hover"], customdata=g["hover_sle"], hovertemplate="%{text}<extra></extra>",
                    ), row=k, col=1)
                    dim_only_idx[label].append(idx)

                    # Median marker: visibility tied to the grouping dropdown
                    # (like its KDE/points siblings, part of dim_only_idx), but
                    # opacity toggled independently by the "Medians:" dropdown.
                    idx = len(fig.data)
                    fig.add_trace(go.Scatter(
                        x=[g["rate"].median()], y=[rate_median_y], mode="markers",
                        marker=dict(symbol="triangle-down", size=11, color=color, line=dict(width=1, color="black")),
                        opacity=1, name=f"{cat} median", legendgroup=f"median:{dim}:{cat}", showlegend=False,
                        visible=False, hovertemplate=f"{cat} median<br>" + "Rate: %{x:.0f} Gt/yr<extra></extra>",
                    ), row=k, col=1)
                    dim_only_idx[label].append(idx)
                    median_trace_idx.append(idx)

            # Gray jittered points -- "All simulations" only.
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=sim_df["rate"], y=sim_df["jitter"], mode="markers",
                marker=dict(color="gray", size=4, opacity=0.6, line_width=0),
                name="Simulation", legendgroup="all_pts", showlegend=(k == 1), visible=True,
                text=sim_df["hover"], customdata=sim_df["hover_sle"], hovertemplate="%{text}<extra></extra>",
            ), row=k, col=1)
            all_only_idx.append(idx)

            # Median marker for "All simulations" (gray).
            idx = len(fig.data)
            fig.add_trace(go.Scatter(
                x=[sim_df["rate"].median()], y=[rate_median_y], mode="markers",
                marker=dict(symbol="triangle-down", size=11, color="gray", line=dict(width=1, color="black")),
                opacity=1, name="All simulations median", legendgroup="median:all", showlegend=(k == 1),
                visible=True, hovertemplate="All simulations median<br>" + "Rate: %{x:.0f} Gt/yr<extra></extra>",
            ), row=k, col=1)
            all_only_idx.append(idx)
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
    buttons = []
    for label, dim, _fixed in GROUP_DIMENSIONS:
        if dim is None:
            visible = [i not in sum(dim_only_idx.values(), []) for i in range(n_traces)]
        else:
            other_idx = sum((v for l, v in dim_only_idx.items() if l != label), []) + all_only_idx
            visible = [i not in other_idx for i in range(n_traces)]
        buttons.append(dict(label=label, method="update", args=[{"visible": visible}]))

    # "Units:" dropdown -- every trace's x-data represents a rate in Gt/yr;
    # the Sea level rise view is just that same data times gt2mmSLE. Rather
    # than rebuild the figure, read each trace's already-built x/text/
    # hovertemplate back out and restyle to the alt-unit version computed
    # from it, so the two controls above (Group by / Distribution medians)
    # keep working unmodified regardless of which unit is selected.
    mass_x_all = [tr.x for tr in fig.data]
    sle_x_all = [
        tuple(v * gt2mmSLE for v in tr.x) if tr.x is not None else None for tr in fig.data
    ]
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

    def _panel_range(row, factor=1.0):
        r = x_range_by_panel.get(row)
        return [r[0] * factor, r[1] * factor] if r is not None else None

    # The "Observed ..." arrow annotation(s) are pinned to the IMBIE slope in
    # data coordinates (annotation.x), so switching units has to move the
    # arrow itself, not just relabel it -- relayout can target a specific
    # annotation's properties via "annotations[i].<prop>".
    mass_annotation_updates, sle_annotation_updates = {}, {}
    for _idx, _slope in observed_annotation_x.items():
        mass_annotation_updates[f"annotations[{_idx}].text"] = "<b>Observed mass change</b>"
        mass_annotation_updates[f"annotations[{_idx}].x"] = _slope
        sle_annotation_updates[f"annotations[{_idx}].text"] = "<b>Observed sea level contribution</b>"
        sle_annotation_updates[f"annotations[{_idx}].x"] = _slope * gt2mmSLE

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
                        {"x": mass_x_all, "text": mass_text_all, "hovertemplate": mass_hovertemplate_all},
                        {"xaxis.title.text": "Rate of ice sheet mass change (Gt/yr)",
                         "xaxis2.title.text": "Rate of ice sheet mass change (Gt/yr)",
                         "xaxis.range": _panel_range(1), "xaxis2.range": _panel_range(2),
                         **mass_annotation_updates},
                    ]),
                    dict(label="Sea level rise", method="update", args=[
                        {"x": sle_x_all, "text": sle_text_all, "hovertemplate": sle_hovertemplate_all},
                        {"xaxis.title.text": "Contribution to sea level rise (mm/yr)",
                         "xaxis2.title.text": "Contribution to sea level rise (mm/yr)",
                         "xaxis.range": _panel_range(1, gt2mmSLE), "xaxis2.range": _panel_range(2, gt2mmSLE),
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
# pandas/scipy path here specifically.
DIM_COLOR_MAPS = _dim_color_maps(ismip6)
ROWS_CACHE = _precompute_rows_cache(ismip6, YEAR_MIN, YEAR_MAX, MIN_YEAR_SPAN)

FIG = plot_interactive_rate_comparison(
    simulated=ismip6, observed=imbie, year_start=YEAR_DEFAULT[0], year_end=YEAR_DEFAULT[1],
    show_title=False, show_subtitle=False,
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
        html.H2(
            TITLE_TEXT,
            style={
                "color": "#2a3f5f", "fontFamily": "Arial, sans-serif",
                "fontWeight": "normal", "fontSize": "26px",
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
                    # minWidth guards the "Units:"/"Distribution medians:"/"Group by:" dropdown
                    # row -- those controls have roughly fixed pixel footprints regardless of the
                    # figure's paper-coordinate width, so below ~1000px they start to visually
                    # collide. Letting the page scroll horizontally on narrow windows beats a
                    # broken control row.
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
)
def _update_figure(year_range):
    """Rebuilds the figure for the selected (year_start, year_end), enforcing
    a minimum MIN_YEAR_SPAN-year window.

    Unlike the "Units:"/"Group by:"/"Distribution medians:" dropdowns (pure
    client-side Plotly updatemenus toggling pre-built traces), the KDEs,
    jitter, medians, and IMBIE slope all depend on year_start/year_end, so
    this still reruns plot_interactive_rate_comparison on the server on
    every change -- mirroring the notebook's ipywidgets.IntRangeSlider
    callback -- but the per-simulation rates it needs are a ROWS_CACHE
    lookup (precomputed for every valid window at startup) rather than a
    fresh scipy/pandas computation, so only the Plotly trace-building itself
    still happens per request.

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

    fig = plot_interactive_rate_comparison(
        simulated=ismip6, observed=imbie, year_start=lo, year_end=hi,
        show_title=False, show_subtitle=False,
        precomputed_dim_color_maps=DIM_COLOR_MAPS, precomputed_rows=ROWS_CACHE[(lo, hi)],
    )
    return fig, [lo, hi]


server = app.server  # WSGI entry point, e.g. `gunicorn app:server`

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
