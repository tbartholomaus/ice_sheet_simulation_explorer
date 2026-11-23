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

import metadata as M

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


def _imbie_arrays(df):
    """A loaded IMBIE frame as plain arrays: a few hundred monthly values,
    small enough that analysis.py recomputes from them per request."""
    return {"year": df["Year"].to_numpy(dtype=float), "mass": df[MASS_COL].to_numpy(dtype=float),
            "sigma": df["Cumulative ice sheet mass change uncertainty (Gt)"].to_numpy(dtype=float)}


# IMBIE[obs][ice_sheet] -> {"year", "mass" (cumulative Gt), "sigma" (its 1-sigma, Gt)}
IMBIE = {
    "imbie3": {"AIS": _imbie_arrays(_load_imbie2026("antarctica")), "GIS": _imbie_arrays(_load_imbie2026("greenland"))},
    "imbie2": {"AIS": _imbie_arrays(_load_imbie2023("antarctica")), "GIS": _imbie_arrays(_load_imbie2023("greenland"))},
}


# ─────────────────────────────────────────────────────────────────────────────
# Simulation metadata lives in dash_app/metadata/ and is read through
# metadata.py -- one documented table per source, shared with the notebook,
# keyed exactly by (ice_sheet, group, model) for models and by
# (publication, experiment) for experiments. ISMIP6-AIS (Seroussi et al.
# 2020) has been reviewed against the paper and the wider literature;
# everything else is carried over unchanged and labelled "(unreviewed)" /
# "Not yet reviewed" until its paper is reviewed.
# ─────────────────────────────────────────────────────────────────────────────

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


# Source-specific experiment columns carried into the hover text. A study's
# own "ocean_sensitivity" (Rahlves 2025) is kept as legacy_ocean_sensitivity,
# so it can't be mistaken for the reviewed ISMIP6 field of the same name.
_EXTRA_EXP_COLS = ["ocean_sensitivity", "basal_melt_param", "retreat_percentile", "mici_params"]


def _exp_meta_from_df(df):
    cols = ["Exp", "climate_model", "scenario", "protocol"] + [c for c in _EXTRA_EXP_COLS if c in df.columns]
    unique = df[cols].drop_duplicates("Exp")
    out = {}
    for row in unique.to_dict(orient="records"):
        if "ocean_sensitivity" in row:
            row["legacy_ocean_sensitivity"] = row.pop("ocean_sensitivity")
        out[str(row.pop("Exp"))] = {k: (str(v) if v == v else "") for k, v in row.items()}
    return out


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
ar6_ais = _read_csv("external_sources_ar6_ais.csv.gz")
ar6_gis = _read_csv("external_sources_ar6_gis.csv.gz")


# Colors (2026-10-09) were chosen to be clearly distinct for normal colour
# vision (CIEDE2000 >= 20 between any two studies, including within each ice
# sheet's panel); Edwards Main/GIS and AR6 AIS/GIS keep one color per study,
# and Edwards Risk Averse is a darker relative of Edwards Main.
# Display order everywhere (sidebar, Study legends/boxes, README): ISMIP6
# first (SOURCE_LABELS below), then the other ice-sheet-model studies by
# publication date, then the Edwards et al. (2021) emulator studies and the
# IPCC AR6 assessment last, since they aren't ice sheet model runs themselves
# (user decisions, 2026-10-08 and 2026-10-09).
EXTRA_SOURCES = [
    {"label": "DeConto & Pollard 2016", "df": deconto2016_ais, "color": "#7f3b08"},
    {"label": "Aschwanden 2019", "df": aschwanden2022_gis, "color": "#5e3c99"},
    {"label": "DeConto 2021", "df": deconto2021_ais, "color": "#e6ab02"},
    {"label": "Coulon 2024", "df": coulon2024_ais, "color": "#117733"},
    {"label": "Rahlves 2025", "df": rahlves2025_gis, "color": "#ff7f0e"},  # TC 19, 1205 (Mar 2025)
    {"label": "Goelzer 2025 (PROTECT GIS)", "df": goelzer2025_gis, "color": "#0072b2"},  # TC 19, 6887 (later 2025)
    # Edwards et al. (2021) covers BOTH ice sheets, so its labels carry an
    # ice-sheet suffix to stay unambiguous. Same color for AIS Main and GIS --
    # they render on separate AIS/GIS panels, so there's no legend collision.
    {"label": "Edwards 2021 (AIS Main)", "df": edwards2021_ais, "color": "#e7298a"},
    {"label": "Edwards 2021 (AIS Risk Averse)", "df": edwards2021_ais_risk, "color": "#8e0152"},
    {"label": "Edwards 2021 (GIS)", "df": edwards2021_gis, "color": "#e7298a"},
    # IPCC AR6 WG1 (2021) assessed projections, also not ice sheet model runs:
    # last, after Edwards. Each "run" is one percentile path (1st-99th) of the
    # medium- or low-confidence distribution; see load_ar6_icesheet().
    {"label": "IPCC AR6 (AIS)", "df": ar6_ais, "color": "#17becf"},
    {"label": "IPCC AR6 (GIS)", "df": ar6_gis, "color": "#17becf"},
]

for _src, _is in [
    ("Rahlves 2025", "GIS"), ("Coulon 2024", "AIS"), ("Aschwanden 2019", "GIS"),
    ("Goelzer 2025 (PROTECT GIS)", "GIS"), ("Edwards 2021 (AIS Main)", "AIS"),
    ("Edwards 2021 (AIS Risk Averse)", "AIS"), ("Edwards 2021 (GIS)", "GIS"), ("DeConto & Pollard 2016", "AIS"),
    ("DeConto 2021", "AIS"), ("IPCC AR6 (AIS)", "AIS"), ("IPCC AR6 (GIS)", "GIS"),
]:
    next(s for s in EXTRA_SOURCES if s["label"] == _src)["df"]["IS"] = _is

PUBLICATION_LABEL = {"AIS": "Seroussi 2020 (ISMIP6 AIS)", "GIS": "Goelzer 2020 (ISMIP6 GIS)"}
ISMIP6_COLOR = "#636363"
SOURCE_LABELS = [PUBLICATION_LABEL["AIS"], PUBLICATION_LABEL["GIS"]] + [s["label"] for s in EXTRA_SOURCES]
SOURCE_DEFAULT_CHECKED = [PUBLICATION_LABEL["AIS"], PUBLICATION_LABEL["GIS"]]
# Studies that aren't runs of a physics-based ice sheet model: statistical
# emulators and assessments. They're listed after the physics-based models,
# under their own heading in the sidebar.
EMULATION_SOURCES = {"Edwards 2021 (AIS Main)", "Edwards 2021 (AIS Risk Averse)", "Edwards 2021 (GIS)",
                     "IPCC AR6 (AIS)", "IPCC AR6 (GIS)"}
assert SOURCE_LABELS[-len(EMULATION_SOURCES):] == [s for s in SOURCE_LABELS if s in EMULATION_SOURCES], \
    "emulation/synthesis studies must come last in EXTRA_SOURCES"
SOURCE_COLOR = {PUBLICATION_LABEL["AIS"]: ISMIP6_COLOR, PUBLICATION_LABEL["GIS"]: ISMIP6_COLOR,
                **{s["label"]: s["color"] for s in EXTRA_SOURCES}}

# Experiment metadata keyed by (publication, experiment id): ISMIP6 from the
# metadata tables, every other study from its own bundled CSV. Keying by
# study means two studies can never overwrite each other's experiment ids.
EXP_META = {(PUBLICATION_LABEL[ice], exp): v for (ice, exp), v in M.ismip6_experiments().items()}
for _s in EXTRA_SOURCES:
    for _exp, _v in _exp_meta_from_df(_s["df"]).items():
        EXP_META[(_s["label"], _exp)] = _v

# (key, label) -- key is a run-table column; None = no grouping. The model
# characteristics hold the classes decided in the metadata review
# (metadata.CATEGORY_ORDER lists each vocabulary).
GROUP_DIMENSIONS = [
    (None, "All simulations"),
    ("publication", "Study"),
    ("initialization", "Initialization"),
    ("ice_model", "Ice sheet model"),
    ("sliding_law", "Sliding law"),
    ("stress_balance", "Stress balance"),
    ("resolution", "Resolution (finest grid)"),
    ("ocean_forcing", "Ocean forcing"),
    ("ocean_sensitivity", "Ocean sensitivity"),
    ("scenario", "Climate scenario"),
    ("climate_model", "GCM"),
]
ANOVA_CHARACTERISTICS = ["ice_model", "sliding_law", "initialization", "stress_balance", "resolution",
                         "ocean_forcing", "ocean_sensitivity", "scenario", "climate_model"]
CATEGORY_ORDER = M.CATEGORY_ORDER
DIM_LABEL = {k: v for k, v in GROUP_DIMENSIONS if k}

# "Combine RCPs and SSPs": each RCP pooled with the SSP of the same nominal
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
    "initialization": {"Spin-up / equilibrium": "#003466", "Hybrid (thickness-nudged)": "#b8860b",
                       "Data assimilation (velocity inversion)": "#1a7f5e"},
    "scenario": {"RCP2.6": "#003466", "RCP4.5": "#b8860b", "RCP8.5": "#990002", "SSP1-2.6": "#1a7f5e",
                 "SSP2-4.5": "#8a6d3a", "SSP5-8.5": "#8b1a00", "SSP1-1.9": "#4daf4a", "SSP3-7.0": "#d95f02",
                 "NDC (current pledges)": "#7570b3", "Control": "#555555", "Unknown": "#888888"},
    "publication": SOURCE_COLOR,
}
_FIXED_COLORS[COMPOSITE_DIM] = {
    **{k: v for k, v in _FIXED_COLORS["scenario"].items() if k not in SCENARIO_COMPOSITE},
    **{c: _FIXED_COLORS["scenario"][ssp] for ssp, c in SCENARIO_COMPOSITE.items() if ssp.startswith("SSP")},
}


def _km(row):
    lo, hi = row.get("res_min_km", ""), row.get("res_max_km", "")
    if not lo:
        return ""
    return f" ({lo} km)" if lo == hi else f" ({lo}–{hi} km)"


def _build_hover(group, model, exp, ism_m, exp_m, publication):
    code = lambda k: f" [{ism_m[k]}]" if ism_m.get(k) else ""  # noqa: E731 -- verbatim paper code
    sliding = ism_m["sliding_law"]
    if ism_m.get("sliding_law_detail"):
        sliding += f" — {ism_m['sliding_law_detail']}, {ism_m['sliding_exponent']}"
    lines = [
        f"<b>{group} / {model}</b> ({publication})",
        f"Experiment: {exp}",
        f"Ice model: {ism_m['ice_model']}",
        f"Initialization: {ism_m['initialization']}{code('init_code')}",
        f"Sliding law: {sliding}",
        f"Stress balance: {ism_m['stress_balance']}{code('stress_balance_code')}",
        f"Resolution: {ism_m['resolution']}{_km(ism_m)}",
        f"Climate model: {exp_m['climate_model']}",
        f"Scenario: {exp_m['scenario']}",
        f"Ocean forcing: {exp_m['ocean_forcing']}",
        f"Ocean sensitivity: {exp_m['ocean_sensitivity']}",
    ]
    if publication not in PUBLICATION_LABEL.values():
        lines.append(f"Protocol: {exp_m['protocol']}")
    for key, label in (("basal_melt_param", "Basal melt param"), ("legacy_ocean_sensitivity", "Study's ocean sensitivity"),
                       ("retreat_percentile", "Retreat percentile"), ("mici_params", None)):
        if exp_m.get(key):
            lines.append(f"{label}: {exp_m[key]}" if label else exp_m[key])
    if exp_m["ice_shelf_fracture"]:
        lines.append("<b>Ice-shelf collapse imposed: excluded from the bias analysis</b>")
    if ism_m["review_status"] == "not yet reviewed":
        lines.append("<i>Metadata not yet reviewed against the source paper</i>")
    return "<br>".join(lines)


def _build_runs():
    """One row per simulation across every source, plus lists of each run's
    raw (Year, cumulative Gt) arrays and an (n_runs x YEAR_GRID) matrix on
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
            # Exact lookups only: a run without a metadata row is an error,
            # never a silent borrow from a similar-looking model.
            ism_m = M.model_meta(ice_sheet, group, model)
            exp_m = EXP_META[(publication, exp)]
            if publication not in PUBLICATION_LABEL.values():
                exp_m = {**exp_m, **M.external_exp_fields(ism_m)}
            meta.append({
                "ice_sheet": ice_sheet, "publication": publication,
                "group": group, "model": model, "exp": exp,
                **{k: ism_m[k] for k in M.CLASS_FIELDS},
                **{k: exp_m[k] for k in M.EXP_CLASS_FIELDS},
                "scenario": exp_m["scenario"],
                "climate_model": exp_m["climate_model"],
                "ice_shelf_fracture": bool(exp_m["ice_shelf_fracture"]),
                "review_status": ism_m["review_status"],
                "hover": _build_hover(group, model, exp, ism_m, exp_m, publication),
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
# Runs have different time steps and spans, so rather than one small array
# per run, every run's points are stacked end to end in run order, with
# RUN_ID[k] naming the run (RUNS row) point k belongs to. Lets
# analysis.run_rates fit all runs at once instead of looping over them.
RUN_ID = np.repeat(np.arange(len(RUN_YEARS)), [len(y) for y in RUN_YEARS])
RUN_YEARS, RUN_VALS = np.concatenate(RUN_YEARS), np.concatenate(RUN_VALS)
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
del edwards2021_ais, edwards2021_ais_risk, edwards2021_gis, deconto2016_ais, deconto2021_ais, ar6_ais, ar6_gis
import gc  # noqa: E402
gc.collect()
IS_ISMIP6 = RUNS["publication"].isin(PUBLICATION_LABEL.values()).to_numpy()
# Runs that enter the bias analysis (rate bias, low-bias shares, ANOVA).
# Imposed ice-shelf-collapse experiments (ISMIP6-AIS exp11/exp12) are left
# out there but still drawn everywhere else (decision 2026-11-23).
BIAS_ELIGIBLE = ~RUNS["ice_shelf_fracture"].to_numpy(dtype=bool)


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
