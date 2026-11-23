"""
Simulation metadata: one shared source of truth for the Dash app (data.py)
and the exploration notebook (exploring_fits/analyze_slr_predictions_interactive.ipynb).

Everything is read from the CSV tables in dash_app/metadata/ (comment lines
start with "#"):

  ismip6_ais_models.csv               reviewed: Seroussi et al. (2020) Table 3,
                                      Appendix C, plus the sliding-law literature check
  ismip6_ais_experiments.csv          reviewed: Seroussi et al. (2020) Table 1
  unreviewed_models.csv               ISMIP6-GIS and the non-ISMIP6 studies, carried
                                      over unchanged until each paper is reviewed
  ismip6_gis_experiments_unreviewed.csv

The review and every decision behind the class columns are documented in
metadata_review/01_Seroussi2020_ISMIP6_AIS.md (in the ismip6_check_261123
review folder).

Model metadata is keyed by (ice_sheet, group, model): every model
configuration/version is independent, and AIS and GIS rows can never borrow
from each other. Lookups are exact -- an unknown key raises KeyError rather
than falling back to a similar-looking group (the old substring fallback
silently mislabeled runs).
"""

import os

import pandas as pd

META_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metadata")

NOT_REVIEWED = "Not yet reviewed"
UNREVIEWED_SUFFIX = " (unreviewed)"

# Categorical fields a run carries, used for grouping and the bias ANOVA.
# initialization / sliding_law / stress_balance / resolution hold the CLASS
# (decided scheme); the verbatim paper codes and details live in other fields.
CLASS_FIELDS = ["ice_model", "initialization", "sliding_law", "stress_balance", "resolution"]
EXP_CLASS_FIELDS = ["ocean_forcing", "ocean_sensitivity"]

# Display order of each class vocabulary (anything else sorts after, alphabetically).
CATEGORY_ORDER = {
    "initialization": ["Spin-up / equilibrium", "Hybrid (thickness-nudged)", "Data assimilation (velocity inversion)"],
    "sliding_law": ["Linear", "Power law (Weertman)", "Budd", "Coulomb-limited"],
    "stress_balance": ["Shallow", "Higher-order"],
    "resolution": ["≤5 km", ">5–10 km", ">10–20 km", ">20 km"],
    "ocean_forcing": ["Open", "Standard"],
    "ocean_sensitivity": ["Low", "Medium", "High", "PIGL"],
}


def _read(name, **kw):
    return pd.read_csv(os.path.join(META_DIR, name), comment="#", keep_default_na=False, **kw)


def resolution_bin(res_min_km):
    """Bin of the finest grid resolution; an edge belongs to the lower bin."""
    r = float(res_min_km)
    return "≤5 km" if r <= 5 else ">5–10 km" if r <= 10 else ">10–20 km" if r <= 20 else ">20 km"


def _ais_models():
    df = _read("ismip6_ais_models.csv", dtype=str)
    out = {}
    for row in df.to_dict(orient="records"):
        assert row["resolution"] == resolution_bin(row["res_min_km"]), row
        row["review_status"] = "reviewed (Seroussi et al. 2020)"
        out[("AIS", row["group"], row["model"])] = row
    return out


def _unreviewed_models():
    df = _read("unreviewed_models.csv", dtype=str)
    out = {}
    for row in df.to_dict(orient="records"):
        na = row["review_status"] == "not applicable"
        if na:
            # Emulators/assessments: no ice sheet model physics to classify.
            fill = row["legacy_initialization"]
            row.update(initialization=fill, sliding_law=row["legacy_sliding_law"],
                       stress_balance=fill, resolution=fill)
        else:
            row.update(initialization=row["legacy_initialization"] + UNREVIEWED_SUFFIX,
                       sliding_law=row["legacy_sliding_law"] + UNREVIEWED_SUFFIX,
                       stress_balance=NOT_REVIEWED, resolution=NOT_REVIEWED)
        # The Edwards (2021) AIS rows are two studies sharing one model key;
        # they agree, so keep the first.
        out.setdefault((row["ice_sheet"], row["group"], row["model"]), row)
    return out


MODEL_META = {**_unreviewed_models(), **_ais_models()}


def model_meta(ice_sheet, group, model):
    """Metadata dict for one model configuration. Exact match only."""
    try:
        return MODEL_META[(ice_sheet, group, model)]
    except KeyError:
        raise KeyError(f"No metadata row for ice_sheet={ice_sheet!r}, group={group!r}, model={model!r}; "
                       "add it to dash_app/metadata/ (no fallback by design).") from None


def ismip6_experiments():
    """{(ice_sheet, exp): dict} for the ISMIP6 experiments."""
    ais = _read("ismip6_ais_experiments.csv", dtype=str)
    out = {}
    for row in ais.to_dict(orient="records"):
        row["ice_shelf_fracture"] = row["ice_shelf_fracture"] == "True"
        row["protocol"] = row["ocean_forcing"]  # kept for older callers
        out[("AIS", row["exp"])] = row
    gis = _read("ismip6_gis_experiments_unreviewed.csv", dtype=str)
    for row in gis.to_dict(orient="records"):
        row.update(ocean_forcing=NOT_REVIEWED, ice_shelf_fracture=False,
                   legacy_ocean_sensitivity=row.pop("ocean_sensitivity"), ocean_sensitivity=NOT_REVIEWED)
        out[("GIS", row["exp"])] = row
    return out


def external_exp_fields(model_row):
    """Experiment-level class fields for a non-ISMIP6 study, until reviewed."""
    fill = model_row["stress_balance"] if model_row["review_status"] == "not applicable" else NOT_REVIEWED
    return {"ocean_forcing": fill, "ocean_sensitivity": fill, "ice_shelf_fracture": False}


def category_order(dim, cats):
    """`cats` in the field's display order, unknown categories alphabetically after."""
    order = CATEGORY_ORDER.get(dim, [])
    known = [c for c in order if c in cats]
    return known + sorted(c for c in cats if c not in order)
