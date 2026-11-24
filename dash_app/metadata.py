"""
Simulation metadata: one shared source of truth for the Dash app (data.py)
and the exploration notebook (exploring_fits/analyze_slr_predictions_interactive.ipynb).

Everything lives in CSV files in dash_app/metadata/ (lines starting with "#"
are comments), so the encoding can be checked without reading any code:

  <Study>_models.csv        one row per model configuration (Group/Model). For
                            each characteristic, *_described says how the sources
                            describe it, *_source where, and *_encoded the value
                            the app uses.
  <Study>_experiments.csv   one row per experiment, transcribed from the study
                            (verbatim column names and values where the study
                            has an experiment table).
  experiment_columns.csv    which experiment column fills which app field, and
                            how verbatim values map to encoded ones.
  encoding_rules.csv        the allowed encoded values of each field, in display
                            order, and the rule that assigns each.

Model metadata is keyed by (ice_sheet, group, model): every model
configuration is independent, and AIS and GIS rows never borrow from each
other. Lookups are exact; an unknown key raises KeyError. Every file is
validated when this module is imported.
"""

import os

import pandas as pd

META_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metadata")

NOT_REVIEWED = "Not yet reviewed"
UNREVIEWED_SUFFIX = " (unreviewed)"

# The studies, in the order their files are read, with the ice sheet of
# single-ice-sheet studies (the others have an "Ice sheet" column).
STUDIES = {
    "Seroussi2020": "AIS", "Goelzer2020": "GIS", "DeConto2016": "AIS", "Aschwanden2019": "GIS",
    "DeConto2021": "AIS", "Coulon2024": "AIS", "Rahlves2025": "GIS", "Goelzer2025": "GIS",
    "Edwards2021": None, "IPCCAR6": None,
}

# Model characteristics, as run-table columns (the encoded value of each).
CLASS_FIELDS = ["ice_model", "initialization", "sliding_law", "stress_balance", "resolution"]
# Experiment characteristics the app groups by.
EXP_FIELDS = ["climate_model", "scenario", "ocean_forcing", "ocean_sensitivity", "ice_shelf_fracture"]


def _read(name):
    return pd.read_csv(os.path.join(META_DIR, name), comment="#", dtype=str, keep_default_na=False)


def resolution_bin(res_min_km):
    """Bin of the finest grid resolution; a bin edge belongs to the lower bin."""
    r = float(res_min_km)
    return "≤5 km" if r <= 5 else ">5–10 km" if r <= 10 else ">10–20 km" if r <= 20 else ">20 km"


# ── encoding rules ────────────────────────────────────────────────────────

RULES = _read("encoding_rules.csv")
VOCAB = {f: list(g["encoded_value"]) for f, g in RULES[~RULES["field"].str.startswith("(")].groupby("field", sort=False)}
CATEGORY_ORDER = VOCAB  # display order of each field's values


def _check_encoded(field, value, status, where):
    if field == "ice_model":
        ok = bool(value)
    elif status == "reviewed":
        ok = value in VOCAB[field]
    elif status == "not yet reviewed":
        ok = value == NOT_REVIEWED or value.endswith(UNREVIEWED_SUFFIX)
    else:  # not applicable: emulators and assessments
        ok = value.startswith("Not applicable")
    if not ok:
        raise ValueError(f"{where}: {field} = {value!r} isn't allowed for review_status {status!r} "
                         "(see encoding_rules.csv)")


# ── models ────────────────────────────────────────────────────────────────

def _load_models():
    out = {}
    for stem, ice in STUDIES.items():
        for row in _read(f"{stem}_models.csv").to_dict(orient="records"):
            row["study"] = stem
            row["ice_sheet"] = row.pop("Ice sheet", ice)
            where = f"{stem}_models.csv {row['Group']}/{row['Model']}"
            for f in CLASS_FIELDS:
                row[f] = row[f"{f}_encoded"]
                _check_encoded(f, row[f], row["review_status"], where)
            if row["review_status"] == "reviewed" and row["resolution"] != resolution_bin(row["res_min_km"]):
                raise ValueError(f"{where}: resolution_encoded doesn't match res_min_km")
            # Optional per-configuration ocean forcing (e.g. ISMIP6-GIS runs made in
            # the open framework); it overrides the experiment's value in run_meta.
            if row.get("ocean_forcing_encoded") and row["ocean_forcing_encoded"] not in VOCAB["ocean_forcing"]:
                raise ValueError(f"{where}: ocean_forcing_encoded = {row['ocean_forcing_encoded']!r} isn't in encoding_rules.csv")
            key = (row["ice_sheet"], row["Group"], row["Model"])
            if key in out:
                raise ValueError(f"{where}: duplicate model key {key}")
            out[key] = row
    return out


MODEL_META = _load_models()


def model_meta(ice_sheet, group, model):
    """Metadata of one model configuration. Exact match only."""
    try:
        return MODEL_META[(ice_sheet, group, model)]
    except KeyError:
        raise KeyError(f"No metadata row for ice_sheet={ice_sheet!r}, group={group!r}, model={model!r}; "
                       "add it to the study's file in dash_app/metadata/.") from None


# ── experiments ───────────────────────────────────────────────────────────

COLUMNS = _read("experiment_columns.csv")


def _value_map(text):
    return dict(pair.split("=", 1) for pair in (p.strip() for p in text.split(";")) if pair)


def _load_experiments():
    out = {}
    for stem, ice in STUDIES.items():
        spec = COLUMNS[COLUMNS["study"] == stem]
        to_field = dict(zip(spec["file_column"], spec["app_field"]))
        maps = {c: _value_map(v) for c, v in zip(spec["file_column"], spec["value_map"]) if v}
        key_col = next(c for c, f in to_field.items() if f == "(key)" and c != "Ice sheet")
        df = _read(f"{stem}_experiments.csv")
        undescribed = set(df.columns) - set(to_field)
        if undescribed:
            raise ValueError(f"{stem}_experiments.csv: columns {sorted(undescribed)} aren't in experiment_columns.csv")
        for row in df.to_dict(orient="records"):
            meta = {"study": stem, "extras": {}}
            for col, val in row.items():
                field = to_field[col]
                val = maps.get(col, {}).get(val, val)
                if field == "(hover)":
                    if val:
                        meta["extras"][col.replace("_", " ")] = val
                elif field != "(key)":
                    meta[field] = val
            if "ice_shelf_fracture" in meta:
                meta["ice_shelf_fracture"] = meta["ice_shelf_fracture"] == "True"
            for f in ("ocean_forcing", "ocean_sensitivity"):
                if f in meta and meta[f] not in VOCAB[f]:
                    raise ValueError(f"{stem}_experiments.csv {row[key_col]}: {f} = {meta[f]!r} "
                                     "isn't in encoding_rules.csv")
            key = (stem, row.get("Ice sheet", ice), row[key_col])
            if key in out:
                raise ValueError(f"{stem}_experiments.csv: duplicate experiment {key}")
            out[key] = meta
    return out


EXP_META = _load_experiments()


def run_meta(ice_sheet, group, model, exp):
    """(model metadata, experiment metadata) of one run. Experiment metadata
    is completed with the app fields its study doesn't provide: "Not yet
    reviewed", or "Not applicable (...)" for emulators and assessments."""
    m = model_meta(ice_sheet, group, model)
    try:
        e = EXP_META[(m["study"], ice_sheet, str(exp))]
    except KeyError:
        raise KeyError(f"No row for experiment {exp!r} ({ice_sheet}) in {m['study']}_experiments.csv") from None
    fill = m["stress_balance"] if m["review_status"] == "not applicable" else NOT_REVIEWED
    e = {"ocean_forcing": fill, "ocean_sensitivity": fill, "ice_shelf_fracture": False, **e}
    if m.get("ocean_forcing_encoded"):
        e["ocean_forcing"] = m["ocean_forcing_encoded"]
    return m, e


def category_order(dim, cats):
    """`cats` in the field's display order, any others alphabetically after."""
    order = CATEGORY_ORDER.get(dim, [])
    return [c for c in order if c in cats] + sorted(c for c in cats if c not in order)
